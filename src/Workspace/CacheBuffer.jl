"""
A write-back, read-through staging layer in front of the DuckDB cache.

The scan and processing controllers never touch the cache directly: they *deposit* interpreted and
processed item data here and *read it back* here. Deposited data is immediately readable from memory,
so a just-interpreted or just-processed item can be served before it is durable. A dedicated flusher
task drains the staged writes through the cache's single locked writer on a fixed cadence (sooner
under memory pressure), then evicts the now-durable entries so later reads fall through to DuckDB.

This unifies every per-item write — interpreted scan output, processed payloads, item statistics, and
item failures — onto one path, and decouples the scan from cache-write latency: depositing never
blocks on the database. Because reads consult the staged entries first, flush cadence affects only
durability and memory, never interactive latency.

The cache keeps one serialized writer (see [`CacheDB`](@ref)); the buffer is what makes that writer
cheap to feed — it amortizes many deposits into one transaction per flush.
"""

"""Durability/eviction cadence: stage in memory, commit a transaction at most this far apart."""
const CACHE_FLUSH_INTERVAL_NS = UInt64(2_000_000_000)

"""
Backpressure ceiling on staged-but-undurable payload rows.

Background producers block once the buffer holds this many rows, bounding memory if processing
outruns the flusher. Priority (user-selected) deposits bypass the ceiling so interaction never
stalls.
"""
const CACHE_BUFFER_ROW_CEILING = 2_000_000

"""One interpreted source-item batch awaiting its flush: records and their aligned payloads."""
struct InterpretedDeposit
    records::Vector{ItemRecord}
    data::Vector{Any}            # aligned to `records`; an `AbstractDataItem`, or `nothing`
    rows::Int
end

"""Row count of one item payload, for backpressure accounting (0 when not a DataFrame)."""
function _payload_rows(item::Union{Nothing,AbstractDataItem})::Int
    item === nothing && return 0
    data = item_data(item)
    return data isa AbstractDataFrame ? nrow(data) : 0
end

"""Total payload rows across one interpreted deposit's data."""
function _interpreted_rows(data::Vector{Any})::Int
    total = 0
    for value in data
        value isa AbstractDataItem || continue
        payload = item_data(value)
        payload isa AbstractDataFrame && (total += nrow(payload))
    end
    return total
end

"""
In-memory staging in front of one [`CacheDB`](@ref).

Owns its lock and condition; `staged` is the read-through store, the `pending_*` collections are the
flusher's work list, and `staged_rows` drives backpressure. Holds no item objects after their flush
commits.
"""
mutable struct CacheBuffer
    cache::CacheDB
    profiler::Profiling.ProfileSession
    metrics::BuildMetrics
    lock::ReentrantLock
    condition::Base.Threads.Condition
    # Read-through store: every interpreted/processed item held in memory. Cacheable items are dropped
    # once their flush makes them durable on disk; non-cacheable items stay resident (disk never holds
    # them, and memory is the only place a worker can read them back).
    staged::Dict{Tuple{Symbol,String},AbstractDataItem}
    pending_interpreted::Vector{InterpretedDeposit}
    pending_processed::Vector{ProcessedWriteRequest}
    pending_failures::Vector{Pair{ItemRecord,Union{Nothing,String}}}
    staged_rows::Int
    row_ceiling::Int
    flush_interval_ns::UInt64
    last_flush_ns::UInt64
    flusher::Union{Nothing,Task}
    flush_requested::Bool
    draining::Bool
    closed::Bool
end

"""Construct an idle buffer; call [`start_cache_buffer!`](@ref) to launch its flusher."""
function CacheBuffer(
    cache::CacheDB,
    profiler::Profiling.ProfileSession,
    metrics::BuildMetrics;
    row_ceiling::Integer=CACHE_BUFFER_ROW_CEILING,
    flush_interval_ns::Integer=CACHE_FLUSH_INTERVAL_NS,
)::CacheBuffer
    buffer_lock = ReentrantLock()
    return CacheBuffer(
        cache,
        profiler,
        metrics,
        buffer_lock,
        Base.Threads.Condition(buffer_lock),
        Dict{Tuple{Symbol,String},AbstractDataItem}(),
        InterpretedDeposit[],
        ProcessedWriteRequest[],
        Pair{ItemRecord,Union{Nothing,String}}[],
        0,
        Int(row_ceiling),
        UInt64(flush_interval_ns),
        time_ns(),
        nothing,
        false,
        false,
        false,
    )
end

"""Whether the flusher has nothing left to write."""
function _pending_empty(buffer::CacheBuffer)::Bool
    return isempty(buffer.pending_interpreted) &&
           isempty(buffer.pending_processed) &&
           isempty(buffer.pending_failures)
end

# --------------------------------------------------------------------------------------------------
# Lifecycle
# --------------------------------------------------------------------------------------------------

"""Launch the dedicated flusher task that drains deposits into the cache writer."""
function start_cache_buffer!(buffer::CacheBuffer)::Nothing
    buffer.flusher === nothing || error("the cache buffer flusher is already running")
    lock(buffer.lock) do
        buffer.closed = false
        buffer.last_flush_ns = time_ns()
    end
    buffer.flusher = Base.Threads.@spawn cache_buffer_flusher!(buffer)
    return nothing
end

"""Drain everything still staged, stop the flusher, and wait for it to exit."""
function stop_cache_buffer!(buffer::CacheBuffer)::Nothing
    flusher = lock(buffer.lock) do
        buffer.closed = true
        notify(buffer.condition; all=true)
        buffer.flusher
    end
    flusher === nothing || wait(flusher)
    buffer.flusher = nothing
    return nothing
end

"""Whether the flusher still holds staged writes (a build is not idle until this is false)."""
function buffer_has_pending_writes(buffer::CacheBuffer)::Bool
    return lock(buffer.lock) do
        !_pending_empty(buffer) || buffer.draining
    end
end

"""Snapshot staged-write counts (queued items and payload rows) for memory diagnostics."""
function buffer_pending_counts(buffer::CacheBuffer)
    return lock(buffer.lock) do
        items = length(buffer.pending_processed)
        for deposit in buffer.pending_interpreted
            items += length(deposit.records)
        end
        (items=items, rows=buffer.staged_rows)
    end
end

"""
Block until staged interpreted writes are durable — the barrier before scan finalization.

Waits only for interpreted source-item rows, not for in-flight processed/stats writes: finalization
needs the source rows present, and processing keeps depositing long after the scan loop ends, so
waiting for everything would needlessly re-serialize the pipeline.
"""
function wait_cache_flushed!(buffer::CacheBuffer)::Nothing
    lock(buffer.lock) do
        buffer.flush_requested = true
        notify(buffer.condition; all=true)
        while !isempty(buffer.pending_interpreted) || buffer.draining
            wait(buffer.condition)
        end
    end
    return nothing
end

# --------------------------------------------------------------------------------------------------
# Depositing (the single write path)
# --------------------------------------------------------------------------------------------------

"""Block a background producer while the buffer is over its row ceiling."""
function _await_capacity!(buffer::CacheBuffer)::Nothing
    while buffer.staged_rows >= buffer.row_ceiling && !buffer.closed
        wait(buffer.condition)
    end
    return nothing
end

"""
Stage one interpreted source-item batch and make every interpreted item readable at once.

`records` and `data` align; each real item — cacheable or not — enters the read-through store, while a
`nothing` payload (a failed or invalidated interpretation) is staged only so the flush deletes any
stale on-disk copy. Blocks only when the buffer is over its ceiling.
"""
function stage_interpreted!(
    buffer::CacheBuffer,
    records::Vector{ItemRecord},
    data::Vector{Any},
)::Nothing
    length(records) == length(data) ||
        throw(DimensionMismatch("interpreted records and data must align"))
    rows = _interpreted_rows(data)
    lock(buffer.lock) do
        _await_capacity!(buffer)
        push!(buffer.pending_interpreted, InterpretedDeposit(records, data, rows))
        for (record, value) in zip(records, data)
            value isa AbstractDataItem || continue
            buffer.staged[(:interpreted, record.id)] = value
        end
        buffer.staged_rows += rows
        notify(buffer.condition; all=true)
    end
    return nothing
end

"""
Stage one processed result: an optional payload to cache plus optional statistics.

`priority` deposits (user-selected work) bypass backpressure. A cached payload is also placed in the
read-through store so the next view of that item is served from memory.
"""
function stage_processed!(
    buffer::CacheBuffer,
    record::ItemRecord,
    item::Union{Nothing,AbstractDataItem},
    stats::Union{Nothing,MetadataDict},
    write_payload::Bool;
    priority::Bool=false,
)::Nothing
    rows = write_payload ? _payload_rows(item) : 0
    request = ProcessedWriteRequest(record, item, stats, write_payload, rows)
    lock(buffer.lock) do
        priority || _await_capacity!(buffer)
        push!(buffer.pending_processed, request)
        if write_payload && item !== nothing
            buffer.staged[(:processed, record.id)] = item
        end
        buffer.staged_rows += rows
        notify(buffer.condition; all=true)
    end
    return nothing
end

"""Stage one item failure (or its clearance) for the next flush."""
function stage_failure!(
    buffer::CacheBuffer,
    record::ItemRecord,
    message::Union{Nothing,String},
)::Nothing
    lock(buffer.lock) do
        push!(buffer.pending_failures, record => message)
        notify(buffer.condition; all=true)
    end
    return nothing
end

# --------------------------------------------------------------------------------------------------
# Reading (read-through)
# --------------------------------------------------------------------------------------------------

"""
Read item data for one stage, serving staged (not-yet-durable) entries from memory first.

A drop-in for [`read_cached_item_data`](@ref): returns a vector aligned to `records`, each element an
`AbstractDataItem` bound to its record or `nothing` for a miss. Staged hits never touch DuckDB; misses
are fetched in one batched cache read.
"""
function buffer_read_item_data(
    buffer::CacheBuffer,
    records::Vector{ItemRecord};
    stage::Symbol=:processed,
)::Vector{Any}
    isempty(records) && return Any[]
    loaded = Vector{Any}(undef, length(records))
    misses = ItemRecord[]
    miss_positions = Int[]
    lock(buffer.lock) do
        for (position, record) in pairs(records)
            staged = get(buffer.staged, (stage, record.id), nothing)
            if staged === nothing
                loaded[position] = nothing
                push!(misses, record)
                push!(miss_positions, position)
            else
                # A cacheable item is returned as a disk read would reconstruct it (a DataItem); a
                # non-cacheable item is returned as-is so its concrete type survives — memory is the
                # only place its data will ever live.
                loaded[position] = cacheable(staged) ? DataItem(record, item_data(staged)) : staged
            end
        end
    end
    if !isempty(misses)
        fetched = read_cached_item_data(buffer.cache, misses; stage)
        for (index, position) in pairs(miss_positions)
            loaded[position] = fetched[index]
        end
    end
    return loaded
end

# --------------------------------------------------------------------------------------------------
# Flushing
# --------------------------------------------------------------------------------------------------

"""Wait for the next flush trigger: the interval boundary, a signal, or shutdown."""
function _wait_for_flush!(buffer::CacheBuffer)::Nothing
    if _pending_empty(buffer)
        # Nothing to write: sleep until a deposit, a flush request, or close signals us.
        wait(buffer.condition)
        return nothing
    end
    elapsed = time_ns() - buffer.last_flush_ns
    remaining = elapsed >= buffer.flush_interval_ns ? 0.0 :
        Float64(buffer.flush_interval_ns - elapsed) / 1e9
    timer = Timer(max(remaining, 0.0)) do _
        lock(buffer.lock) do
            notify(buffer.condition; all=true)
        end
    end
    try
        wait(buffer.condition)
    finally
        close(timer)
    end
    return nothing
end

"""Flusher loop: accumulate for one interval (or until pressure), then drain in one pass."""
function cache_buffer_flusher!(buffer::CacheBuffer)::Nothing
    while true
        action = lock(buffer.lock) do
            while true
                _pending_empty(buffer) && buffer.closed && return :stop
                due = !_pending_empty(buffer) && (
                    buffer.closed ||
                    buffer.flush_requested ||
                    buffer.staged_rows >= buffer.row_ceiling ||
                    (time_ns() - buffer.last_flush_ns) >= buffer.flush_interval_ns
                )
                due && return :flush
                _wait_for_flush!(buffer)
            end
        end
        action === :stop && break
        _flush_once!(buffer)
    end
    return nothing
end

"""Drain all currently pending deposits in one writer pass, then evict the durable entries."""
function _flush_once!(buffer::CacheBuffer)::Nothing
    interpreted, processed, failures, rows = lock(buffer.lock) do
        snapshot = (
            buffer.pending_interpreted,
            buffer.pending_processed,
            buffer.pending_failures,
            buffer.staged_rows,
        )
        buffer.pending_interpreted = InterpretedDeposit[]
        buffer.pending_processed = ProcessedWriteRequest[]
        buffer.pending_failures = Pair{ItemRecord,Union{Nothing,String}}[]
        buffer.flush_requested = false
        buffer.last_flush_ns = time_ns()
        buffer.draining = true
        snapshot
    end
    try
        # Each write is isolated so one failing group neither aborts the others nor — critically —
        # kills the flusher task and deadlocks every backpressured producer behind it.
        _try_write(() -> _write_interpreted!(buffer, interpreted), :interpreted, length(interpreted))
        _try_write(() -> _write_processed!(buffer, processed), :processed, length(processed))
        _try_write(() -> _write_failures!(buffer, failures), :failures, length(failures))
    finally
        lock(buffer.lock) do
            # Evict regardless of success: on failure the data is simply not durable, and a later read
            # recomputes it through the source fallback rather than serving a stale staged copy forever.
            _evict_interpreted!(buffer, interpreted)
            _evict_processed!(buffer, processed)
            buffer.staged_rows = max(0, buffer.staged_rows - rows)
            buffer.draining = false
            notify(buffer.condition; all=true)
        end
    end
    return nothing
end

"""Run one flush write, logging and swallowing failures so the flusher can never die."""
function _try_write(work::Function, stage::Symbol, count::Int)::Nothing
    try
        work()
    catch error
        error isa InterruptException && rethrow()
        @error(
            "Cache buffer flush failed; dropping a staged write batch",
            stage,
            items=count,
            exception=(error, catch_backtrace()),
        )
    end
    return nothing
end

"""Commit staged interpreted batches through the existing source-reconcile transaction."""
function _write_interpreted!(buffer::CacheBuffer, deposits::Vector{InterpretedDeposit})::Nothing
    isempty(deposits) && return nothing
    record_batches = Vector{Vector{ItemRecord}}(undef, length(deposits))
    data_batches = Vector{Vector{Any}}(undef, length(deposits))
    for (index, deposit) in pairs(deposits)
        record_batches[index] = deposit.records
        data_batches[index] = deposit.data
    end
    started = time_ns()
    reconcile_source_items!(buffer.cache, record_batches, data_batches; stage=:interpreted)
    record_cache_phase!(
        buffer.metrics.interpreted_write_ns, buffer.metrics.interpreted_writes, started)
    return nothing
end

"""Commit staged processed payloads and item statistics through the existing writers."""
function _write_processed!(buffer::CacheBuffer, batch::Vector{ProcessedWriteRequest})::Nothing
    isempty(batch) && return nothing
    payload = ProcessedWriteRequest[request for request in batch if request.write_payload]
    item_stats = Dict{String,MetadataDict}(
        request.record.id => (request.stats::MetadataDict)
        for request in batch if request.stats !== nothing
    )
    if !isempty(payload)
        started = time_ns()
        write_cached_item_data!(
            buffer.cache,
            ItemRecord[request.record for request in payload],
            Any[request.item for request in payload];
            stage=:processed,
        )
        record_cache_phase!(
            buffer.metrics.processed_write_ns, buffer.metrics.processed_writes, started)
    end
    if !isempty(item_stats)
        started = time_ns()
        persist_stats!(
            buffer.cache, item_stats, Dict{Tuple{Vararg{String}},Dict{Symbol,Any}}())
        record_cache_phase!(buffer.metrics.stats_write_ns, buffer.metrics.stats_writes, started)
    end
    return nothing
end

"""Commit staged item failures."""
function _write_failures!(
    buffer::CacheBuffer,
    failures::Vector{Pair{ItemRecord,Union{Nothing,String}}},
)::Nothing
    for (record, message) in failures
        persist_item_failure!(buffer.cache, record, message)
    end
    return nothing
end

"""
Drop interpreted read-through entries that are now durable on disk (identity-checked).

Only cacheable items are evicted: their data is on disk, so reads fall through to DuckDB. Non-cacheable
items will never be persisted, so they stay memory-resident — evicting them would force a redundant
source re-read the instant a worker asked for them.
"""
function _evict_interpreted!(buffer::CacheBuffer, deposits::Vector{InterpretedDeposit})::Nothing
    for deposit in deposits
        for (record, value) in zip(deposit.records, deposit.data)
            value isa AbstractDataItem || continue
            cacheable(value) || continue
            key = (:interpreted, record.id)
            get(buffer.staged, key, nothing) === value && delete!(buffer.staged, key)
        end
    end
    return nothing
end

"""Drop processed read-through entries whose payload we just committed (identity-checked)."""
function _evict_processed!(buffer::CacheBuffer, batch::Vector{ProcessedWriteRequest})::Nothing
    for request in batch
        (request.write_payload && request.item !== nothing) || continue
        key = (:processed, request.record.id)
        get(buffer.staged, key, nothing) === request.item && delete!(buffer.staged, key)
    end
    return nothing
end
