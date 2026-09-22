"""Record a discovered source item's identity and fingerprint, independently of its stage results."""
function store_source_identity!(
    cache::CacheDB, item::AbstractDataSourceItem, display_label::AbstractString=label(item),
)::Nothing
    source_key = source_item_key!(cache, id(item); mint=true)
    edit!(cache.source_items, id(item), SourceItemRow(
        id(item), source_key, String(display_label), _serialize_hex(fingerprint(item)),
        source_item_path(item), source_item_timestamp(item)))
    _stage_ledger_source!(cache.stage_ledger, source_key, true)
    return nothing
end

"""Record successful completion of a source stage after publishing its output."""
function store_source_result!(cache::CacheDB, stage::PipelineStage, source_key::Int64)::Nothing
    state = CachedKeyedResultState(Int8(stage), source_key, Int8(RESULT_READY), source_key, nothing)
    edit!(cache.keyed_result_states, (Int8(stage), source_key), state)
    _stage_ledger_result!(cache.stage_ledger, (Int8(stage), source_key), state)
    return nothing
end

"""
Cache a `read` result in memory and publish its completion.

Read results have no disk persistence by default. They use the shared memory cache with a 256 MiB
FIFO budget measured by `Base.summarysize`, retaining a single oversized value. The estimate covers
Julia-owned memory; memory owned by external libraries is outside that estimate. Eviction removes
the payload, not the fact that reading succeeded; consumers that need the value again must reread. A scheduled consumer holds
its own reference until it finishes, independently of cache eviction. Values, including `nothing`,
are passed through unchanged; the cache does not serialize or close user-owned handles.
"""
function store_source_read!(cache::CacheDB, source_key::Int64, value)::Nothing
    append!(cache.source_reads, source_key, value)
    store_source_result!(cache, SOURCE_READ, source_key)
    return nothing
end

"""Read a source-stage payload by private source key; `Some(value)` is a hit, `nothing` a miss."""
function read_payload(cache::CacheDB, source_key::Integer; stage::PipelineStage=SOURCE_READ)
    stage === SOURCE_READ || throw(ArgumentError("Source payloads are available at SOURCE_READ"))
    return read_hit(cache.source_reads, Int64(source_key))
end

"""Whether the memory cache retains a source read result, without loading it."""
function has_payload(cache::CacheDB, source_key::Integer; stage::PipelineStage=SOURCE_READ)::Bool
    stage === SOURCE_READ || return false
    return haskey(cache.source_reads, Int64(source_key))
end

"""Whether a source has finished interpretation or failed at either source stage."""
function source_complete(cache::CacheDB, source_key::Integer)::Bool
    interpreted = cached_result_state(cache, SOURCE_INTERPRET, source_key)
    interpreted === nothing || return true
    loaded = cached_result_state(cache, SOURCE_READ, source_key)
    return loaded !== nothing && CacheResultStatus(loaded.status) === RESULT_FAILED
end

"""Whether the cache knows a source identity, including an unfinished or failed source."""
function cache_knows_source(cache::CacheDB, source_key::Integer)::Bool
    return lock(cache.stage_ledger.lock) do
        Int64(source_key) in cache.stage_ledger.source_items
    end
end

"""Delete a source identity and all its stage results, including its memory-only read result."""
function delete_source_item!(cache::CacheDB, source_key::Integer,
        old_records::Vector{ItemRecord})::Nothing
    delete_source_output!(cache, source_key, old_records)
    delete!(cache.source_reads, Int64(source_key))
    clear_cached_result_state!(cache, SOURCE_READ, source_key)
    _delete_source_identity!(cache, Int64(source_key))
    _stage_ledger_source!(cache.stage_ledger, source_key, false)
    return nothing
end

function _delete_source_identity!(cache::CacheDB, source_key::Int64)::Nothing
    delete!(cache.source_items, source_item_id(cache, source_key))
    return nothing
end
