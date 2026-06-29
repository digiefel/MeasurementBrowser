# --------------------------------------------------------------------------------------------------
# Per-table cache buffer
#
# One buffer sits in front of each database table (see docs/cache.md). Each disk-backed buffer owns
# its own write connection, its own read connection, and its own background flush task, so no two
# buffers ever write the same table and a write–write conflict is impossible by construction. One
# memory-only buffer holds interpreted payloads, which are never written to disk.
#
# Scope (first cut): the streaming write path is APPEND-ONLY. `edit`/`delete` are the intended buffer
# verbs (docs/cache.md "Mutations, deferred") but are not built yet; a build clears the cache once up
# front and then only appends. The pending store is keyed by row identity, so a re-stored key
# coalesces in place — the append-only path is the first step of that coalescing model, not a
# throwaway.
# --------------------------------------------------------------------------------------------------

# Rows that actually cost memory (payload data rows) bound the buffer; bookkeeping rows are tiny and
# effectively never reach this. Tunable against the benchmark in Step 4.
const CACHE_BUFFER_ROW_CEILING = 2_000_000
const CACHE_BUFFER_BOOKKEEPING_CEILING = typemax(Int)
const CACHE_BUFFER_FLUSH_INTERVAL = 2.0

# Internal measurement columns ride alongside the user's columns; the `__mb_` prefix is reserved.
const MB_SEQ_COLUMN = "__mb_seq"
const MB_ROW_COLUMN = "__mb_row"

# --------------------------------------------------------------------------------------------------
# Generic SQL and payload-table helpers (mechanism: identifiers, shapes, table lifecycle)
# --------------------------------------------------------------------------------------------------

"""Quote an SQL identifier without treating user-provided column names as SQL."""
function _quote_identifier(identifier::AbstractString)::String
    return "\"" * replace(String(identifier), "\"" => "\"\"") * "\""
end

"""Return a stable id for one stage and a DataFrame's ordered column names and element types."""
function _dataframe_storage_id(data::AbstractDataFrame, stage::Symbol)::String
    schema = IOBuffer()
    print(schema, String(stage), ':')
    for name in names(data)
        type_name = string(eltype(data[!, name]))
        print(schema, ncodeunits(name), ':', name, ncodeunits(type_name), ':', type_name)
    end
    return bytes2hex(sha1(take!(schema)))
end

"""Return the physical table name for a validated DataFrame schema id."""
function _dataframe_table_name(storage_id::AbstractString)::String
    id = String(storage_id)
    occursin(r"^[0-9a-f]{40}$", id) ||
        error("Invalid DataFrame cache storage id '$id'; expected a SHA-1 digest")
    return "dataframe_$id"
end

"""Read the known DataFrame schemas into memory once for a cache operation."""
function _load_dataframe_schemas(connection)::Dict{String,Vector{String}}
    schemas = Dict{String,Vector{String}}()
    for row in DBInterface.execute(
        connection, "SELECT storage_id, column_names FROM dataframe_schemas")
        schemas[String(row.storage_id)] = String[name for name in row.column_names]
    end
    return schemas
end

"""Drop a physical DataFrame table after its final item-data entry has been removed."""
function _drop_unused_dataframe_storage!(
    connection,
    storage_id::AbstractString,
)::Nothing
    statement = DBInterface.prepare(
        connection, "SELECT count(*) AS count FROM item_data WHERE storage_id = ?")
    count = only(DBInterface.execute(statement, (String(storage_id),))).count
    count == 0 || return nothing
    DBInterface.execute(
        connection, "DROP TABLE $(_quote_identifier(_dataframe_table_name(storage_id)))")
    DBInterface.execute(
        DBInterface.prepare(connection, "DELETE FROM dataframe_schemas WHERE storage_id = ?"),
        (String(storage_id),),
    )
    return nothing
end

"""Delete all cached data derived from one source item."""
function _delete_cached_source_item_data!(
    connection,
    source_item_id::AbstractString,
)::Nothing
    storage_ids = String[
        String(row.storage_id)
        for row in DBInterface.execute(
            DBInterface.prepare(connection, """
                SELECT DISTINCT storage_id FROM item_data WHERE source_item_id = ?
            """),
            (String(source_item_id),),
        )
    ]
    for storage_id in storage_ids
        table = _quote_identifier(_dataframe_table_name(storage_id))
        DBInterface.execute(
            DBInterface.prepare(connection, """
                DELETE FROM $table
                WHERE __mb_seq IN (
                    SELECT seq FROM item_data
                    WHERE source_item_id = ? AND storage_id = ?
                )
            """),
            (String(source_item_id), storage_id),
        )
    end
    DBInterface.execute(
        DBInterface.prepare(connection, "DELETE FROM item_data WHERE source_item_id = ?"),
        (String(source_item_id),),
    )
    for storage_id in storage_ids
        _drop_unused_dataframe_storage!(connection, storage_id)
    end
    return nothing
end

"""Delete every cached DataFrame and its schema/index rows."""
function _delete_all_cached_item_data!(connection)::Nothing
    for storage_id in keys(_load_dataframe_schemas(connection))
        DBInterface.execute(
            connection, "DROP TABLE $(_quote_identifier(_dataframe_table_name(storage_id)))")
    end
    DBInterface.execute(connection, "DELETE FROM item_data")
    DBInterface.execute(connection, "DELETE FROM dataframe_schemas")
    return nothing
end

# --------------------------------------------------------------------------------------------------
# Per-table row descriptions (the only table-aware code)
# --------------------------------------------------------------------------------------------------

"""One `source_items` row: one per discovered source item, keyed by its id."""
struct SourceItemRow
    id::String
    fingerprint_hex::Union{Nothing,String}
    fingerprint_hash::Union{Nothing,String}
    path::Union{Nothing,String}
    timestamp::Union{Nothing,DateTime}
end

"""One `items` row: one data-less logical-item record, keyed by item id."""
struct ItemRow
    id::String
    source_item_id::String
    item_label::String
    kind::String
    collection::Vector{String}
    item_fingerprint_hex::Union{Nothing,String}
end

"""One `metadata` (EAV) row, keyed by `(scope, entity, key)`."""
struct MetaRow
    scope::Int8
    entity::String
    key::String
    value::MetadataValue
end

"""One `item_failures` row, keyed by item id."""
struct FailureRow
    item_id::String
    source_item_id::String
    message::String
end

"""One cacheable processed DataFrame payload awaiting its measurement table + pointer, keyed by id."""
struct PayloadEntry
    record::ItemRecord
    data::AbstractDataFrame
    storage_id::String
    column_names::Vector{String}
end

# The memory weight one stored row contributes to the ceiling: payload data rows for the heavy ones,
# a single unit for tiny bookkeeping rows.
_row_weight(::Any)::Int = 1
_row_weight(entry::PayloadEntry)::Int = nrow(entry.data)
function _row_weight(item::AbstractDataItem)::Int
    data = item_data(item)
    return data isa AbstractDataFrame ? nrow(data) : 1
end

# --------------------------------------------------------------------------------------------------
# The generic table buffer
# --------------------------------------------------------------------------------------------------

"""
Bounded memory in front of one database table (or, when `memory_only`, in front of nothing).

`entries` is the coalescing pending store keyed by row identity; a re-stored key overwrites in place.
A disk-backed buffer flushes `entries` to its table on its own schedule and serves reads from memory
or that table; a memory-only buffer never flushes and drops on overflow instead of blocking.
"""
mutable struct CacheTableBuffer{R}
    name::Symbol
    write_conn::Union{Nothing,DuckDB.Connection}
    read_conn::Union{Nothing,DuckDB.Connection}
    read_lock::ReentrantLock
    known_schemas::Set{String}          # measurement-table storage ids already created (payload only)
    entries::Dict{Any,R}
    pending_items::Int
    pending_rows::Int
    cond::Base.Threads.Condition
    ceiling::Int
    counts_rows::Bool                   # ceiling counts payload rows (true) or whole entries (false)
    memory_only::Bool
    drop_on_error::Bool                 # drop a failed flush batch (recomputable) vs. retry it
    interval::Float64
    last_flush::Float64
    flush_requested::Bool
    closing::Bool
    task::Union{Nothing,Task}
end

function CacheTableBuffer{R}(
    name::Symbol;
    memory_only::Bool=false,
    counts_rows::Bool=false,
    drop_on_error::Bool=false,
    ceiling::Int=CACHE_BUFFER_BOOKKEEPING_CEILING,
)::CacheTableBuffer{R} where {R}
    return CacheTableBuffer{R}(
        name, nothing, nothing, ReentrantLock(), Set{String}(), Dict{Any,R}(), 0, 0,
        Base.Threads.Condition(), ceiling, counts_rows, memory_only, drop_on_error,
        CACHE_BUFFER_FLUSH_INTERVAL, time(), false, false, nothing)
end

_load(buffer::CacheTableBuffer)::Int = buffer.counts_rows ? buffer.pending_rows : buffer.pending_items

"""
Store one row, coalescing on `key`. Disk-backed buffers backpressure over the ceiling; memory-only
buffers never block and drop the incoming row instead. Returns whether the row was retained.
"""
function _store!(buffer::CacheTableBuffer{R}, key, row::R)::Bool where {R}
    lock(buffer.cond)
    try
        if buffer.memory_only
            haskey(buffer.entries, key) || _load(buffer) < buffer.ceiling || return false
        else
            while _load(buffer) >= buffer.ceiling && !haskey(buffer.entries, key)
                wait(buffer.cond)
            end
        end
        existing = get(buffer.entries, key, nothing)
        if existing === nothing
            buffer.pending_items += 1
        else
            buffer.pending_rows -= _row_weight(existing)
        end
        buffer.entries[key] = row
        buffer.pending_rows += _row_weight(row)
        notify(buffer.cond)
        return true
    finally
        unlock(buffer.cond)
    end
end

"""Whether the buffer should flush now (called while holding `cond`)."""
function _flush_due(buffer::CacheTableBuffer)::Bool
    isempty(buffer.entries) && return false
    return buffer.closing || buffer.flush_requested ||
        _load(buffer) >= buffer.ceiling ||
        (time() - buffer.last_flush) >= buffer.interval
end

"""Run one disk-backed buffer's flush task until it is closed and drained."""
function _run_flusher!(buffer::CacheTableBuffer)::Nothing
    while true
        lock(buffer.cond)
        try
            while !_flush_due(buffer)
                if buffer.closing && isempty(buffer.entries)
                    return nothing
                end
                # Wake at the latest when the flush interval elapses; a store notifies sooner.
                timer = Timer(buffer.interval) do _
                    lock(() -> notify(buffer.cond), buffer.cond)
                end
                try
                    wait(buffer.cond)
                finally
                    close(timer)
                end
            end
        finally
            unlock(buffer.cond)
        end
        _flush_once!(buffer)
        lock(buffer.cond) do
            buffer.flush_requested = false
            notify(buffer.cond)
        end
    end
    return nothing
end

"""Flush every pending row in one transaction, then evict the rows that committed."""
function _flush_once!(buffer::CacheTableBuffer{R})::Nothing where {R}
    snapshot = Tuple{Any,R}[]
    lock(buffer.cond) do
        for (key, row) in buffer.entries
            push!(snapshot, (key, row))
        end
        buffer.last_flush = time()
    end
    isempty(snapshot) && return nothing

    rows = R[row for (_, row) in snapshot]
    connection = buffer.write_conn
    committed = false
    DBInterface.execute(connection, "BEGIN TRANSACTION")
    try
        _flush_rows!(rows, connection, buffer)
        DBInterface.execute(connection, "COMMIT")
        committed = true
    catch error
        error isa InterruptException && rethrow()
        try
            DBInterface.execute(connection, "ROLLBACK")
        catch
        end
        # A dead flush task would strand blocked stores, so never let one out (rule 6). Recomputable
        # data (payloads) is dropped and re-derived on the next read; durable bookkeeping is left
        # pending to retry — the rollback discarded it, so re-appending cannot duplicate.
        @warn "Cache flush failed; $(buffer.drop_on_error ? "dropping" : "retrying") batch" buffer=buffer.name exception=error
    end

    (committed || buffer.drop_on_error) || return nothing
    lock(buffer.cond) do
        for (key, row) in snapshot
            buffer.entries[key] === row || continue   # a newer store replaced it; keep that one
            delete!(buffer.entries, key)
            buffer.pending_items -= 1
            buffer.pending_rows -= _row_weight(row)
        end
        notify(buffer.cond)
    end
    return nothing
end

"""Request a flush and block until the buffer is empty (a one-shot drain barrier)."""
function _flush_and_wait!(buffer::CacheTableBuffer)::Nothing
    lock(buffer.cond) do
        buffer.flush_requested = true
        notify(buffer.cond)
        while !isempty(buffer.entries)
            wait(buffer.cond)
        end
    end
    return nothing
end

# --------------------------------------------------------------------------------------------------
# Per-table flush (the generic machinery never names a table; these do)
# --------------------------------------------------------------------------------------------------

"""Append rows to `table` through one appender, closed defensively."""
function _append_table!(emit, connection, table::AbstractString)::Nothing
    appender = DuckDB.Appender(connection, table)
    try
        emit(appender)
        DuckDB.flush(appender)
    finally
        # Closing flushes; on an aborted transaction that itself throws and would mask the real error.
        try
            DuckDB.close(appender)
        catch
        end
    end
    return nothing
end

function _flush_rows!(rows::Vector{SourceItemRow}, connection, ::CacheTableBuffer)::Nothing
    _append_table!(connection, "source_items") do appender
        for row in rows
            for value in (row.id, row.fingerprint_hex, row.fingerprint_hash,
                          row.path, row.timestamp, nothing)
                DuckDB.append(appender, value)
            end
            DuckDB.end_row(appender)
        end
    end
    return nothing
end

function _flush_rows!(rows::Vector{ItemRow}, connection, ::CacheTableBuffer)::Nothing
    _append_table!(connection, "items") do appender
        for row in rows
            for value in (row.id, row.source_item_id, row.item_label, row.kind,
                          row.collection, row.item_fingerprint_hex)
                DuckDB.append(appender, value)
            end
            DuckDB.end_row(appender)
        end
    end
    return nothing
end

function _flush_rows!(rows::Vector{MetaRow}, connection, ::CacheTableBuffer)::Nothing
    _append_table!(connection, "metadata") do appender
        for row in rows
            codec = _meta_codec(row.value)
            stored = codec.to_db(row.value)
            DuckDB.append(appender, row.scope)
            DuckDB.append(appender, row.entity)
            DuckDB.append(appender, row.key)
            DuckDB.append(appender, Int8(codec.vtype))
            for column in META_VALUE_COLUMNS
                DuckDB.append(appender, column === codec.column ? stored : missing)
            end
            DuckDB.end_row(appender)
        end
    end
    return nothing
end

function _flush_rows!(rows::Vector{FailureRow}, connection, ::CacheTableBuffer)::Nothing
    _append_table!(connection, "item_failures") do appender
        for row in rows
            for value in (row.item_id, row.source_item_id, row.message)
                DuckDB.append(appender, value)
            end
            DuckDB.end_row(appender)
        end
    end
    return nothing
end

"""Map a Julia element type to the DuckDB column type used in a measurement table's `CREATE TABLE`."""
function _duckdb_sql_type(T::Type)::String
    T = Base.nonmissingtype(T)
    T <: Bool && return "BOOLEAN"
    T <: Int8 && return "TINYINT"
    T <: Int16 && return "SMALLINT"
    T <: Int32 && return "INTEGER"
    T <: Int64 && return "BIGINT"
    T <: Int128 && return "HUGEINT"
    T <: UInt8 && return "UTINYINT"
    T <: UInt16 && return "USMALLINT"
    T <: UInt32 && return "UINTEGER"
    T <: UInt64 && return "UBIGINT"
    T <: Float32 && return "FLOAT"
    T <: Float64 && return "DOUBLE"
    T <: AbstractString && return "VARCHAR"
    T <: DateTime && return "TIMESTAMP"
    T <: Date && return "DATE"
    error("No DuckDB column type for Julia element type $T")
end

"""Mint one fresh integer surrogate per id from the catalog sequence in a single round trip."""
function _mint_seqs(connection, ids::Vector{String})::Dict{String,Int64}
    seqs = Dict{String,Int64}()
    isempty(ids) && return seqs
    index = 0
    for row in DBInterface.execute(
        connection, "SELECT nextval('item_data_seq') AS seq FROM range($(length(ids)))")
        index += 1
        seqs[ids[index]] = Int64(row.seq)
    end
    return seqs
end

"""
Write cacheable processed payloads and their pointer rows in one transaction.

The owning buffer holds both the per-shape measurement tables and the `item_data` pointer table, so a
reader can never see a pointer whose payload is not yet on disk. A shape's table is created on first
sight; thereafter rows are appended directly.
"""
function _flush_rows!(rows::Vector{PayloadEntry}, connection, buffer::CacheTableBuffer)::Nothing
    seqs = _mint_seqs(connection, String[entry.record.id for entry in rows])
    groups = Dict{String,Vector{PayloadEntry}}()
    for entry in rows
        push!(get!(() -> PayloadEntry[], groups, entry.storage_id), entry)
    end
    for (storage_id, entries) in groups
        columns = first(entries).column_names
        table = _dataframe_table_name(storage_id)
        if !(storage_id in buffer.known_schemas)
            definitions = join(
                ("$(_quote_identifier(name)) " *
                 "$(_duckdb_sql_type(eltype(first(entries).data[!, name])))" for name in columns),
                ", ")
            DBInterface.execute(connection, "CREATE TABLE $(_quote_identifier(table)) (" *
                "$MB_SEQ_COLUMN BIGINT, $MB_ROW_COLUMN BIGINT, $definitions)")
            DBInterface.execute(
                DBInterface.prepare(connection, "INSERT INTO dataframe_schemas VALUES (?, ?)"),
                (storage_id, columns))
            push!(buffer.known_schemas, storage_id)
        end
        _append_table!(connection, table) do appender
            for entry in entries
                seq = seqs[entry.record.id]
                column_vectors = [entry.data[!, name] for name in columns]
                for row_index in 1:nrow(entry.data)
                    DuckDB.append(appender, seq)
                    DuckDB.append(appender, Int64(row_index))
                    for column in column_vectors
                        DuckDB.append(appender, column[row_index])
                    end
                    DuckDB.end_row(appender)
                end
            end
        end
    end
    _append_table!(connection, "item_data") do appender
        for entry in rows
            record = entry.record
            for value in (
                record.id, "processed", record.source_item_id,
                _fingerprint_hash(record.source_item_fingerprint),
                _fingerprint_hash(record.item_fingerprint),
                entry.storage_id, Int64(nrow(entry.data)), seqs[record.id],
            )
                DuckDB.append(appender, value)
            end
            DuckDB.end_row(appender)
        end
    end
    return nothing
end

# --------------------------------------------------------------------------------------------------
# Reading payloads back (memory first, then the buffer's read connection)
# --------------------------------------------------------------------------------------------------

"""Whether a stored fingerprint hash matches the requested one, treating SQL NULL as `nothing`."""
function _hash_matches(stored, requested)::Bool
    stored = _null_to_nothing(stored)
    stored === nothing && return requested === nothing
    return requested !== nothing && String(stored) == requested
end

"""Reconstruct processed payloads not held in memory from the measurement tables on disk."""
function _read_payloads_from_disk(
    buffer::CacheTableBuffer,
    records::Vector{ItemRecord},
)::Vector{Any}
    loaded = Any[nothing for _ in records]
    requested = Dict{String,Vector{Int}}()
    ids = String[]
    for (position, record) in pairs(records)
        record.source_item_fingerprint === nothing && continue
        positions = get!(requested, record.id) do
            push!(ids, record.id)
            Int[]
        end
        push!(positions, position)
    end
    isempty(ids) && return loaded

    lock(buffer.read_lock) do
        connection = buffer.read_conn
        placeholders = join(fill("?", length(ids)), ", ")
        by_storage = Dict{String,Tuple{Vector{String},Vector{Tuple{Int64,Vector{Int},Int}}}}()
        statement = DBInterface.prepare(connection, """
            SELECT d.item_id, d.sif_hash, d.if_hash, d.storage_id, d.row_count, d.seq, s.column_names
            FROM item_data d JOIN dataframe_schemas s ON s.storage_id = d.storage_id
            WHERE d.stage = 'processed' AND d.item_id IN ($placeholders)
        """)
        for row in DBInterface.execute(statement, Tuple(ids))
            id = String(row.item_id)
            positions = Int[
                position
                for position in get(requested, id, Int[])
                if _hash_matches(row.sif_hash,
                       _fingerprint_hash(records[position].source_item_fingerprint)) &&
                   _hash_matches(row.if_hash,
                       _fingerprint_hash(records[position].item_fingerprint))
            ]
            isempty(positions) && continue
            columns = String[String(name) for name in row.column_names]
            _, entries = get!(by_storage, String(row.storage_id)) do
                (columns, Tuple{Int64,Vector{Int},Int}[])
            end
            push!(entries, (Int64(row.seq), positions, Int(row.row_count)))
        end

        for (storage_id, (columns, entries)) in by_storage
            table = _quote_identifier(_dataframe_table_name(storage_id))
            selected = join((_quote_identifier(name) for name in columns), ", ")
            placeholders = join(fill("?", length(entries)), ", ")
            statement = DBInterface.prepare(connection, """
                SELECT $MB_SEQ_COLUMN, $selected FROM $table
                WHERE $MB_SEQ_COLUMN IN ($placeholders)
                ORDER BY $MB_SEQ_COLUMN, $MB_ROW_COLUMN
            """)
            raw = DataFrame(DBInterface.execute(statement, Tuple(entry[1] for entry in entries)))
            expected = Dict(entry[1] => (entry[2], entry[3]) for entry in entries)

            rebuild(range) = DataFrame(
                (Symbol(name) => view(raw[!, Symbol(name)], range) for name in columns)...;
                copycols=false)

            seq_column = raw[!, Symbol(MB_SEQ_COLUMN)]
            seen = Set{Int64}()
            first_row = 1
            total = nrow(raw)
            while first_row <= total
                seq = Int64(seq_column[first_row])
                last_row = first_row
                while last_row < total && Int64(seq_column[last_row + 1]) == seq
                    last_row += 1
                end
                positions, _ = expected[seq]
                data = rebuild(first_row:last_row)
                for position in positions
                    loaded[position] = DataItem(records[position], data)
                end
                push!(seen, seq)
                first_row = last_row + 1
            end
            # Zero-row payloads never appear in the table; rebuild them from the typed empty columns.
            for (seq, (positions, row_count)) in expected
                seq in seen && continue
                row_count == 0 || error("Cached payload for seq $seq is missing its $row_count rows")
                data = rebuild(1:0)
                for position in positions
                    loaded[position] = DataItem(records[position], data)
                end
            end
        end
    end
    return loaded
end

# --------------------------------------------------------------------------------------------------
# The owner: one buffer per table behind a single front door
# --------------------------------------------------------------------------------------------------

"""
The single door to the cache: a buffer per table plus the memory-only interpreted buffer.

Every per-item write and read in the workspace goes through this owner; nothing else opens a
connection or runs SQL against the cache during a build. Holds only the raw `DuckDB.DB` handle — it
knows nothing of `CacheDB` or the pipeline stages; `ProjectCache` owns and drives it.
"""
mutable struct CacheBuffer
    db::DuckDB.DB
    metrics::BuildMetrics
    source_items::CacheTableBuffer{SourceItemRow}
    items::CacheTableBuffer{ItemRow}
    metadata::CacheTableBuffer{MetaRow}
    failures::CacheTableBuffer{FailureRow}
    payload::CacheTableBuffer{PayloadEntry}
    interpreted::CacheTableBuffer{AbstractDataItem}
    # Processed payloads that `cacheable` opts out of disk storage for. Like interpreted data they are
    # memory-only and best-effort: held so a re-selection skips reprocessing, bounded by the row
    # ceiling, dropped on overflow, and recomputed from source on a miss (docs/cache.md).
    processed_memory::CacheTableBuffer{AbstractDataItem}
end

function CacheBuffer(
    db::DuckDB.DB,
    metrics::BuildMetrics,
)::CacheBuffer
    return CacheBuffer(
        db, metrics,
        CacheTableBuffer{SourceItemRow}(:source_items),
        CacheTableBuffer{ItemRow}(:items),
        CacheTableBuffer{MetaRow}(:metadata),
        CacheTableBuffer{FailureRow}(:item_failures),
        CacheTableBuffer{PayloadEntry}(
            :item_data; counts_rows=true, drop_on_error=true, ceiling=CACHE_BUFFER_ROW_CEILING),
        CacheTableBuffer{AbstractDataItem}(
            :interpreted; memory_only=true, counts_rows=true, ceiling=CACHE_BUFFER_ROW_CEILING),
        CacheTableBuffer{AbstractDataItem}(
            :processed_memory; memory_only=true, counts_rows=true, ceiling=CACHE_BUFFER_ROW_CEILING))
end

_disk_buffers(buffer::CacheBuffer) =
    (buffer.source_items, buffer.items, buffer.metadata, buffer.failures, buffer.payload)

"""Open each disk-backed buffer's connections and start its flush task."""
function start_cache_buffer!(buffer::CacheBuffer)::Nothing
    database = buffer.db
    for table in _disk_buffers(buffer)
        table.write_conn = DBInterface.connect(database)
        table.last_flush = time()
        table.task = Base.Threads.@spawn _run_flusher!(table)
    end
    # The payload buffer is the only one that reads back, and it must know which measurement-table
    # shapes already exist so it creates each exactly once.
    buffer.payload.read_conn = DBInterface.connect(database)
    for row in DBInterface.execute(buffer.payload.write_conn,
                                   "SELECT storage_id FROM dataframe_schemas")
        push!(buffer.payload.known_schemas, String(row.storage_id))
    end
    return nothing
end

"""Drain every disk-backed buffer, stop its task, and close its connections (not the database)."""
function stop_cache_buffer!(buffer::CacheBuffer)::Nothing
    for table in _disk_buffers(buffer)
        lock(table.cond) do
            table.closing = true
            notify(table.cond)
        end
        table.task === nothing || wait(table.task)
        table.write_conn === nothing || DBInterface.close!(table.write_conn)
        table.read_conn === nothing || DBInterface.close!(table.read_conn)
    end
    return nothing
end

"""Block until every disk-backed buffer has flushed its pending rows to disk."""
function wait_cache_flushed!(buffer::CacheBuffer)::Nothing
    for table in _disk_buffers(buffer)
        _flush_and_wait!(table)
    end
    return nothing
end

"""Whether any durable write is still staged."""
function buffer_has_pending_writes(buffer::CacheBuffer)::Bool
    return any(_disk_buffers(buffer)) do table
        lock(() -> !isempty(table.entries), table.cond)
    end
end

"""
Run a one-shot maintenance mutation on a transient write connection the buffer opens and closes.

The buffer is the single owner of every connection to the cache file. Per-table streaming writes go
through the per-table buffers above; the rare bulk/DDL mutations that span tables — stamping identity,
wiping for a rebuild, deleting the source items a scan dropped, the closing checkpoint — run here
instead of on a second long-lived writer. Each only ever runs while the per-table buffers are quiescent
(identity and the deletes are opening moves before any append; the checkpoint runs after the buffers
stop), so a fresh connection sees committed state and can never race a flush. No connection outlives the
call, so there is no persistent second write door.
"""
function with_maintenance(work::Function, buffer::CacheBuffer)::Any
    connection = DBInterface.connect(buffer.db)
    try
        return work(connection)
    finally
        DBInterface.close!(connection)
    end
end

"""Aggregate pending counts: `items` staged durable rows, `rows` staged payload data rows."""
function buffer_pending_counts(buffer::CacheBuffer)
    items = 0
    rows = 0
    for table in _disk_buffers(buffer)
        lock(table.cond) do
            items += table.pending_items
            rows += table.pending_rows
        end
    end
    return (items=items, rows=rows)
end

# The domain write/read surface (stage mapping, ItemRecord↔rows, stats) lives in ProjectCache, which
# owns this buffer. CacheBuffer keeps only the generic per-table machinery above.
