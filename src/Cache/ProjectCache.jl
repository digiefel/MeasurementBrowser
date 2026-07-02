const PROJECT_CACHE_SCHEMA_VERSION = 9

"""
DuckDB buffer-pool limit (MiB) for cache connections.

DuckDB otherwise defaults to most of system RAM and retains committed table blocks for the workspace's
lifetime; cache writes are bounded batches, so a smaller buffer keeps large caches out of swap without
reducing scan parallelism. Settable at runtime — e.g. from a GUI slider — via
[`set_cache_memory_limit!`](@ref). Changes apply to caches opened afterward.
"""
const CACHE_MEMORY_LIMIT_MIB = Ref(1024)

_memory_limit_sql(mib::Integer)::String = "SET memory_limit = '$(Int(mib)) MiB'"

"""Set the default cache buffer-pool limit (MiB) for caches opened later."""
function set_cache_memory_limit!(mib::Integer)::Int
    mib >= 1 || throw(ArgumentError("cache memory limit must be at least 1 MiB"))
    return CACHE_MEMORY_LIMIT_MIB[] = Int(mib)
end

"""Which dict one metadata (EAV) row belongs to. Stored in the `scope` column."""
@enum MetaScope::Int8 begin
    SCOPE_ITEM_PARAMETERS = 0
    SCOPE_ITEM_STATS = 1
    SCOPE_NODE_PARAMETERS = 2
    SCOPE_NODE_STATS = 3
end

"""Discriminator naming which [`MetadataValue`](@ref) variant a metadata row holds."""
@enum MetaVType::Int8 begin
    VT_BOOL = 1
    VT_INT = 2
    VT_FLOAT = 3
    VT_STRING = 4
    VT_SYMBOL = 5
    VT_DATE = 6
    VT_DATETIME = 7
    VT_MISSING = 8
    VT_VBOOL = 9
    VT_VINT = 10
    VT_VFLOAT = 11
    VT_VSTRING = 12
end

"""Error raised when a project cache is missing or does not match the current project."""
struct ProjectCacheError <: Exception
    path::String
    message::String
end

"""Print the cache path and validation failure."""
Base.showerror(io::IO, err::ProjectCacheError)::Nothing =
    print(io, "Invalid project cache $(err.path): $(err.message)")

"""The project, source, and DuckDB file belonging to one cache."""
struct ProjectCacheIdentity
    project_name::String
    source_id::String
    source_label::String
    cache_path::String
end

@enum CacheResultKind::Int8 begin
    PROCESSING_RESULT = 1
    ITEM_STATS_RESULT = 2
    COLLECTION_STATS_RESULT = 3
end

@enum CacheResultStatus::Int8 begin
    RESULT_READY = 1
    RESULT_FAILED = 2
end

struct CacheResultKey
    kind::CacheResultKind
    entity::String
end

Base.:(==)(left::CacheResultKey, right::CacheResultKey)::Bool =
    left.kind == right.kind && left.entity == right.entity
Base.isequal(left::CacheResultKey, right::CacheResultKey)::Bool =
    left.kind == right.kind && isequal(left.entity, right.entity)
Base.hash(key::CacheResultKey, seed::UInt)::UInt = hash(key.entity, hash(key.kind, seed))

"""One persisted result-state row. `CacheResultKey` is rebuilt at the domain edge from `(kind, entity)`."""
struct CachedResultState
    kind::Int8
    entity::String
    status::Int8
    source_item_id::String
    message::Union{Nothing,String}
end

CachedResultState(row)::CachedResultState = CachedResultState(
    Int8(row.kind),
    String(row.entity),
    Int8(row.status),
    String(row.source_item_id),
    _null_to_nothing(row.message),
)

"""All data-less content needed to restore and compare a project cache."""
struct ProjectCacheIndex
    identity::ProjectCacheIdentity
    source::SourceScan
    analysis_errors::Dict{String,String}
    result_states::Dict{CacheResultKey,CachedResultState}
end

"""Summary of the cache/source comparison shown by workspace status."""
struct ProjectCacheStatus
    total_source_items::Int
    cached_source_items::Int
    fresh_source_items::Int
    stale_source_items::Int
    new_source_items::Int
    deleted_source_items::Int
    error_source_items::Int
end

# Storage rows (field names are the table columns; `R(database_row)` is the pure decoder)

struct SourceItemRow
    id::String
    fingerprint_hex::Union{Nothing,String}
    fingerprint_hash::Union{Nothing,String}
    path::Union{Nothing,String}
    timestamp::Union{Nothing,DateTime}
end

SourceItemRow(row)::SourceItemRow = SourceItemRow(
    String(row.id),
    _null_to_nothing(row.fingerprint_hex),
    _null_to_nothing(row.fingerprint_hash),
    _null_to_nothing(row.path),
    _null_to_nothing(row.timestamp),
)

struct ItemRow
    id::String
    source_item_id::String
    item_label::String
    kind::String
    collection::Vector{String}
    item_fingerprint_hex::Union{Nothing,String}
end

ItemRow(row)::ItemRow = ItemRow(
    String(row.id),
    String(row.source_item_id),
    String(row.item_label),
    String(row.kind),
    Vector{String}(row.collection),
    _null_to_nothing(row.item_fingerprint_hex),
)

struct MetaRow
    scope::Int8
    entity::String
    key::String
    vtype::Int8
    b::Union{Missing,Bool}
    i::Union{Missing,Int64}
    d::Union{Missing,Float64}
    s::Union{Missing,String}
    ts::Union{Missing,DateTime}
    lb::Union{Missing,Vector{Bool}}
    li::Union{Missing,Vector{Int64}}
    ld::Union{Missing,Vector{Float64}}
    ls::Union{Missing,Vector{String}}
end

MetaRow(row)::MetaRow = MetaRow(
    Int8(row.scope),
    String(row.entity),
    String(row.key),
    Int8(row.vtype),
    row.b,
    row.i,
    row.d,
    row.s,
    row.ts,
    _meta_column(Bool, row.lb),
    _meta_column(Int64, row.li),
    _meta_column(Float64, row.ld),
    _meta_column(String, row.ls),
)

struct FailureRow
    item_id::String
    source_item_id::String
    message::String
end

FailureRow(row)::FailureRow =
    FailureRow(String(row.item_id), String(row.source_item_id), String(row.message))

# Identity

"""Bind one project name to a source and its package-owned cache path. The name is used verbatim as one directory component; invalid names fail rather than being silently rewritten."""
function project_cache_identity(
    project_name::AbstractString,
    source,
)::ProjectCacheIdentity
    name = String(project_name)
    isempty(name) && throw(ArgumentError("project name cannot be empty"))
    name in (".", "..") && throw(ArgumentError("project name '$name' is not a safe directory name"))
    (isabspath(name) || basename(name) != name || occursin('\\', name) || occursin('\0', name)) &&
        throw(ArgumentError(
            "project name '$name' must be one directory name without path separators",
        ))
    depot = isempty(DEPOT_PATH) ? homedir() : DEPOT_PATH[1]
    return ProjectCacheIdentity(
        name,
        source_id(source),
        source_label(source),
        joinpath(depot, "measurementbrowser", name, "cache.duckdb"),
    )
end

# Connection + schema

"""One open DuckDB cache file for a workspace: typed `RowStore`/`TabularFamilyStore` owners plus memory-only payload buffers. Only schema, meta-header, and checkpoint work open short-lived connections."""
mutable struct CacheDB
    identity::ProjectCacheIdentity
    db::DuckDB.DB
    metrics::BuildMetrics
    source_items::RowStore{String,SourceItemRow}
    items::RowStore{String,ItemRow}
    metadata::RowStore{Tuple{Int8,String,String},MetaRow}
    failures::RowStore{Tuple{String,String},FailureRow}
    result_states::RowStore{Tuple{Int8,String},CachedResultState}
    payload::TabularFamilyStore
    interpreted::MemoryStore{String,AbstractDataItem}
    processed_memory::MemoryStore{String,AbstractDataItem}
    profiler::Profiling.ProfileSession
end

"""Delete one generated cache file. Call only from an explicit rebuild path."""
function _remove_cache_files!(path::AbstractString)::Nothing
    cache_path = String(path)
    for candidate in (cache_path, cache_path * ".wal")
        ispath(candidate) && rm(candidate; force=true)
    end
    return nothing
end

"""Open one generated cache. `rebuild=true` first discards the old generated file."""
function open_cache_db(
    identity::ProjectCacheIdentity,
    profiler::Profiling.ProfileSession,
    metrics::BuildMetrics=BuildMetrics(),
    ;
    rebuild::Bool=false,
)::CacheDB
    mkpath(dirname(identity.cache_path))
    rebuild && _remove_cache_files!(identity.cache_path)
    db = DBInterface.connect(DuckDB.DB, identity.cache_path)
    connection = DBInterface.connect(db)
    try
        try
            DBInterface.execute(connection, _memory_limit_sql(CACHE_MEMORY_LIMIT_MIB[]))
            has_meta = only(Int[
                row.n
                for row in DBInterface.execute(
                    connection,
                    "SELECT COUNT(*) AS n FROM information_schema.tables " *
                    "WHERE table_schema = 'main' AND table_name = 'meta'",
                )
            ]) > 0
            if has_meta
                schema_versions = String[
                    String(row.value)
                    for row in DBInterface.execute(
                        connection,
                        "SELECT value FROM meta WHERE key = 'schema_version'",
                    )
                ]
                if !isempty(schema_versions) &&
                   only(schema_versions) != string(PROJECT_CACHE_SCHEMA_VERSION)
                    throw(ProjectCacheError(
                        identity.cache_path,
                        "cache schema is out of date; pass rebuild=true to discard and rebuild it",
                    ))
                end
            end
            ensure_schema!(connection)
        finally
            DBInterface.close!(connection)
        end
        payload = TabularFamilyStore(db; profiler, row_limit=CACHE_BUFFER_ROW_LIMIT)
        return CacheDB(
            identity,
            db,
            metrics,
            RowStore{String,SourceItemRow}(db, "source_items", ("id",); profiler),
            RowStore{String,ItemRow}(db, "items", ("id",); profiler),
            RowStore{Tuple{Int8,String,String},MetaRow}(
                db, "metadata", ("scope", "entity", "key"); profiler),
            RowStore{Tuple{String,String},FailureRow}(
                db, "item_failures", ("item_id", "source_item_id"); profiler),
            RowStore{Tuple{Int8,String},CachedResultState}(
                db, "result_states", ("kind", "entity"); profiler),
            payload,
            MemoryStore{String,AbstractDataItem}(; row_limit=CACHE_BUFFER_ROW_LIMIT),
            MemoryStore{String,AbstractDataItem}(; row_limit=CACHE_BUFFER_ROW_LIMIT),
            profiler,
        )
    catch
        DBInterface.close!(db)
        rethrow()
    end
end

"""Open a cache with internal profiling disabled for direct cache operations and tests."""
function open_cache_db(identity::ProjectCacheIdentity; rebuild::Bool=false)::CacheDB
    return open_cache_db(
        identity,
        Profiling.ProfileSession(false, false, nothing, nothing);
        rebuild,
    )
end

"""Cache stores are live when opened; this keeps Workspace on a semantic lifecycle call."""
start_cache!(::CacheDB)::Nothing = nothing

"""Cache stores close with `close_cache_db!`; this keeps Workspace on a semantic lifecycle call."""
stop_cache!(::CacheDB)::Nothing = nothing

function _buffer_rows(item::AbstractDataItem)::Int
    data = item_data(item)
    return data isa AbstractDataFrame ? nrow(data) : 1
end

"""Create the cache tables if they do not already exist."""
function ensure_schema!(connection)::Nothing
    DBInterface.execute(connection, """
        CREATE TABLE IF NOT EXISTS meta(
            key TEXT PRIMARY KEY,
            value TEXT)
    """)
    DBInterface.execute(connection, create_table_sql(SourceItemRow, "source_items", (:id,)))
    DBInterface.execute(connection, create_table_sql(ItemRow, "items", (:id,)))
    DBInterface.execute(connection, create_table_sql(MetaRow, "metadata", (:scope, :entity, :key)))
    DBInterface.execute(connection, create_table_sql(FailureRow, "item_failures", (:item_id, :source_item_id)))
    DBInterface.execute(connection, create_table_sql(CachedResultState, "result_states", (:kind, :entity)))
    DBInterface.execute(connection, """
        CREATE TABLE IF NOT EXISTS item_data(
            item_id TEXT PRIMARY KEY,
            storage_id USMALLINT,
            seq UINTEGER,
            UNIQUE (storage_id, seq))
    """)
    DBInterface.execute(connection, """
        CREATE TABLE IF NOT EXISTS dataframe_schemas(
            storage_id USMALLINT PRIMARY KEY,
            column_names VARCHAR[],
            column_types VARCHAR[])
    """)
    return nothing
end

# The domain write/read surface over the generic stores: pipeline stages and `ItemRecord`s map onto
# RowStore rows and TabularFamilyStore batches. The stores know none of this; they buffer typed rows.

"""Store one source item's interpreted result: data-less records to disk-backed buffers, payloads to memory. Item statistics are published independently."""
function store_interpreted!(cache::CacheDB, records::Vector{ItemRecord},
        data::Vector{<:AbstractDataItem})::Nothing
    length(records) == length(data) || throw(ArgumentError(
        "Cannot store interpreted data: received $(length(records)) records and " *
        "$(length(data)) items; the vectors must align one-to-one",
    ))
    started = time_ns()
    for record in records
        hex, hash = _encode_fingerprint(record.source_item_fingerprint)
        append!(cache.source_items, record.source_item_id,
            SourceItemRow(record.source_item_id, hex, hash,
                record.source_item_path, record.source_item_timestamp))
        append!(cache.items, record.id,
            ItemRow(record.id, record.source_item_id, record.item_label, String(record.kind),
                record.collection, _serialize_hex(record.item_fingerprint)))
        for (key, value) in record.parameters
            append!(cache.metadata, (Int8(SCOPE_ITEM_PARAMETERS), record.id, String(key)),
                metadata_value_to_row(Int8(SCOPE_ITEM_PARAMETERS), record.id, String(key), value))
        end
    end
    for (record, item) in zip(records, data)
        append!(cache.interpreted, record.id, item)
    end
    record_cache_phase!(
        cache.metrics.interpreted_write_ns,
        cache.metrics.interpreted_writes,
        started,
    )
    return nothing
end

"""Store one processed result in the payload family when cacheable, else in the memory-only processed buffer. Non-cacheable items are held like interpreted ones (docs/cache.md "Two backing modes")."""
function store_processed!(
    cache::CacheDB,
    record::ItemRecord,
    item::AbstractDataItem,
)::Nothing
    payload = item_data(item)
    disk = cacheable(item) && payload isa AbstractDataFrame && !isempty(names(payload))
    if disk
        started = time_ns()
        delete!(cache.processed_memory, record.id)
        append!(cache.payload, record.id, payload)
        record_cache_phase!(
            cache.metrics.processed_write_ns,
            cache.metrics.processed_writes,
            started,
        )
    else
        delete!(cache.payload, record.id)
        append!(cache.processed_memory, record.id, item)
    end
    delete!(cache.result_states, (Int8(PROCESSING_RESULT), record.id))
    return nothing
end

"""Persist one independent failed work result."""
function store_result_failure!(cache::CacheDB, key::CacheResultKey,
        source_item_id::AbstractString, message::AbstractString)::Nothing
    edit!(cache.result_states, (Int8(key.kind), key.entity),
        CachedResultState(Int8(key.kind), key.entity, Int8(RESULT_FAILED),
            String(source_item_id), String(message)))
    return nothing
end

"""Persist one item-stat result, including a successful empty result."""
function store_item_stats!(cache::CacheDB, record::ItemRecord, stats::AbstractDict)::Nothing
    started = time_ns()
    for (key, value) in metadata_dict(stats)
        edit!(cache.metadata, (Int8(SCOPE_ITEM_STATS), record.id, String(key)),
            metadata_value_to_row(Int8(SCOPE_ITEM_STATS), record.id, String(key), value))
    end
    result_key = CacheResultKey(ITEM_STATS_RESULT, record.id)
    edit!(cache.result_states, (Int8(result_key.kind), result_key.entity),
        CachedResultState(Int8(result_key.kind), result_key.entity, Int8(RESULT_READY),
            record.source_item_id, nothing))
    record_cache_phase!(
        cache.metrics.stats_write_ns,
        cache.metrics.stats_writes,
        started,
    )
    return nothing
end

"""Persist one collection-stat result, including a successful empty result."""
function store_collection_stats!(
    cache::CacheDB, collection_key::AbstractString, stats::AbstractDict,
)::Nothing
    entity = String(collection_key)
    for (key, value) in metadata_dict(stats)
        edit!(cache.metadata, (Int8(SCOPE_NODE_STATS), entity, String(key)),
            metadata_value_to_row(Int8(SCOPE_NODE_STATS), entity, String(key), value))
    end
    result_key = CacheResultKey(COLLECTION_STATS_RESULT, entity)
    edit!(cache.result_states, (Int8(result_key.kind), result_key.entity),
        CachedResultState(Int8(result_key.kind), result_key.entity, Int8(RESULT_READY), "", nothing))
    return nothing
end

"""Delete every cached descendant of one source item through typed buffer mutations."""
function delete_source_item!(
    cache::CacheDB,
    source_item_id::AbstractString,
    item_ids::Vector{String},
    metadata_keys::Vector{Tuple{Int8,String,String}},
    failure_keys::Vector{Tuple{String,String}},
    result_keys::Vector{Tuple{Int8,String}},
)::Nothing
    source_id = String(source_item_id)
    delete!(cache.source_items, source_id)
    for item_id in item_ids
        delete!(cache.items, item_id)
        delete!(cache.payload, item_id)
        delete!(cache.interpreted, item_id)
        delete!(cache.processed_memory, item_id)
    end
    for key in metadata_keys
        delete!(cache.metadata, key)
    end
    for key in failure_keys
        delete!(cache.failures, key)
    end
    for key in result_keys
        delete!(cache.result_states, key)
    end
    return nothing
end

"""Delete every cached result owned by one published source item."""
function delete_source_item!(
    cache::CacheDB,
    source_item_id::AbstractString,
    old_records::Vector{ItemRecord},
)::Nothing
    source_id = String(source_item_id)
    item_ids = String[record.id for record in old_records]
    metadata_keys = Tuple{Int8,String,String}[]
    failure_keys = Tuple{String,String}[]
    result_keys = Tuple{Int8,String}[]
    for record in old_records
        for key in keys(record.parameters)
            push!(metadata_keys, (Int8(SCOPE_ITEM_PARAMETERS), record.id, String(key)))
        end
        for key in keys(record.stats)
            push!(metadata_keys, (Int8(SCOPE_ITEM_STATS), record.id, String(key)))
        end
        push!(failure_keys, (record.id, source_id))
        push!(result_keys, (Int8(PROCESSING_RESULT), record.id))
        push!(result_keys, (Int8(ITEM_STATS_RESULT), record.id))
    end
    push!(failure_keys, (source_id, source_id))
    return delete_source_item!(
        cache, source_id, item_ids, metadata_keys, failure_keys, result_keys)
end

"""Delete cached statistics for collections whose published stats are invalid."""
function delete_collection_stats!(
    cache::CacheDB,
    collection_keys::Vector{String},
)::Nothing
    for collection_key in unique(collection_keys)
        delete_by_key_prefix!(cache.metadata, (Int8(SCOPE_NODE_STATS), collection_key))
        delete!(cache.result_states, (Int8(COLLECTION_STATS_RESULT), collection_key))
    end
    return nothing
end

"""Aggregate staged durable keys and payload rows."""
function cache_pending_counts(cache::CacheDB)::NamedTuple{(:items, :rows),Tuple{Int,Int}}
    items = 0
    rows = 0
    for buffer in
        (cache.source_items, cache.items, cache.metadata, cache.failures,
         cache.result_states, cache.payload)
        lock(buffer.flush_condition) do
            items += length(buffer.queued) + length(buffer.writing)
            rows += buffer.queued_rows + buffer.writing_rows
        end
    end
    return (items=items, rows=rows)
end

"""Whether any cache store has queued or writing mutations."""
function cache_has_pending_writes(cache::CacheDB)::Bool
    pending = cache_pending_counts(cache)
    return pending.items > 0 || pending.rows > 0
end

"""Read item data through the buffer (memory first, then the database for processed). Returns a vector aligned to `records` of loaded items or `nothing` on a miss; the caller falls back to the source. Interpreted payloads are memory-only, so a memory miss is final."""
function read_item_data(cache::CacheDB, records::Vector{ItemRecord}; stage::Symbol)::Vector{Any}
    stage in (:interpreted, :processed) ||
        throw(ArgumentError("unknown item-data cache stage '$stage'"))
    isempty(records) && return Any[]

    if stage === :interpreted
        return Any[read(cache.interpreted, record.id) for record in records]
    end

    loaded = Vector{Any}(undef, length(records))
    disk_ids = String[]
    for (index, record) in pairs(records)
        held = read(cache.processed_memory, record.id)
        if held !== nothing
            loaded[index] = held
            continue
        end
        push!(disk_ids, record.id)
        loaded[index] = nothing
    end

    disk_data = isempty(disk_ids) ?
        Dict{String,Union{Nothing,AbstractDataFrame}}() :
        read(cache.payload, unique(disk_ids))
    for (index, record) in pairs(records)
        loaded[index] !== nothing && continue
        data = get(disk_data, record.id, nothing)
        loaded[index] = data === nothing ? nothing : DataItem(record, data)
    end
    return loaded
end

"""Close a workspace cache's database file, checkpointing first."""
function close_cache_db!(cachedb::CacheDB)::Nothing
    failure::Union{Nothing,Exception} = nothing
    for buffer in
        (cachedb.source_items, cachedb.items, cachedb.metadata, cachedb.failures,
         cachedb.result_states, cachedb.payload, cachedb.interpreted, cachedb.processed_memory)
        try
            close!(buffer)
        catch error
            failure === nothing && (failure = error)
        end
    end
    connection = DBInterface.connect(cachedb.db)
    try
        DBInterface.execute(connection, "CHECKPOINT")
    finally
        DBInterface.close!(connection)
        DBInterface.close!(cachedb.db)
    end
    failure === nothing || throw(failure)
    return nothing
end

# Encoding helpers

"""Serialize any value to bytes — the basis for both content hashing and hex-text storage."""
function _serialize_bytes(value)::Vector{UInt8}
    io = IOBuffer()
    serialize(io, value)
    return take!(io)
end

"""Hex-encode any value for a TEXT column, or `nothing` (SQL NULL) when absent."""
_serialize_hex(value) = value === nothing ? nothing : bytes2hex(_serialize_bytes(value))

"""Deserialize a hex-text cell back to its value (or `nothing` for SQL NULL)."""
function _deserialize_hex(value)
    (value === nothing || ismissing(value)) && return nothing
    return deserialize(IOBuffer(hex2bytes(String(value))))
end

"""Map a SQL NULL (`missing`) to `nothing`, passing other values through."""
_null_to_nothing(value) = ismissing(value) ? nothing : value

"""Coerce one metadata list column to its stored element type, preserving SQL NULL as `missing`."""
_meta_column(::Type{T}, value) where {T} = ismissing(value) ? missing : Vector{T}(value)

"""Hex-text storage and content hash for one fingerprint, as `(hex, hash)` (`nothing` when absent)."""
function _encode_fingerprint(fingerprint)::Tuple{Union{Nothing,String},Union{Nothing,String}}
    fingerprint === nothing && return (nothing, nothing)
    bytes = _serialize_bytes(fingerprint)
    return (bytes2hex(bytes), bytes2hex(sha1(bytes)))
end

"""Encode one [`MetadataValue`](@ref) into a flat [`MetaRow`](@ref): the active typed column is set, the rest stay SQL NULL. Decode with [`metadata_row_to_value`](@ref)."""
function metadata_value_to_row(scope::Int8, entity::String, key::String, value)::MetaRow
    b = missing; i = missing; d = missing; s = missing; ts = missing
    lb = missing; li = missing; ld = missing; ls = missing
    vtype = if value isa Bool
        b = value; VT_BOOL
    elseif value isa Int64
        i = value; VT_INT
    elseif value isa Float64
        d = value; VT_FLOAT
    elseif value isa String
        s = value; VT_STRING
    elseif value isa Symbol
        s = String(value); VT_SYMBOL
    elseif value isa Date
        ts = DateTime(value); VT_DATE
    elseif value isa DateTime
        ts = value; VT_DATETIME
    elseif value === missing
        VT_MISSING
    elseif value isa Vector{Bool}
        lb = value; VT_VBOOL
    elseif value isa Vector{Int64}
        li = value; VT_VINT
    elseif value isa Vector{Float64}
        ld = value; VT_VFLOAT
    elseif value isa Vector{String}
        ls = value; VT_VSTRING
    else
        error("No metadata encoding for $(typeof(value))")
    end
    return MetaRow(scope, entity, key, Int8(vtype), b, i, d, s, ts, lb, li, ld, ls)
end

"""Reconstruct one [`MetadataValue`](@ref) from a flat [`MetaRow`](@ref) at the domain edge."""
function metadata_row_to_value(row::MetaRow)::MetadataValue
    vtype = MetaVType(row.vtype)
    vtype === VT_BOOL && return row.b
    vtype === VT_INT && return row.i
    vtype === VT_FLOAT && return row.d
    vtype === VT_STRING && return row.s
    vtype === VT_SYMBOL && return Symbol(row.s)
    vtype === VT_DATE && return Date(row.ts)
    vtype === VT_DATETIME && return row.ts
    vtype === VT_MISSING && return missing
    vtype === VT_VBOOL && return Vector{Bool}(row.lb)
    vtype === VT_VINT && return Vector{Int64}(row.li)
    vtype === VT_VFLOAT && return Vector{Float64}(row.ld)
    vtype === VT_VSTRING && return Vector{String}(row.ls)
    error("Unknown metadata vtype $vtype")
end

# Index construction + comparison (in memory; unchanged from the previous cache)

"""Build the complete cache index for one finished source scan. Analysis failures stay attached to their source items; successful items stay available."""
function ProjectCacheIndex(
    identity::ProjectCacheIdentity,
    source::SourceScan,
)::ProjectCacheIndex
    errors = Dict{String,Vector{String}}()
    for failure in source.analysis_failures
        owner = isempty(failure.id) ? failure.source_item_id : failure.id
        push!(get!(() -> String[], errors, owner), failure.message)
    end
    return ProjectCacheIndex(
        identity,
        source,
        Dict(id => join(messages, "\n") for (id, messages) in errors),
        Dict{CacheResultKey,CachedResultState}(),
    )
end

# Incremental writing (per source item, as the scan streams)

"""Delete every index and item-data row, leaving an empty (but schema-valid) cache for a rebuild."""
function clear_cache_index!(cachedb::CacheDB)::Nothing
    for buffer in
        (cachedb.source_items, cachedb.items, cachedb.metadata, cachedb.failures,
         cachedb.result_states, cachedb.payload, cachedb.interpreted, cachedb.processed_memory)
        clear!(buffer)
    end
    connection = DBInterface.connect(cachedb.db)
    try
        DBInterface.execute(connection, "DELETE FROM meta")
    finally
        DBInterface.close!(connection)
    end
    return nothing
end

"""Stamp the identity `meta` rows so a scan's incremental writes are loadable before it finishes. Without them [`load_cache_index`](@ref) treats the cache as unbuilt. Scan-wide summaries are derived on load, not stored."""
function write_meta_header!(cachedb::CacheDB)::Nothing
    identity = cachedb.identity
    connection = DBInterface.connect(cachedb.db)
    try
        statement = DBInterface.prepare(connection,
            "INSERT INTO meta VALUES (?, ?) ON CONFLICT (key) DO UPDATE SET value = excluded.value")
        for (key, value) in (
            "schema_version" => string(PROJECT_CACHE_SCHEMA_VERSION),
            "project_name" => identity.project_name,
            "source_id" => identity.source_id,
            "source_label" => identity.source_label,
        )
            DBInterface.execute(statement, (key, value))
        end
    finally
        DBInterface.close!(connection)
    end
    return nothing
end

# Loading the index

"""Read and validate the sole unbuffered cache header."""
function _load_meta(cache::CacheDB)::Dict{String,String}
    identity = cache.identity
    meta = Dict{String,String}()
    connection = DBInterface.connect(cache.db)
    try
        for row in DBInterface.execute(connection, "SELECT key, value FROM meta")
            meta[String(row.key)] = String(row.value)
        end
    finally
        DBInterface.close!(connection)
    end
    isempty(meta) &&
        throw(ProjectCacheError(identity.cache_path, "cache has not been built yet"))
    get(meta, "schema_version", "") == string(PROJECT_CACHE_SCHEMA_VERSION) ||
        throw(ProjectCacheError(
            identity.cache_path,
            "cache schema is out of date; pass rebuild=true to discard and rebuild it",
        ))
    cached_project = get(meta, "project_name", "")
    cached_project == identity.project_name || throw(ProjectCacheError(
        identity.cache_path,
        "cache belongs to project '$cached_project', not '$(identity.project_name)'",
    ))
    cached_source = get(meta, "source_id", "")
    cached_source == identity.source_id || throw(ProjectCacheError(
        identity.cache_path,
        "project '$(identity.project_name)' is cached for source '$cached_source', not " *
        "'$(identity.source_id)'; use Rebuild Cache to replace it",
    ))
    return meta
end

"""Whether the cache has an identity header; full identity validation happens on load."""
function cache_built(cache::CacheDB)::Bool
    connection = DBInterface.connect(cache.db)
    try
        count = only(DBInterface.execute(
            connection,
            "SELECT count(*) AS count FROM meta",
        )).count
        return count > 0
    finally
        DBInterface.close!(connection)
    end
end

"""Group the metadata buffer's complete logical view."""
function _load_metadata(cache::CacheDB)::Dict{Tuple{MetaScope,String},MetadataDict}
    grouped = Dict{Tuple{MetaScope,String},MetadataDict}()
    for row in values(read(cache.metadata))
        scope = MetaScope(row.scope)
        dict = get!(() -> MetadataDict(), grouped, (scope, row.entity))
        dict[Symbol(row.key)] = metadata_row_to_value(row)
    end
    return grouped
end

"""Reconstruct the full `SourceScan` stored in one DuckDB cache."""
function _load_source_scan(cache::CacheDB)::SourceScan
    identity = cache.identity
    _load_meta(cache)
    source_rows = read(cache.source_items)
    item_rows = read(cache.items)
    failure_rows = read(cache.failures)
    metadata = _load_metadata(cache)

    all_source_ids = String[String(id) for id in keys(source_rows)]
    fingerprints = Dict{String,Any}()
    locations = Dict{String,Tuple{Union{Nothing,String},Union{Nothing,DateTime}}}()
    failures = ItemFailure[]
    for row in values(source_rows)
        fingerprints[row.id] = _deserialize_hex(row.fingerprint_hex)
        locations[row.id] = (row.path, row.timestamp)
    end
    for row in values(failure_rows)
        push!(failures, ItemFailure(row.source_item_id, row.item_id, row.message))
    end

    records = ItemRecord[]
    for row in values(item_rows)
        path, timestamp = get(locations, row.source_item_id, (nothing, nothing))
        push!(records, ItemRecord(;
            id=row.id,
            source_item_id=row.source_item_id,
            source_item_fingerprint=get(fingerprints, row.source_item_id, nothing),
            source_item_path=path,
            source_item_timestamp=timestamp,
            item_label=row.item_label,
            kind=Symbol(row.kind),
            collection=row.collection,
            parameters=get(metadata, (SCOPE_ITEM_PARAMETERS, row.id), MetadataDict()),
            stats=get(metadata, (SCOPE_ITEM_STATS, row.id), MetadataDict()),
            item_fingerprint=_deserialize_hex(row.item_fingerprint_hex),
        ))
    end

    sort!(records; by=record -> (record.source_item_id, record.id))
    # Scan-wide summaries are derived, not stored: a source item with no items was skipped or failed
    # (both produce zero records), and collection parameters exist iff any node carries them.
    source_ids_with_items = Set(record.source_item_id for record in records)
    skipped_count = count(id -> !(id in source_ids_with_items), all_source_ids)
    has_collection_parameters = any(key -> key[1] == SCOPE_NODE_PARAMETERS, keys(metadata))
    hierarchy = Hierarchy(identity.source_id, has_collection_parameters, skipped_count)
    for record in records
        insert_item!(hierarchy, record)
    end
    for ((scope, entity), dict) in metadata
        (scope == SCOPE_NODE_PARAMETERS || scope == SCOPE_NODE_STATS) || continue
        node = get(hierarchy.index, collection_path_tuple(entity), nothing)
        node === nothing && continue
        target = scope == SCOPE_NODE_PARAMETERS ? node.parameters : node.stats
        merge!(target, dict)
    end
    sort!(hierarchy)

    return SourceScan(
        identity.source_id,
        identity.source_label,
        fingerprints,
        hierarchy,
        failures,
    )
end

"""Load the complete browser index from a workspace's open cache, snapshotting committed state. Reuses the already-open [`CacheDB`](@ref) instead of opening a second handle."""
function load_cache_index(
    cachedb::CacheDB;
    on_progress::Union{Nothing,Function}=nothing,
)::ProjectCacheIndex
    identity = cachedb.identity
    index = @profile_span cachedb.profiler :cache :load_index ProfileAttributes(
        source_id=identity.source_id,
    ) begin
        base = ProjectCacheIndex(identity, _load_source_scan(cachedb))
        records = Dict(record.id => record for record in base.source.hierarchy.all_items)
        result_states = Dict{CacheResultKey,CachedResultState}()
        for item_id in cached_item_ids(cachedb.payload)
            record = get(records, item_id, nothing)
            record === nothing && continue
            key = CacheResultKey(PROCESSING_RESULT, item_id)
            result_states[key] =
                CachedResultState(Int8(PROCESSING_RESULT), item_id, Int8(RESULT_READY),
                    record.source_item_id, nothing)
        end
        for state in values(read(cachedb.result_states))
            result_states[CacheResultKey(CacheResultKind(state.kind), state.entity)] = state
        end
        ProjectCacheIndex(base.identity, base.source, base.analysis_errors, result_states)
    end
    emit_progress(
        on_progress;
        phase=:cache_load,
        total_source_items=1,
        processed_source_items=1,
        loaded_items=length(index.source.hierarchy.all_items),
        skipped_source_items=index.source.hierarchy.skipped_count,
        current_source_item=index.identity.cache_path,
    )
    return index
end
