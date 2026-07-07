const PROJECT_CACHE_SCHEMA_VERSION = 12

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

"""Error raised when a project cache is missing or does not match the current project."""
struct ProjectCacheError <: Exception
    path::String
    message::String
end

"""Print the cache path and validation failure."""
Base.showerror(io::IO, err::ProjectCacheError)::Nothing =
    print(io, "Invalid project cache $(err.path): $(err.message)")

"""Error raised when the generated cache file exists but uses an old schema."""
struct ProjectCacheSchemaError <: Exception
    path::String
    message::String
end

Base.showerror(io::IO, err::ProjectCacheSchemaError)::Nothing =
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
    ITEM_ANALYSIS_RESULT = 2
    COLLECTION_ANALYSIS_RESULT = 3
    COLLECTION_PROCESS_RESULT = 4
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
    # The computed layers (analyze output + collection-process overwrite) per item, restored for
    # items whose governing result is valid: delivered effective minus the entries layer.
    item_metadata::Dict{String,MetadataDict}
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

"""Persisted cache counts surfaced through workspace status."""
struct CacheStageSummary
    cached_sources::Int
    interpreted_items::Int
    processed::Int
    analyzed::Int
    collection_processed::Int
    collection_analyzed::Int
    failed_interpret::Int
    failed_process::Int
    failed_analyze::Int
    failed_collection::Int
end

CacheStageSummary()::CacheStageSummary =
    CacheStageSummary(0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

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
    item_key::Int64
    source_item_id::String
    item_label::String
    kind::String
    collection::Vector{String}
    item_fingerprint_hex::Union{Nothing,String}
end

ItemRow(row)::ItemRow = ItemRow(
    String(row.id),
    Int64(row.item_key),
    String(row.source_item_id),
    String(row.item_label),
    String(row.kind),
    Vector{String}(row.collection),
    _null_to_nothing(row.item_fingerprint_hex),
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

abstract type AbstractCacheDB end

"""One open DuckDB cache file for a workspace: typed `RowStore`/`TabularFamilyStore` owners plus memory-only payload buffers. Only schema, meta-header, and checkpoint work open short-lived connections."""
mutable struct CacheDB <: AbstractCacheDB
    identity::ProjectCacheIdentity
    db::DuckDB.DB
    metrics::BuildMetrics
    source_items::RowStore{String,SourceItemRow}
    items::RowStore{String,ItemRow}
    source_item_metadata::WideRowStore{Int64}
    source_collection_metadata::WideRowStore{String}
    analyzed_item_metadata::WideRowStore{Int64}
    analyzed_collection_metadata::WideRowStore{String}
    failures::RowStore{Tuple{String,String},FailureRow}
    result_states::RowStore{Tuple{Int8,String},CachedResultState}
    payload::TabularFamilyStore
    interpreted::MemoryStore{String,AbstractDataItem}
    processed_memory::MemoryStore{String,AbstractDataItem}
    # One integer surrogate per item id, shared with the payload store; minted at interpretation.
    item_keys::Dict{String,Int64}
    next_item_key::Int64
    key_lock::ReentrantLock
    # The signature of the analyzed_item_metadata columns backing the query view; the view is
    # recreated only when it changes.
    query_view_signature::Vector{Symbol}
    profiler::Profiling.ProfileSession
end

"""Session-only cache state used when DuckDB persistence is disabled or unavailable."""
mutable struct MemoryCacheDB <: AbstractCacheDB
    identity::ProjectCacheIdentity
    metrics::BuildMetrics
    interpreted::MemoryStore{String,AbstractDataItem}
    processed_memory::MemoryStore{String,AbstractDataItem}
    collection_processed_memory::MemoryStore{String,AbstractDataItem}
    source_items::Set{String}
    items::Set{String}
    failures::Dict{Tuple{String,String},String}
    result_states::Dict{Tuple{Int8,String},CachedResultState}
    lock::ReentrantLock
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
                    throw(ProjectCacheSchemaError(
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
        item_keys, next_item_key = _load_item_keys(db)
        return CacheDB(
            identity,
            db,
            metrics,
            RowStore{String,SourceItemRow}(db, "source_items", ("id",); profiler),
            RowStore{String,ItemRow}(db, "items", ("id",); profiler),
            WideRowStore{Int64}(db, "source_item_metadata", "item_key"; profiler),
            WideRowStore{String}(db, "source_collection_metadata", "collection_key"; profiler),
            WideRowStore{Int64}(db, "analyzed_item_metadata", "item_key"; profiler),
            WideRowStore{String}(db, "analyzed_collection_metadata", "collection_key"; profiler),
            RowStore{Tuple{String,String},FailureRow}(
                db, "item_failures", ("item_id", "source_item_id"); profiler),
            RowStore{Tuple{Int8,String},CachedResultState}(
                db, "result_states", ("kind", "entity"); profiler),
            payload,
            MemoryStore{String,AbstractDataItem}(; row_limit=CACHE_BUFFER_ROW_LIMIT),
            MemoryStore{String,AbstractDataItem}(; row_limit=CACHE_BUFFER_ROW_LIMIT),
            item_keys,
            next_item_key,
            ReentrantLock(),
            Symbol[],
            profiler,
        )
    catch
        DBInterface.close!(db)
        rethrow()
    end
end

"""Load the committed item-id → surrogate map and the next free surrogate from the `items` table."""
function _load_item_keys(db::DuckDB.DB)::Tuple{Dict{String,Int64},Int64}
    connection = DBInterface.connect(db)
    try
        item_keys = Dict{String,Int64}()
        for row in DBInterface.execute(connection, "SELECT id, item_key FROM items")
            item_keys[String(row.id)] = Int64(row.item_key)
        end
        next = isempty(item_keys) ? Int64(1) : maximum(values(item_keys)) + 1
        return item_keys, next
    finally
        DBInterface.close!(connection)
    end
end

"""Mint or resolve the integer surrogate for one item id. Minting happens only at interpretation."""
function item_key!(cache::CacheDB, id::AbstractString; mint::Bool=false)::Int64
    key = String(id)
    return lock(cache.key_lock) do
        existing = get(cache.item_keys, key, nothing)
        existing !== nothing && return existing
        mint || error("no item_key minted for item '$key'")
        minted = cache.next_item_key
        cache.next_item_key += 1
        cache.item_keys[key] = minted
        minted
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

function open_memory_cache_db(
    identity::ProjectCacheIdentity,
    profiler::Profiling.ProfileSession,
    metrics::BuildMetrics=BuildMetrics(),
)::MemoryCacheDB
    return MemoryCacheDB(
        identity,
        metrics,
        MemoryStore{String,AbstractDataItem}(; row_limit=CACHE_BUFFER_ROW_LIMIT),
        MemoryStore{String,AbstractDataItem}(; row_limit=CACHE_BUFFER_ROW_LIMIT),
        MemoryStore{String,AbstractDataItem}(; row_limit=CACHE_BUFFER_ROW_LIMIT),
        Set{String}(),
        Set{String}(),
        Dict{Tuple{String,String},String}(),
        Dict{Tuple{Int8,String},CachedResultState}(),
        ReentrantLock(),
        profiler,
    )
end

"""Cache stores are live when opened; this keeps Workspace on a semantic lifecycle call."""
start_cache!(::AbstractCacheDB)::Nothing = nothing

"""Cache stores close with `close_cache_db!`; this keeps Workspace on a semantic lifecycle call."""
stop_cache!(::CacheDB)::Nothing = nothing
stop_cache!(::MemoryCacheDB)::Nothing = nothing

set_cache_memory_limit!(::AbstractCacheDB, mib::Integer)::Int =
    set_cache_memory_limit!(mib)

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
    DBInterface.execute(connection, create_table_sql(FailureRow, "item_failures", (:item_id, :source_item_id)))
    DBInterface.execute(connection, create_table_sql(CachedResultState, "result_states", (:kind, :entity)))
    DBInterface.execute(connection, """
        CREATE TABLE IF NOT EXISTS wide_columns(
            table_name TEXT,
            column_name TEXT,
            vtype TINYINT,
            PRIMARY KEY (table_name, column_name))
    """)
    DBInterface.execute(connection, """
        CREATE TABLE IF NOT EXISTS source_item_metadata(item_key BIGINT PRIMARY KEY)
    """)
    DBInterface.execute(connection, """
        CREATE TABLE IF NOT EXISTS source_collection_metadata(collection_key TEXT PRIMARY KEY)
    """)
    DBInterface.execute(connection, """
        CREATE TABLE IF NOT EXISTS analyzed_item_metadata(item_key BIGINT PRIMARY KEY)
    """)
    DBInterface.execute(connection, """
        CREATE TABLE IF NOT EXISTS analyzed_collection_metadata(collection_key TEXT PRIMARY KEY)
    """)
    DBInterface.execute(connection, """
        CREATE TABLE IF NOT EXISTS item_data(
            item_key BIGINT,
            stage TINYINT,
            storage_id USMALLINT,
            seq UINTEGER,
            PRIMARY KEY (item_key, stage),
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

"""
Store one source item's interpretation index rows: its identity row plus each record's rows.

The identity row belongs to the interpretation, not to its records: a source item that interprets
into zero records still persists its fingerprint, so a reopen recognizes it as unchanged instead of
rediscovering it as new. Returns type-conflict messages for keys dropped from the metadata rows.
"""
function store_interpreted_records!(
    cache::CacheDB,
    source_item::AbstractDataSourceItem,
    records::Vector{ItemRecord},
    effective::Vector{<:AbstractDataItem},
)::Vector{String}
    started = time_ns()
    id = source_item_id(source_item)
    hex, hash = _encode_fingerprint(fingerprint(source_item))
    append!(cache.source_items, id, SourceItemRow(
        id, hex, hash, source_item_path(source_item), source_item_timestamp(source_item)))
    dropped = WideConflict[]
    for (record, _) in zip(records, effective)
        key = item_key!(cache, record.id; mint=true)
        append!(cache.items, record.id,
            ItemRow(record.id, key, record.source_item_id, record.item_label,
                String(record.kind), record.collection, _serialize_hex(record.item_fingerprint)))
        append!(dropped, edit!(cache.source_item_metadata, key, record.metadata))
    end
    record_cache_phase!(
        cache.metrics.interpreted_write_ns,
        cache.metrics.interpreted_writes,
        started,
    )
    return _conflict_messages(unique(dropped))
end

"""Store one source item's interpreted data in the memory-only interpreted buffer."""
function store_interpreted_data!(
    cache::AbstractCacheDB,
    records::Vector{ItemRecord},
    data::Vector{<:AbstractDataItem},
)::Nothing
    length(records) == length(data) || throw(ArgumentError(
        "Cannot store interpreted data: received $(length(records)) records and " *
        "$(length(data)) items; the vectors must align one-to-one",
    ))
    for (record, item) in zip(records, data)
        append!(cache.interpreted, record.id, item)
    end
    return nothing
end

"""
Store one source item's interpreted result: records to disk-backed buffers, payloads to memory.

`data` carries each record's interpreted payload; entries land in `source_item_metadata`. Returns the
interpretation-time type-conflict messages.
"""
function store_interpreted!(cache::AbstractCacheDB, source_item::AbstractDataSourceItem,
        records::Vector{ItemRecord}, data::Vector{<:AbstractDataItem})::Vector{String}
    conflicts = store_interpreted_records!(cache, source_item, records, data)
    store_interpreted_data!(cache, records, data)
    return conflicts
end

function store_interpreted_records!(
    cache::MemoryCacheDB,
    source_item::AbstractDataSourceItem,
    records::Vector{ItemRecord},
    ::Vector{<:AbstractDataItem},
)::Vector{String}
    lock(cache.lock) do
        push!(cache.source_items, source_item_id(source_item))
        for record in records
            push!(cache.items, record.id)
        end
    end
    return String[]
end

"""
Store one processed payload for a record at a payload stage (`:processed` or `:collection_processed`).

Cacheable DataFrame payloads land in the payload family; non-cacheable ones stay in the memory-only
processed buffer. Only the `:processed` stage tracks the `PROCESSING_RESULT` ledger.
"""
function store_processed!(
    cache::CacheDB,
    record::ItemRecord,
    item::AbstractDataItem;
    stage::Symbol=:processed,
)::Nothing
    key = item_key!(cache, record.id)
    payload = item_data(item)
    disk = cacheable(item) && payload isa AbstractDataFrame && !isempty(names(payload))
    if disk
        started = time_ns()
        stage === :processed && delete!(cache.processed_memory, record.id)
        append!(cache.payload, (key, _payload_stage_code(stage)), payload)
        record_cache_phase!(
            cache.metrics.processed_write_ns,
            cache.metrics.processed_writes,
            started,
        )
        stage === :processed && edit!(cache.result_states, (Int8(PROCESSING_RESULT), record.id),
            CachedResultState(
                Int8(PROCESSING_RESULT),
                record.id,
                Int8(RESULT_READY),
                record.source_item_id,
                nothing,
            ))
    else
        delete!(cache.payload, (key, _payload_stage_code(stage)))
        stage === :processed && append!(cache.processed_memory, record.id, item)
        stage === :processed && delete!(cache.result_states, (Int8(PROCESSING_RESULT), record.id))
    end
    return nothing
end

function store_processed!(
    cache::MemoryCacheDB,
    record::ItemRecord,
    item::AbstractDataItem;
    stage::Symbol=:processed,
)::Nothing
    if stage === :processed
        append!(cache.processed_memory, record.id, item)
        lock(cache.lock) do
            cache.result_states[(Int8(PROCESSING_RESULT), record.id)] = CachedResultState(
                Int8(PROCESSING_RESULT),
                record.id,
                Int8(RESULT_READY),
                record.source_item_id,
                nothing,
            )
        end
    elseif stage === :collection_processed
        append!(cache.collection_processed_memory, record.id, item)
    end
    return nothing
end

"""Persist one independent failed work result."""
function store_result_failure!(cache::CacheDB, kind::CacheResultKind,
        entity::AbstractString, source_item_id::AbstractString,
        message::AbstractString)::Nothing
    edit!(cache.result_states, (Int8(kind), String(entity)),
        CachedResultState(Int8(kind), String(entity), Int8(RESULT_FAILED),
            String(source_item_id), String(message)))
    return nothing
end

"""Persist one failed source interpretation."""
function store_source_item_failure!(
    cache::CacheDB,
    source_item_id::AbstractString,
    message::AbstractString,
)::Nothing
    source_id = String(source_item_id)
    append!(
        cache.failures,
        (source_id, source_id),
        FailureRow(source_id, source_id, String(message)),
    )
    return nothing
end

"""Format wide-column type conflicts as the messages the failure channels surface."""
_conflict_messages(dropped::Vector{WideConflict})::Vector{String} = String[
    "metadata :$name expected $(_vtype_julia(registered)), got $(_vtype_julia(incoming)); " *
    "value dropped"
    for (name, registered, incoming) in dropped
]

"""Persist one item's delivered metadata. Returns type-conflict messages for dropped keys."""
function store_item_metadata!(
    cache::CacheDB, record::ItemRecord, effective::AbstractDict,
)::Vector{String}
    started = time_ns()
    key = item_key!(cache, record.id)
    dropped = edit!(cache.analyzed_item_metadata, key, metadata_dict(effective))
    edit!(cache.result_states, (Int8(ITEM_ANALYSIS_RESULT), record.id),
        CachedResultState(Int8(ITEM_ANALYSIS_RESULT), record.id, Int8(RESULT_READY),
            record.source_item_id, nothing))
    record_cache_phase!(
        cache.metrics.metadata_write_ns,
        cache.metrics.metadata_writes,
        started,
    )
    return _conflict_messages(dropped)
end

"""Persist one collection's analyze output. Returns type-conflict messages for dropped keys."""
function store_collection_metadata!(
    cache::CacheDB, collection_key::AbstractString, analysis::AbstractDict,
)::Vector{String}
    entity = String(collection_key)
    dropped = edit!(cache.analyzed_collection_metadata, entity, metadata_dict(analysis))
    edit!(cache.result_states, (Int8(COLLECTION_ANALYSIS_RESULT), entity),
        CachedResultState(
            Int8(COLLECTION_ANALYSIS_RESULT), entity, Int8(RESULT_READY), "", nothing))
    return _conflict_messages(dropped)
end

"""Persist one collection-process result state (payloads and member metadata land separately)."""
function store_collection_process_result!(
    cache::CacheDB, collection_key::AbstractString,
)::Nothing
    entity = String(collection_key)
    edit!(cache.result_states, (Int8(COLLECTION_PROCESS_RESULT), entity),
        CachedResultState(
            Int8(COLLECTION_PROCESS_RESULT), entity, Int8(RESULT_READY), "", nothing))
    return nothing
end

function store_result_failure!(
    cache::MemoryCacheDB,
    kind::CacheResultKind,
    entity::AbstractString,
    source_item_id::AbstractString,
    message::AbstractString,
)::Nothing
    entity_value = String(entity)
    source_id = String(source_item_id)
    lock(cache.lock) do
        cache.failures[(entity_value, source_id)] = String(message)
        cache.result_states[(Int8(kind), entity_value)] = CachedResultState(
            Int8(kind), entity_value, Int8(RESULT_FAILED), source_id, String(message))
    end
    return nothing
end

function store_source_item_failure!(
    cache::MemoryCacheDB,
    source_item_id::AbstractString,
    message::AbstractString,
)::Nothing
    source_id = String(source_item_id)
    lock(cache.lock) do
        cache.failures[(source_id, source_id)] = String(message)
    end
    return nothing
end

function store_item_metadata!(cache::MemoryCacheDB, record::ItemRecord, ::AbstractDict)::Vector{String}
    lock(cache.lock) do
        cache.result_states[(Int8(ITEM_ANALYSIS_RESULT), record.id)] = CachedResultState(
            Int8(ITEM_ANALYSIS_RESULT),
            record.id,
            Int8(RESULT_READY),
            record.source_item_id,
            nothing,
        )
    end
    return String[]
end

function store_collection_metadata!(cache::MemoryCacheDB, collection_key::AbstractString, ::AbstractDict)::Vector{String}
    entity = String(collection_key)
    lock(cache.lock) do
        cache.result_states[(Int8(COLLECTION_ANALYSIS_RESULT), entity)] = CachedResultState(
            Int8(COLLECTION_ANALYSIS_RESULT), entity, Int8(RESULT_READY), "", nothing)
    end
    return String[]
end

function store_collection_process_result!(cache::MemoryCacheDB, collection_key::AbstractString)::Nothing
    entity = String(collection_key)
    lock(cache.lock) do
        cache.result_states[(Int8(COLLECTION_PROCESS_RESULT), entity)] = CachedResultState(
            Int8(COLLECTION_PROCESS_RESULT), entity, Int8(RESULT_READY), "", nothing)
    end
    return nothing
end

"""Write one source collection-metadata row."""
function edit_source_collection_metadata!(cache::CacheDB, key::AbstractString, metadata::AbstractDict)::Nothing
    edit!(cache.source_collection_metadata, String(key), metadata_dict(metadata))
    return nothing
end

"""Write one source item-metadata (entries) row."""
function edit_source_item_metadata!(cache::CacheDB, key::Int64, metadata::AbstractDict)::Nothing
    edit!(cache.source_item_metadata, key, metadata_dict(metadata))
    return nothing
end

edit_source_collection_metadata!(::MemoryCacheDB, ::AbstractString, ::AbstractDict)::Nothing = nothing
edit_source_item_metadata!(::MemoryCacheDB, ::Int64, ::AbstractDict)::Nothing = nothing

"""Delete every cached result owned by one published source item, keyed by item surrogate."""
function delete_source_item!(
    cache::CacheDB,
    source_item_id::AbstractString,
    old_records::Vector{ItemRecord},
)::Nothing
    source_id = String(source_item_id)
    delete!(cache.source_items, source_id)
    for record in old_records
        key = item_key!(cache, record.id)
        delete!(cache.items, record.id)
        delete!(cache.source_item_metadata, key)
        delete!(cache.analyzed_item_metadata, key)
        delete!(cache.payload, (key, PAYLOAD_STAGE_PROCESSED))
        delete!(cache.payload, (key, PAYLOAD_STAGE_COLLECTION_PROCESSED))
        delete!(cache.interpreted, record.id)
        delete!(cache.processed_memory, record.id)
        delete!(cache.failures, (record.id, source_id))
        delete!(cache.result_states, (Int8(PROCESSING_RESULT), record.id))
        delete!(cache.result_states, (Int8(ITEM_ANALYSIS_RESULT), record.id))
    end
    delete!(cache.failures, (source_id, source_id))
    return nothing
end

function delete_source_item!(
    cache::MemoryCacheDB,
    source_item_id::AbstractString,
    old_records::Vector{ItemRecord},
)::Nothing
    source_id = String(source_item_id)
    lock(cache.lock) do
        delete!(cache.source_items, source_id)
        delete!(cache.failures, (source_id, source_id))
    end
    for record in old_records
        delete!(cache.interpreted, record.id)
        delete!(cache.processed_memory, record.id)
        delete!(cache.collection_processed_memory, record.id)
        lock(cache.lock) do
            delete!(cache.items, record.id)
            delete!(cache.failures, (record.id, source_id))
            delete!(cache.result_states, (Int8(PROCESSING_RESULT), record.id))
            delete!(cache.result_states, (Int8(ITEM_ANALYSIS_RESULT), record.id))
        end
    end
    return nothing
end

"""Delete cached analysis for collections whose published metadata is invalid."""
function delete_collection_metadata!(
    cache::CacheDB,
    collection_keys::Vector{String},
)::Nothing
    for collection_key in unique(collection_keys)
        delete!(cache.analyzed_collection_metadata, collection_key)
        delete!(cache.result_states, (Int8(COLLECTION_ANALYSIS_RESULT), collection_key))
        delete!(cache.result_states, (Int8(COLLECTION_PROCESS_RESULT), collection_key))
    end
    return nothing
end

function delete_collection_metadata!(
    cache::MemoryCacheDB,
    collection_keys::Vector{String},
)::Nothing
    for collection_key in unique(collection_keys)
        key = String(collection_key)
        lock(cache.lock) do
            delete!(cache.result_states, (Int8(COLLECTION_ANALYSIS_RESULT), key))
            delete!(cache.result_states, (Int8(COLLECTION_PROCESS_RESULT), key))
        end
    end
    return nothing
end

"""Aggregate staged durable keys and payload rows."""
function cache_pending_counts(cache::CacheDB)::NamedTuple{(:items, :rows),Tuple{Int,Int}}
    items = 0
    rows = 0
    for buffer in
        (cache.source_items, cache.items, cache.source_item_metadata,
         cache.source_collection_metadata, cache.analyzed_item_metadata,
         cache.analyzed_collection_metadata, cache.failures, cache.result_states, cache.payload)
        lock(buffer.flush_condition) do
            items += length(buffer.queued) + length(buffer.writing)
            rows += buffer.queued_rows + buffer.writing_rows
        end
    end
    return (items=items, rows=rows)
end

cache_pending_counts(::MemoryCacheDB)::NamedTuple{(:items, :rows),Tuple{Int,Int}} =
    (items=0, rows=0)

"""Return persisted stage counts owned by the cache layer."""
function cache_stage_summary(cache::CacheDB)::CacheStageSummary
    source_items = read(cache.source_items)
    items = read(cache.items)
    failures = read(cache.failures)
    states = read(cache.result_states)
    return _cache_stage_summary(source_items, items, failures, states)
end

function cache_stage_summary(cache::MemoryCacheDB)::CacheStageSummary
    return lock(cache.lock) do
        _cache_stage_summary(
            copy(cache.source_items),
            copy(cache.items),
            copy(cache.failures),
            copy(cache.result_states),
        )
    end
end

function _cache_stage_summary(source_items, items, failures, states)::CacheStageSummary
    ready(kind::CacheResultKind)::Int = count(
        state -> CacheResultKind(state.kind) === kind &&
            CacheResultStatus(state.status) === RESULT_READY,
        values(states),
    )
    failed(kind::CacheResultKind)::Int = count(
        state -> CacheResultKind(state.kind) === kind &&
            CacheResultStatus(state.status) === RESULT_FAILED,
        values(states),
    )
    failed_interpret = count(key -> first(key) == last(key), keys(failures))
    return CacheStageSummary(
        length(source_items),
        length(items),
        ready(PROCESSING_RESULT),
        ready(ITEM_ANALYSIS_RESULT),
        ready(COLLECTION_PROCESS_RESULT),
        ready(COLLECTION_ANALYSIS_RESULT),
        failed_interpret,
        failed(PROCESSING_RESULT),
        failed(ITEM_ANALYSIS_RESULT),
        failed(COLLECTION_PROCESS_RESULT) + failed(COLLECTION_ANALYSIS_RESULT),
    )
end

"""Whether any cache store has queued or writing mutations."""
function cache_has_pending_writes(cache::AbstractCacheDB)::Bool
    pending = cache_pending_counts(cache)
    return pending.items > 0 || pending.rows > 0
end

"""
Read item data through the buffer for a stage (`:interpreted`, `:processed`, or
`:collection_processed`). Returns a vector aligned to `records` of loaded items or `nothing` on a
miss; the caller falls back to the source. Interpreted payloads are memory-only, so a memory miss is
final. `:collection_processed` is disk-backed on `CacheDB` and memory-backed on `MemoryCacheDB`.
"""
function read_item_data(cache::CacheDB, records::Vector{ItemRecord}; stage::Symbol)::Vector{Any}
    stage in (:interpreted, :processed, :collection_processed) ||
        throw(ArgumentError("unknown item-data cache stage '$stage'"))
    isempty(records) && return Any[]

    if stage === :interpreted
        return Any[read(cache.interpreted, record.id) for record in records]
    end

    # Raw payload reads: at runtime delivery consults the live graph, then reads committed payloads.
    # The `result_states` ledger is authoritative for finished work between sessions.
    stage_code = _payload_stage_code(stage)
    loaded = Vector{Any}(undef, length(records))
    disk_keys = Tuple{Int64,Int8}[]
    for (index, record) in pairs(records)
        if stage === :processed
            held = read(cache.processed_memory, record.id)
            if held !== nothing
                loaded[index] = held
                continue
            end
        end
        push!(disk_keys, (item_key!(cache, record.id), stage_code))
        loaded[index] = nothing
    end

    disk_data = isempty(disk_keys) ?
        Dict{Tuple{Int64,Int8},Union{Nothing,AbstractDataFrame}}() :
        read(cache.payload, unique(disk_keys))
    for (index, record) in pairs(records)
        loaded[index] !== nothing && continue
        data = get(disk_data, (item_key!(cache, record.id), stage_code), nothing)
        loaded[index] = data === nothing ? nothing : DataItem(record, data)
    end
    return loaded
end

function read_item_data(cache::MemoryCacheDB, records::Vector{ItemRecord}; stage::Symbol)::Vector{Any}
    stage in (:interpreted, :processed, :collection_processed) ||
        throw(ArgumentError("unknown item-data cache stage '$stage'"))
    store = if stage === :interpreted
        cache.interpreted
    elseif stage === :collection_processed
        cache.collection_processed_memory
    else
        cache.processed_memory
    end
    return Any[read(store, record.id) for record in records]
end

"""Close a workspace cache's database file, checkpointing first."""
function close_cache_db!(cachedb::CacheDB)::Nothing
    failure::Union{Nothing,Exception} = nothing
    for buffer in
        (cachedb.source_items, cachedb.items, cachedb.source_item_metadata,
         cachedb.source_collection_metadata, cachedb.analyzed_item_metadata,
         cachedb.analyzed_collection_metadata, cachedb.failures,
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

function close_cache_db!(cachedb::MemoryCacheDB)::Nothing
    close!(cachedb.interpreted)
    close!(cachedb.processed_memory)
    close!(cachedb.collection_processed_memory)
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

"""Hex-text storage and content hash for one fingerprint, as `(hex, hash)` (`nothing` when absent)."""
function _encode_fingerprint(fingerprint)::Tuple{Union{Nothing,String},Union{Nothing,String}}
    fingerprint === nothing && return (nothing, nothing)
    bytes = _serialize_bytes(fingerprint)
    return (bytes2hex(bytes), bytes2hex(sha1(bytes)))
end

# Index construction + comparison (in memory)

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
        Dict{String,MetadataDict}(),
        Dict{CacheResultKey,CachedResultState}(),
    )
end

"""Delete every index and item-data row and reset the item-key map, leaving a schema-valid cache."""
function clear_cache_index!(cachedb::CacheDB)::Nothing
    for buffer in
        (cachedb.source_items, cachedb.items, cachedb.source_item_metadata,
         cachedb.source_collection_metadata, cachedb.analyzed_item_metadata,
         cachedb.analyzed_collection_metadata, cachedb.failures,
         cachedb.result_states, cachedb.payload, cachedb.interpreted, cachedb.processed_memory)
        clear!(buffer)
    end
    lock(cachedb.key_lock) do
        empty!(cachedb.item_keys)
        cachedb.next_item_key = Int64(1)
    end
    connection = DBInterface.connect(cachedb.db)
    try
        DBInterface.execute(connection, "DELETE FROM meta")
    finally
        DBInterface.close!(connection)
    end
    return nothing
end

# Incremental writing (per source item, as the scan streams)

function clear_cache_index!(cachedb::MemoryCacheDB)::Nothing
    clear!(cachedb.interpreted)
    clear!(cachedb.processed_memory)
    clear!(cachedb.collection_processed_memory)
    lock(cachedb.lock) do
        empty!(cachedb.source_items)
        empty!(cachedb.items)
        empty!(cachedb.failures)
        empty!(cachedb.result_states)
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

write_meta_header!(::MemoryCacheDB)::Nothing = nothing

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
        throw(ProjectCacheSchemaError(
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

cache_built(::MemoryCacheDB)::Bool = false

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

"""
Load the previously persisted source-item fingerprints from the `source_items` table, keyed by
source-item id. A memory-only cache holds nothing across sessions, so every discovered item
re-interprets.
"""
function _load_source_item_fingerprints(cache::CacheDB)::Dict{String,Any}
    fingerprints = Dict{String,Any}()
    for row in values(read(cache.source_items))
        fingerprints[row.id] = _deserialize_hex(row.fingerprint_hex)
    end
    return fingerprints
end

_load_source_item_fingerprints(::MemoryCacheDB)::Dict{String,Any} = Dict{String,Any}()

"""Reconstruct the full `SourceScan` stored in one DuckDB cache (entries-layer metadata only)."""
function _load_source_scan(
    cache::CacheDB,
    item_metadata_by_key::Dict{Int64,MetadataDict},
    collection_analysis::Dict{String,MetadataDict},
)::SourceScan
    identity = cache.identity
    _load_meta(cache)
    source_rows = read(cache.source_items)
    item_rows = read(cache.items)
    failure_rows = read(cache.failures)

    all_source_ids = String[String(id) for id in keys(source_rows)]
    locations = Dict{String,Tuple{Union{Nothing,String},Union{Nothing,DateTime}}}()
    failures = ItemFailure[]
    for row in values(source_rows)
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
            source_item_path=path,
            source_item_timestamp=timestamp,
            item_label=row.item_label,
            kind=Symbol(row.kind),
            collection=row.collection,
            metadata=get(item_metadata_by_key, row.item_key, MetadataDict()),
            item_fingerprint=_deserialize_hex(row.item_fingerprint_hex),
        ))
    end

    sort!(records; by=record -> (record.source_item_id, record.id))
    # Scan-wide summaries are derived, not stored: a source item with no items was skipped or failed.
    source_ids_with_items = Set(record.source_item_id for record in records)
    skipped_count = count(id -> !(id in source_ids_with_items), all_source_ids)
    hierarchy = Hierarchy(identity.source_id, false, skipped_count)
    for record in records
        insert_item!(hierarchy, record)
    end
    for (entity, dict) in collection_analysis
        node = get(hierarchy.index, collection_path_tuple(entity), nothing)
        node === nothing && continue
        merge!(node.analysis, dict)
    end
    sort!(hierarchy)

    return SourceScan(
        identity.source_id,
        identity.source_label,
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
        entries_by_key = read(cachedb.source_item_metadata)
        analyzed_by_key = read(cachedb.analyzed_item_metadata)
        collection_analysis = read(cachedb.analyzed_collection_metadata)
        key_to_id = Dict{Int64,String}(
            row.item_key => row.id for row in values(read(cachedb.items)))
        # The computed layer is the delivered effective minus the entries layer: what analyze and
        # collection-process added, restored per item id.
        item_metadata = Dict{String,MetadataDict}()
        for (key, analyzed) in analyzed_by_key
            id = get(key_to_id, key, nothing)
            id === nothing && continue
            entries = get(entries_by_key, key, MetadataDict())
            computed = MetadataDict()
            for (name, value) in analyzed
                get(entries, name, nothing) == value || (computed[name] = value)
            end
            isempty(computed) || (item_metadata[id] = computed)
        end
        base = ProjectCacheIndex(
            identity, _load_source_scan(cachedb, entries_by_key, collection_analysis))
        result_states = Dict{CacheResultKey,CachedResultState}()
        for state in values(read(cachedb.result_states))
            result_states[CacheResultKey(CacheResultKind(state.kind), state.entity)] = state
        end
        ProjectCacheIndex(
            base.identity,
            base.source,
            base.analysis_errors,
            item_metadata,
            result_states,
        )
    end
    emit_progress(
        on_progress;
        phase=:cache_load,
        total_source_items=1,
        processed_source_items=1,
        loaded_items=length(all_items(index.source.hierarchy)),
        skipped_source_items=index.source.hierarchy.skipped_count,
        current_source_item=index.identity.cache_path,
    )
    return index
end

# Query surface

"""
Return the ids of items whose delivered effective metadata satisfies a SQL predicate.

Queries committed DB state over `items LEFT JOIN analyzed_item_metadata USING (item_key)`; the
buffer's ≤2s staleness applies. Columns are exactly the metadata names users coded.
"""
function query_items(cache::CacheDB, predicate::AbstractString)::Vector{String}
    _ensure_items_view!(cache)
    connection = DBInterface.connect(cache.db)
    try
        return String[String(row.id) for row in DBInterface.execute(
            connection, "SELECT id FROM mb_item_query WHERE $(predicate)")]
    finally
        DBInterface.close!(connection)
    end
end

query_items(::MemoryCacheDB, ::AbstractString)::Vector{String} = String[]

"""Recreate the query view when the effective-metadata column signature changed."""
function _ensure_items_view!(cache::CacheDB)::Nothing
    connection = DBInterface.connect(cache.db)
    try
        signature = Symbol[]
        for row in DBInterface.execute(connection, """
            SELECT column_name FROM information_schema.columns
            WHERE table_schema = 'main' AND table_name = 'analyzed_item_metadata'
            ORDER BY ordinal_position
        """)
            push!(signature, Symbol(String(row.column_name)))
        end
        signature == cache.query_view_signature && return nothing
        DBInterface.execute(connection, "DROP VIEW IF EXISTS mb_item_query")
        DBInterface.execute(connection, """
            CREATE VIEW mb_item_query AS
            SELECT * FROM items LEFT JOIN analyzed_item_metadata USING (item_key)
        """)
        cache.query_view_signature = signature
    finally
        DBInterface.close!(connection)
    end
    return nothing
end
