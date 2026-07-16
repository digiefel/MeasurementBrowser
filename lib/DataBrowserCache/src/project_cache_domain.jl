const PROJECT_CACHE_SCHEMA_VERSION = 18

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

"""Error raised when generated collection values can no longer be reconstructed safely."""
struct ProjectCacheDataError <: Exception
    path::String
    message::String
end

Base.showerror(io::IO, err::ProjectCacheDataError)::Nothing =
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
    entity::Union{String,Int64}
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

"""One persisted collection result-state row, keyed by the private collection-record key."""
struct CachedCollectionResultState
    kind::Int8
    entity::Int64
    status::Int8
    source_item_id::String
    message::Union{Nothing,String}
end

CachedCollectionResultState(row)::CachedCollectionResultState = CachedCollectionResultState(
    Int8(row.kind),
    Int64(row.entity),
    Int8(row.status),
    String(row.source_item_id),
    _null_to_nothing(row.message),
)

const AnyCachedResultState = Union{CachedResultState,CachedCollectionResultState}

"""All data-less content needed to restore and compare a project cache."""
struct ProjectCacheIndex
    identity::ProjectCacheIdentity
    source::SourceScan
    analysis_errors::Dict{String,String}
    # The computed layers (analyze output + collection-process overwrite) per item, restored for
    # items whose governing result is valid: delivered effective minus the entries layer.
    item_metadata::Dict{String,MetadataDict}
    result_states::Dict{CacheResultKey,AnyCachedResultState}
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

mutable struct CacheStageLedger
    summary::CacheStageSummary
    source_items::Set{String}
    items::Set{String}
    failures::Set{Tuple{String,String}}
    result_states::Dict{CacheResultKey,AnyCachedResultState}
    lock::ReentrantLock
end

CacheStageLedger()::CacheStageLedger = CacheStageLedger(
    CacheStageSummary(),
    Set{String}(),
    Set{String}(),
    Set{Tuple{String,String}}(),
    Dict{CacheResultKey,AnyCachedResultState}(),
    ReentrantLock(),
)

function CacheStageLedger(source_items, items, failures, states...)::CacheStageLedger
    result_states = Dict{CacheResultKey,AnyCachedResultState}()
    for store in states, state in values(store)
        result_states[CacheResultKey(CacheResultKind(state.kind), state.entity)] = state
    end
    return CacheStageLedger(
        _cache_stage_summary(source_items, items, failures, states...),
        Set{String}(String(id) for id in keys(source_items)),
        Set{String}(String(id) for id in keys(items)),
        Set{Tuple{String,String}}((String(first(key)), String(last(key))) for key in keys(failures)),
        result_states,
        ReentrantLock(),
    )
end

function _stage_summary_delta(
    summary::CacheStageSummary;
    cached_sources::Int=0,
    interpreted_items::Int=0,
    processed::Int=0,
    analyzed::Int=0,
    collection_processed::Int=0,
    collection_analyzed::Int=0,
    failed_interpret::Int=0,
    failed_process::Int=0,
    failed_analyze::Int=0,
    failed_collection::Int=0,
)::CacheStageSummary
    return CacheStageSummary(
        summary.cached_sources + cached_sources,
        summary.interpreted_items + interpreted_items,
        summary.processed + processed,
        summary.analyzed + analyzed,
        summary.collection_processed + collection_processed,
        summary.collection_analyzed + collection_analyzed,
        summary.failed_interpret + failed_interpret,
        summary.failed_process + failed_process,
        summary.failed_analyze + failed_analyze,
        summary.failed_collection + failed_collection,
    )
end

function _stage_result_delta(
    summary::CacheStageSummary,
    state::AnyCachedResultState,
    sign::Int,
)::CacheStageSummary
    kind = CacheResultKind(state.kind)
    status = CacheResultStatus(state.status)
    if status === RESULT_READY
        kind === PROCESSING_RESULT &&
            return _stage_summary_delta(summary; processed=sign)
        kind === ITEM_ANALYSIS_RESULT &&
            return _stage_summary_delta(summary; analyzed=sign)
        kind === COLLECTION_PROCESS_RESULT &&
            return _stage_summary_delta(summary; collection_processed=sign)
        kind === COLLECTION_ANALYSIS_RESULT &&
            return _stage_summary_delta(summary; collection_analyzed=sign)
    elseif status === RESULT_FAILED
        kind === PROCESSING_RESULT &&
            return _stage_summary_delta(summary; failed_process=sign)
        kind === ITEM_ANALYSIS_RESULT &&
            return _stage_summary_delta(summary; failed_analyze=sign)
        kind in (COLLECTION_PROCESS_RESULT, COLLECTION_ANALYSIS_RESULT) &&
            return _stage_summary_delta(summary; failed_collection=sign)
    end
    return summary
end

function _stage_ledger_source!(ledger::CacheStageLedger, id::AbstractString, present::Bool)::Nothing
    source_id = String(id)
    lock(ledger.lock) do
        exists = source_id in ledger.source_items
        if present && !exists
            push!(ledger.source_items, source_id)
            ledger.summary = _stage_summary_delta(ledger.summary; cached_sources=1)
        elseif !present && exists
            delete!(ledger.source_items, source_id)
            ledger.summary = _stage_summary_delta(ledger.summary; cached_sources=-1)
        end
    end
    return nothing
end

function _stage_ledger_item!(ledger::CacheStageLedger, id::AbstractString, present::Bool)::Nothing
    item_id = String(id)
    lock(ledger.lock) do
        exists = item_id in ledger.items
        if present && !exists
            push!(ledger.items, item_id)
            ledger.summary = _stage_summary_delta(ledger.summary; interpreted_items=1)
        elseif !present && exists
            delete!(ledger.items, item_id)
            ledger.summary = _stage_summary_delta(ledger.summary; interpreted_items=-1)
        end
    end
    return nothing
end

function _stage_ledger_failure!(
    ledger::CacheStageLedger,
    key::Tuple{String,String},
    present::Bool,
)::Nothing
    lock(ledger.lock) do
        exists = key in ledger.failures
        delta = first(key) == last(key) ? 1 : 0
        if present && !exists
            push!(ledger.failures, key)
            delta == 1 &&
                (ledger.summary = _stage_summary_delta(ledger.summary; failed_interpret=1))
        elseif !present && exists
            delete!(ledger.failures, key)
            delta == 1 &&
                (ledger.summary = _stage_summary_delta(ledger.summary; failed_interpret=-1))
        end
    end
    return nothing
end

function _stage_ledger_result!(
    ledger::CacheStageLedger,
    key::CacheResultKey,
    state::Union{Nothing,AnyCachedResultState},
)::Nothing
    lock(ledger.lock) do
        previous = get(ledger.result_states, key, nothing)
        previous === nothing ||
            (ledger.summary = _stage_result_delta(ledger.summary, previous, -1))
        if state === nothing
            delete!(ledger.result_states, key)
        else
            ledger.result_states[key] = state
            ledger.summary = _stage_result_delta(ledger.summary, state, 1)
        end
    end
    return nothing
end

function _stage_ledger_result!(
    ledger::CacheStageLedger,
    key::Tuple{Int8,T},
    state::Union{Nothing,AnyCachedResultState},
)::Nothing where {T<:Union{AbstractString,Integer}}
    entity = key[2] isa AbstractString ? String(key[2]) : Int64(key[2])
    return _stage_ledger_result!(
        ledger, CacheResultKey(CacheResultKind(key[1]), entity), state)
end

function _reset_stage_ledger!(ledger::CacheStageLedger)::Nothing
    lock(ledger.lock) do
        empty!(ledger.source_items)
        empty!(ledger.items)
        empty!(ledger.failures)
        empty!(ledger.result_states)
        ledger.summary = CacheStageSummary()
    end
    return nothing
end

# Storage rows (field names are the table columns; `R(database_row)` is the pure decoder)

struct SourceItemRow
    id::String
    fingerprint_hex::Union{Nothing,String}
    path::Union{Nothing,String}
    timestamp::Union{Nothing,DateTime}
end

SourceItemRow(row)::SourceItemRow = SourceItemRow(
    String(row.id),
    _null_to_nothing(row.fingerprint_hex),
    _null_to_nothing(row.path),
    _null_to_nothing(row.timestamp),
)

"""One persisted package-owned collection record; no live user collection value is serialized."""
struct CollectionRow
    collection_key::Int64
    parent_key::Union{Nothing,Int64}
    id::String
    label::String
    metadata_hex::String
    registration_name::Union{Nothing,String}
end

CollectionRow(row)::CollectionRow = CollectionRow(
    Int64(row.collection_key),
    _null_to_nothing(row.parent_key),
    String(row.id),
    String(row.label),
    String(row.metadata_hex),
    _null_to_nothing(row.registration_name),
)

struct ItemRow
    id::String
    item_key::Int64
    source_item_id::String
    label::String
    kind::String
    collection_key::Union{Nothing,Int64}
end

ItemRow(row)::ItemRow = ItemRow(
    String(row.id),
    Int64(row.item_key),
    String(row.source_item_id),
    String(row.label),
    String(row.kind),
    _null_to_nothing(row.collection_key),
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
        joinpath(depot, "databrowser", name, "cache.duckdb"),
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
    collections::RowStore{Int64,CollectionRow}
    items::RowStore{String,ItemRow}
    source_item_metadata::WideRowStore{Int64}
    analyzed_item_metadata::WideRowStore{Int64}
    analyzed_collection_metadata::WideRowStore{Int64}
    failures::RowStore{Tuple{String,String},FailureRow}
    result_states::RowStore{Tuple{Int8,String},CachedResultState}
    collection_result_states::RowStore{Tuple{Int8,Int64},CachedCollectionResultState}
    payload::TabularFamilyStore
    interpreted::MemoryStore{String,AbstractDataItem}
    processed_memory::MemoryStore{String,AbstractDataItem}
    collection_processed_memory::MemoryStore{String,AbstractDataItem}
    # One integer surrogate per item id, shared with the payload store; minted at interpretation.
    item_keys::Dict{String,Int64}
    next_item_key::Int64
    persisted_collection_keys::Set{Int64}
    key_lock::ReentrantLock
    stage_ledger::CacheStageLedger
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
    collection_result_states::Dict{Tuple{Int8,Int64},CachedCollectionResultState}
    stage_ledger::CacheStageLedger
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
        source_items = RowStore{String,SourceItemRow}(db, "source_items", ("id",); profiler)
        collections = RowStore{Int64,CollectionRow}(
            db, "collections", ("collection_key",); profiler)
        items = RowStore{String,ItemRow}(db, "items", ("id",); profiler)
        failures = RowStore{Tuple{String,String},FailureRow}(
            db, "item_failures", ("item_id", "source_item_id"); profiler)
        result_states = RowStore{Tuple{Int8,String},CachedResultState}(
            db, "result_states", ("kind", "entity"); profiler)
        collection_result_states = RowStore{Tuple{Int8,Int64},CachedCollectionResultState}(
            db, "collection_result_states", ("kind", "entity"); profiler)
        payload = TabularFamilyStore(db; profiler, row_limit=CACHE_BUFFER_ROW_LIMIT)
        item_keys, next_item_key = _load_item_keys(db)
        stage_ledger = CacheStageLedger(
            read(source_items),
            read(items),
            read(failures),
            read(result_states),
            read(collection_result_states),
        )
        return CacheDB(
            identity,
            db,
            metrics,
            source_items,
            collections,
            items,
            WideRowStore{Int64}(db, "source_item_metadata", "item_key"; profiler),
            WideRowStore{Int64}(db, "analyzed_item_metadata", "item_key"; profiler),
            WideRowStore{Int64}(db, "analyzed_collection_metadata", "collection_key"; profiler),
            failures,
            result_states,
            collection_result_states,
            payload,
            MemoryStore{String,AbstractDataItem}(; row_limit=CACHE_BUFFER_ROW_LIMIT),
            MemoryStore{String,AbstractDataItem}(; row_limit=CACHE_BUFFER_ROW_LIMIT),
            MemoryStore{String,AbstractDataItem}(; row_limit=CACHE_BUFFER_ROW_LIMIT),
            item_keys,
            next_item_key,
            Set{Int64}(keys(read(collections))),
            ReentrantLock(),
            stage_ledger,
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
        Dict{Tuple{Int8,Int64},CachedCollectionResultState}(),
        CacheStageLedger(),
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
    return Tables.istable(data) ? _payload_rows(data) : 1
end

"""Create the cache tables if they do not already exist."""
function ensure_schema!(connection)::Nothing
    DBInterface.execute(connection, """
        CREATE TABLE IF NOT EXISTS meta(
            key TEXT PRIMARY KEY,
            value TEXT)
    """)
    DBInterface.execute(connection, create_table_sql(SourceItemRow, "source_items", (:id,)))
    DBInterface.execute(connection,
        create_table_sql(CollectionRow, "collections", (:collection_key,)))
    DBInterface.execute(connection, create_table_sql(ItemRow, "items", (:id,)))
    DBInterface.execute(connection, create_table_sql(FailureRow, "item_failures", (:item_id, :source_item_id)))
    DBInterface.execute(connection, create_table_sql(CachedResultState, "result_states", (:kind, :entity)))
    DBInterface.execute(connection, create_table_sql(
        CachedCollectionResultState, "collection_result_states", (:kind, :entity)))
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
        CREATE TABLE IF NOT EXISTS analyzed_item_metadata(item_key BIGINT PRIMARY KEY)
    """)
    DBInterface.execute(connection, """
        CREATE TABLE IF NOT EXISTS analyzed_collection_metadata(collection_key BIGINT PRIMARY KEY)
    """)
    DBInterface.execute(connection, """
        CREATE TABLE IF NOT EXISTS item_data(
            item_key BIGINT,
            stage TINYINT,
            storage_id USMALLINT,
            seq UINTEGER,
            container TEXT,
            PRIMARY KEY (item_key, stage),
            UNIQUE (storage_id, seq))
    """)
    DBInterface.execute(connection, """
        CREATE TABLE IF NOT EXISTS payload_schemas(
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
    source_item_id_value = id(source_item)
    hex = _serialize_hex(fingerprint(source_item))
    append!(cache.source_items, source_item_id_value, SourceItemRow(
        source_item_id_value, hex, source_item_path(source_item),
        source_item_timestamp(source_item)))
    _stage_ledger_source!(cache.stage_ledger, source_item_id_value, true)
    dropped = WideConflict[]
    for (record, _) in zip(records, effective)
        key = item_key!(cache, record.id; mint=true)
        append!(cache.items, record.id,
            ItemRow(record.id, key, record.source_item_id, record.label,
                String(record.kind), nothing))
        _stage_ledger_item!(cache.stage_ledger, record.id, true)
        append!(dropped, edit!(cache.source_item_metadata, key, record.metadata))
    end
    record_cache_phase!(
        cache.metrics.interpreted_write_ns,
        cache.metrics.interpreted_writes,
        started,
    )
    return _conflict_messages(unique(dropped))
end

"""
Persist the package-owned collection records used by newly published items and attach each item
row to its leaf. Arbitrary user-defined `AbstractCollection` values never enter the cache.
"""
function store_collection_index!(
    cache::CacheDB,
    collections::CollectionIndex,
    records::Vector{ItemRecord},
)::Nothing
    for record in records
        for collection_key in collection_path_keys(collections, record.collection_key)
            collection_record = collections.records[collection_key]
            lock(cache.key_lock) do
                row = CollectionRow(
                    collection_record.key,
                    collection_record.parent_key,
                    collection_record.id,
                    collection_record.label,
                    _serialize_hex(collection_record.own_metadata),
                    collection_record.registration_name,
                )
                if !(collection_record.key in cache.persisted_collection_keys)
                    append!(cache.collections, collection_record.key, row)
                    push!(cache.persisted_collection_keys, collection_record.key)
                else
                    edit!(cache.collections, collection_record.key, row)
                end
            end
        end
        key = item_key!(cache, record.id)
        edit!(cache.items, record.id,
            ItemRow(record.id, key, record.source_item_id, record.label,
                String(record.kind), record.collection_key))
    end
    return nothing
end

store_collection_index!(::MemoryCacheDB, ::CollectionIndex, ::Vector{ItemRecord})::Nothing = nothing

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
        push!(cache.source_items, id(source_item))
        for record in records
            push!(cache.items, record.id)
        end
    end
    _stage_ledger_source!(cache.stage_ledger, id(source_item), true)
    for record in records
        _stage_ledger_item!(cache.stage_ledger, record.id, true)
    end
    return String[]
end

"""
Store one processed payload for a record at a payload stage (`:processed` or `:collection_processed`).

Cacheable tabular payloads with DuckDB-storable columns land in the payload family; everything
else stays in the memory-only processed buffer. Only the `:processed` stage tracks the
`PROCESSING_RESULT` ledger.
"""
function store_processed!(
    cache::CacheDB,
    record::ItemRecord,
    item::AbstractDataItem;
    stage::Symbol=:processed,
)::Nothing
    key = item_key!(cache, record.id)
    payload = item_data(item)
    disk = item isa RegisteredDataItem && cacheable(item) && _storable_table(payload)
    if disk
        started = time_ns()
        if stage === :processed
            delete!(cache.processed_memory, record.id)
        else
            delete!(cache.collection_processed_memory, record.id)
        end
        append!(cache.payload, (key, _payload_stage_code(stage)), payload)
        record_cache_phase!(
            cache.metrics.processed_write_ns,
            cache.metrics.processed_writes,
            started,
        )
        if stage === :processed
            state = CachedResultState(
                Int8(PROCESSING_RESULT),
                record.id,
                Int8(RESULT_READY),
                record.source_item_id,
                nothing,
            )
            edit!(cache.result_states, (Int8(PROCESSING_RESULT), record.id), state)
            _stage_ledger_result!(
                cache.stage_ledger, (Int8(PROCESSING_RESULT), record.id), state)
        end
    else
        delete!(cache.payload, (key, _payload_stage_code(stage)))
        if stage === :processed
            append!(cache.processed_memory, record.id, item)
        else
            append!(cache.collection_processed_memory, record.id, item)
        end
        if stage === :processed
            delete!(cache.result_states, (Int8(PROCESSING_RESULT), record.id))
            _stage_ledger_result!(
                cache.stage_ledger, (Int8(PROCESSING_RESULT), record.id), nothing)
        end
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
            state = CachedResultState(
                Int8(PROCESSING_RESULT),
                record.id,
                Int8(RESULT_READY),
                record.source_item_id,
                nothing,
            )
            cache.result_states[(Int8(PROCESSING_RESULT), record.id)] = state
            _stage_ledger_result!(
                cache.stage_ledger, (Int8(PROCESSING_RESULT), record.id), state)
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
    entity_value = String(entity)
    state = CachedResultState(Int8(kind), entity_value, Int8(RESULT_FAILED),
        String(source_item_id), String(message))
    edit!(cache.result_states, (Int8(kind), entity_value), state)
    _stage_ledger_result!(cache.stage_ledger, (Int8(kind), entity_value), state)
    return nothing
end

function store_result_failure!(cache::CacheDB, kind::CacheResultKind,
        entity::Integer, source_item_id::AbstractString,
        message::AbstractString)::Nothing
    entity_value = Int64(entity)
    state = CachedCollectionResultState(Int8(kind), entity_value, Int8(RESULT_FAILED),
        String(source_item_id), String(message))
    edit!(cache.collection_result_states, (Int8(kind), entity_value), state)
    _stage_ledger_result!(cache.stage_ledger, (Int8(kind), entity_value), state)
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
    _stage_ledger_failure!(cache.stage_ledger, (source_id, source_id), true)
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
    state = CachedResultState(Int8(ITEM_ANALYSIS_RESULT), record.id, Int8(RESULT_READY),
        record.source_item_id, nothing)
    edit!(cache.result_states, (Int8(ITEM_ANALYSIS_RESULT), record.id), state)
    _stage_ledger_result!(
        cache.stage_ledger, (Int8(ITEM_ANALYSIS_RESULT), record.id), state)
    record_cache_phase!(
        cache.metrics.metadata_write_ns,
        cache.metrics.metadata_writes,
        started,
    )
    return _conflict_messages(dropped)
end

"""Persist one collection's analyze output. Returns type-conflict messages for dropped keys."""
function store_collection_metadata!(
    cache::CacheDB, collection_key::Integer, analysis::AbstractDict,
)::Vector{String}
    entity = Int64(collection_key)
    dropped = edit!(cache.analyzed_collection_metadata, entity, metadata_dict(analysis))
    state = CachedCollectionResultState(
        Int8(COLLECTION_ANALYSIS_RESULT), entity, Int8(RESULT_READY), "", nothing)
    edit!(cache.collection_result_states, (Int8(COLLECTION_ANALYSIS_RESULT), entity), state)
    _stage_ledger_result!(
        cache.stage_ledger, (Int8(COLLECTION_ANALYSIS_RESULT), entity), state)
    return _conflict_messages(dropped)
end

"""Persist one collection-process result state (payloads and member metadata land separately)."""
function store_collection_process_result!(
    cache::CacheDB, collection_key::Integer,
)::Nothing
    entity = Int64(collection_key)
    state = CachedCollectionResultState(
        Int8(COLLECTION_PROCESS_RESULT), entity, Int8(RESULT_READY), "", nothing)
    edit!(cache.collection_result_states, (Int8(COLLECTION_PROCESS_RESULT), entity), state)
    _stage_ledger_result!(
        cache.stage_ledger, (Int8(COLLECTION_PROCESS_RESULT), entity), state)
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
        state = CachedResultState(
            Int8(kind), entity_value, Int8(RESULT_FAILED), source_id, String(message))
        cache.result_states[(Int8(kind), entity_value)] = state
        _stage_ledger_failure!(cache.stage_ledger, (entity_value, source_id), true)
        _stage_ledger_result!(cache.stage_ledger, (Int8(kind), entity_value), state)
    end
    return nothing
end

function store_result_failure!(
    cache::MemoryCacheDB,
    kind::CacheResultKind,
    entity::Integer,
    source_item_id::AbstractString,
    message::AbstractString,
)::Nothing
    entity_value = Int64(entity)
    source_id = String(source_item_id)
    lock(cache.lock) do
        cache.failures[(string(entity_value), source_id)] = String(message)
        state = CachedCollectionResultState(
            Int8(kind), entity_value, Int8(RESULT_FAILED), source_id, String(message))
        cache.collection_result_states[(Int8(kind), entity_value)] = state
        _stage_ledger_failure!(
            cache.stage_ledger, (string(entity_value), source_id), true)
        _stage_ledger_result!(cache.stage_ledger, (Int8(kind), entity_value), state)
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
    _stage_ledger_failure!(cache.stage_ledger, (source_id, source_id), true)
    return nothing
end

function store_item_metadata!(cache::MemoryCacheDB, record::ItemRecord, ::AbstractDict)::Vector{String}
    lock(cache.lock) do
        state = CachedResultState(
            Int8(ITEM_ANALYSIS_RESULT),
            record.id,
            Int8(RESULT_READY),
            record.source_item_id,
            nothing,
        )
        cache.result_states[(Int8(ITEM_ANALYSIS_RESULT), record.id)] = state
        _stage_ledger_result!(
            cache.stage_ledger, (Int8(ITEM_ANALYSIS_RESULT), record.id), state)
    end
    return String[]
end

function store_collection_metadata!(cache::MemoryCacheDB, collection_key::Integer, ::AbstractDict)::Vector{String}
    entity = Int64(collection_key)
    lock(cache.lock) do
        state = CachedCollectionResultState(
            Int8(COLLECTION_ANALYSIS_RESULT), entity, Int8(RESULT_READY), "", nothing)
        cache.collection_result_states[(Int8(COLLECTION_ANALYSIS_RESULT), entity)] = state
        _stage_ledger_result!(
            cache.stage_ledger, (Int8(COLLECTION_ANALYSIS_RESULT), entity), state)
    end
    return String[]
end

function store_collection_process_result!(cache::MemoryCacheDB, collection_key::Integer)::Nothing
    entity = Int64(collection_key)
    lock(cache.lock) do
        state = CachedCollectionResultState(
            Int8(COLLECTION_PROCESS_RESULT), entity, Int8(RESULT_READY), "", nothing)
        cache.collection_result_states[(Int8(COLLECTION_PROCESS_RESULT), entity)] = state
        _stage_ledger_result!(
            cache.stage_ledger, (Int8(COLLECTION_PROCESS_RESULT), entity), state)
    end
    return nothing
end

"""Write one source item-metadata (entries) row."""
function edit_source_item_metadata!(cache::CacheDB, key::Int64, metadata::AbstractDict)::Nothing
    edit!(cache.source_item_metadata, key, metadata_dict(metadata))
    return nothing
end

edit_source_item_metadata!(::MemoryCacheDB, ::Int64, ::AbstractDict)::Nothing = nothing

"""Delete every cached result owned by one published source item, keyed by item surrogate."""
function delete_source_item!(
    cache::CacheDB,
    source_item_id::AbstractString,
    old_records::Vector{ItemRecord},
)::Nothing
    source_id = String(source_item_id)
    delete!(cache.source_items, source_id)
    _stage_ledger_source!(cache.stage_ledger, source_id, false)
    for record in old_records
        key = item_key!(cache, record.id)
        delete!(cache.items, record.id)
        _stage_ledger_item!(cache.stage_ledger, record.id, false)
        delete!(cache.source_item_metadata, key)
        delete!(cache.analyzed_item_metadata, key)
        delete!(cache.payload, (key, PAYLOAD_STAGE_PROCESSED))
        delete!(cache.payload, (key, PAYLOAD_STAGE_COLLECTION_PROCESSED))
        delete!(cache.interpreted, record.id)
        delete!(cache.processed_memory, record.id)
        delete!(cache.collection_processed_memory, record.id)
        delete!(cache.failures, (record.id, source_id))
        _stage_ledger_failure!(cache.stage_ledger, (record.id, source_id), false)
        delete!(cache.result_states, (Int8(PROCESSING_RESULT), record.id))
        _stage_ledger_result!(
            cache.stage_ledger, (Int8(PROCESSING_RESULT), record.id), nothing)
        delete!(cache.result_states, (Int8(ITEM_ANALYSIS_RESULT), record.id))
        _stage_ledger_result!(
            cache.stage_ledger, (Int8(ITEM_ANALYSIS_RESULT), record.id), nothing)
    end
    delete!(cache.failures, (source_id, source_id))
    _stage_ledger_failure!(cache.stage_ledger, (source_id, source_id), false)
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
    _stage_ledger_source!(cache.stage_ledger, source_id, false)
    _stage_ledger_failure!(cache.stage_ledger, (source_id, source_id), false)
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
        _stage_ledger_item!(cache.stage_ledger, record.id, false)
        _stage_ledger_failure!(cache.stage_ledger, (record.id, source_id), false)
        _stage_ledger_result!(
            cache.stage_ledger, (Int8(PROCESSING_RESULT), record.id), nothing)
        _stage_ledger_result!(
            cache.stage_ledger, (Int8(ITEM_ANALYSIS_RESULT), record.id), nothing)
    end
    return nothing
end

"""Delete collection records that the live index pruned after item removal or replacement."""
function delete_collection_records!(
    cache::CacheDB,
    collection_keys::Vector{Int64},
)::Nothing
    lock(cache.key_lock) do
        for collection_key in unique(collection_keys)
            delete!(cache.collections, collection_key)
            delete!(cache.persisted_collection_keys, collection_key)
        end
    end
    return nothing
end

delete_collection_records!(::MemoryCacheDB, ::Vector{Int64})::Nothing = nothing

"""Delete cached analysis for collections whose published metadata is invalid."""
function delete_collection_metadata!(
    cache::CacheDB,
    collection_keys::Vector{Int64},
)::Nothing
    for collection_key in unique(collection_keys)
        delete!(cache.analyzed_collection_metadata, collection_key)
        delete!(cache.collection_result_states, (Int8(COLLECTION_ANALYSIS_RESULT), collection_key))
        _stage_ledger_result!(
            cache.stage_ledger, (Int8(COLLECTION_ANALYSIS_RESULT), collection_key), nothing)
        delete!(cache.collection_result_states, (Int8(COLLECTION_PROCESS_RESULT), collection_key))
        _stage_ledger_result!(
            cache.stage_ledger, (Int8(COLLECTION_PROCESS_RESULT), collection_key), nothing)
    end
    return nothing
end

function delete_collection_metadata!(
    cache::MemoryCacheDB,
    collection_keys::Vector{Int64},
)::Nothing
    for collection_key in unique(collection_keys)
        lock(cache.lock) do
            delete!(cache.collection_result_states,
                (Int8(COLLECTION_ANALYSIS_RESULT), collection_key))
            delete!(cache.collection_result_states,
                (Int8(COLLECTION_PROCESS_RESULT), collection_key))
        end
        _stage_ledger_result!(
            cache.stage_ledger, (Int8(COLLECTION_ANALYSIS_RESULT), collection_key), nothing)
        _stage_ledger_result!(
            cache.stage_ledger, (Int8(COLLECTION_PROCESS_RESULT), collection_key), nothing)
    end
    return nothing
end

"""Aggregate staged durable keys and payload rows."""
function cache_pending_counts(cache::CacheDB)::NamedTuple{(:items, :rows),Tuple{Int,Int}}
    items = 0
    rows = 0
    for buffer in
        (cache.source_items, cache.collections, cache.items, cache.source_item_metadata,
         cache.analyzed_item_metadata,
         cache.analyzed_collection_metadata, cache.failures, cache.result_states,
         cache.collection_result_states, cache.payload)
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
function cache_stage_summary(cache::AbstractCacheDB)::CacheStageSummary
    return lock(cache.stage_ledger.lock) do
        cache.stage_ledger.summary
    end
end

"""Drop the cached interpretation ledger for one source item."""
function clear_cached_source_state!(cache::CacheDB, source_item_id::AbstractString)::Nothing
    source_id = String(source_item_id)
    delete!(cache.source_items, source_id)
    delete!(cache.failures, (source_id, source_id))
    _stage_ledger_source!(cache.stage_ledger, source_id, false)
    _stage_ledger_failure!(cache.stage_ledger, (source_id, source_id), false)
    return nothing
end

function clear_cached_source_state!(cache::MemoryCacheDB, source_item_id::AbstractString)::Nothing
    source_id = String(source_item_id)
    lock(cache.lock) do
        delete!(cache.source_items, source_id)
        delete!(cache.failures, (source_id, source_id))
    end
    _stage_ledger_source!(cache.stage_ledger, source_id, false)
    _stage_ledger_failure!(cache.stage_ledger, (source_id, source_id), false)
    return nothing
end

"""Drop one cached work-result ledger row."""
function clear_cached_result_state!(
    cache::CacheDB,
    kind::CacheResultKind,
    entity::AbstractString,
)::Nothing
    key = (Int8(kind), String(entity))
    delete!(cache.result_states, key)
    _stage_ledger_result!(cache.stage_ledger, key, nothing)
    return nothing
end

function clear_cached_result_state!(
    cache::MemoryCacheDB,
    kind::CacheResultKind,
    entity::AbstractString,
)::Nothing
    key = (Int8(kind), String(entity))
    lock(cache.lock) do
        delete!(cache.result_states, key)
    end
    _stage_ledger_result!(cache.stage_ledger, key, nothing)
    return nothing
end

function clear_cached_result_state!(
    cache::CacheDB,
    kind::CacheResultKind,
    entity::Integer,
)::Nothing
    key = (Int8(kind), Int64(entity))
    delete!(cache.collection_result_states, key)
    _stage_ledger_result!(cache.stage_ledger, key, nothing)
    return nothing
end

function clear_cached_result_state!(
    cache::MemoryCacheDB,
    kind::CacheResultKind,
    entity::Integer,
)::Nothing
    key = (Int8(kind), Int64(entity))
    lock(cache.lock) do
        delete!(cache.collection_result_states, key)
    end
    _stage_ledger_result!(cache.stage_ledger, key, nothing)
    return nothing
end

function _cache_stage_summary(source_items, items, failures, state_stores...)::CacheStageSummary
    ready(kind::CacheResultKind)::Int = count(
        state -> CacheResultKind(state.kind) === kind &&
            CacheResultStatus(state.status) === RESULT_READY,
        Iterators.flatten(values(store) for store in state_stores),
    )
    failed(kind::CacheResultKind)::Int = count(
        state -> CacheResultKind(state.kind) === kind &&
            CacheResultStatus(state.status) === RESULT_FAILED,
        Iterators.flatten(values(store) for store in state_stores),
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
        held = stage === :processed ?
            read(cache.processed_memory, record.id) :
            read(cache.collection_processed_memory, record.id)
        if held !== nothing
            loaded[index] = held
            continue
        end
        push!(disk_keys, (item_key!(cache, record.id), stage_code))
        loaded[index] = nothing
    end

    disk_data = isempty(disk_keys) ?
        Dict{Tuple{Int64,Int8},Any}() :
        read(cache.payload, unique(disk_keys))
    for (index, record) in pairs(records)
        loaded[index] !== nothing && continue
        data = get(disk_data, (item_key!(cache, record.id), stage_code), nothing)
        loaded[index] = data === nothing ? nothing : RegisteredDataItem(record, data)
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
        (cachedb.source_items, cachedb.collections, cachedb.items, cachedb.source_item_metadata,
         cachedb.analyzed_item_metadata,
         cachedb.analyzed_collection_metadata, cachedb.failures,
         cachedb.result_states, cachedb.collection_result_states, cachedb.payload,
         cachedb.interpreted,
         cachedb.processed_memory, cachedb.collection_processed_memory)
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
        Dict{CacheResultKey,AnyCachedResultState}(),
    )
end

"""Delete every index and item-data row and reset the item-key map, leaving a schema-valid cache."""
function clear_cache_index!(cachedb::CacheDB)::Nothing
    for buffer in
        (cachedb.source_items, cachedb.collections, cachedb.items, cachedb.source_item_metadata,
         cachedb.analyzed_item_metadata,
         cachedb.analyzed_collection_metadata, cachedb.failures,
         cachedb.result_states, cachedb.collection_result_states, cachedb.payload,
         cachedb.interpreted,
         cachedb.processed_memory, cachedb.collection_processed_memory)
        clear!(buffer)
    end
    lock(cachedb.key_lock) do
        empty!(cachedb.item_keys)
        cachedb.next_item_key = Int64(1)
        empty!(cachedb.persisted_collection_keys)
    end
    _reset_stage_ledger!(cachedb.stage_ledger)
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
        empty!(cachedb.collection_result_states)
    end
    _reset_stage_ledger!(cachedb.stage_ledger)
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

"""Rebuild the flat collection index from package-owned persisted records."""
function _load_collection_index(
    rows::Dict{Int64,CollectionRow},
    source_id::String,
    collection_analysis::Dict{Int64,MetadataDict},
)::CollectionIndex
    collections = CollectionIndex(source_id)
    loaded = Set{Int64}()
    visiting = Set{Int64}()

    function load_record(key::Int64)::Nothing
        key in loaded && return nothing
        key in visiting && error("Persisted collection hierarchy contains a parent cycle at key $key")
        row = get(rows, key, nothing)
        row === nothing && error("Persisted collection hierarchy is missing collection key $key")
        push!(visiting, key)
        row.parent_key === nothing || load_record(row.parent_key)
        own_metadata = _deserialize_hex(row.metadata_hex)
        own_metadata isa AbstractDict || error(
            "Persisted collection key $key metadata decoded as $(typeof(own_metadata)), not a dictionary")
        register_collection!(collections, CollectionRecord(
            row.collection_key,
            row.id,
            row.parent_key,
            row.label,
            metadata_dict(own_metadata),
            row.registration_name,
            get(collection_analysis, key, MetadataDict()),
        ))
        delete!(visiting, key)
        push!(loaded, key)
        return nothing
    end

    for key in keys(rows)
        load_record(key)
    end
    return collections
end

"""Reconstruct the full `SourceScan` stored in one DuckDB cache (entries-layer metadata only)."""
function _load_source_scan(
    cache::CacheDB,
    item_metadata_by_key::Dict{Int64,MetadataDict},
    collection_analysis::Dict{Int64,MetadataDict},
)::SourceScan
    identity = cache.identity
    _load_meta(cache)
    source_rows = read(cache.source_items)
    item_rows = read(cache.items)
    failure_rows = read(cache.failures)
    collections = try
        _load_collection_index(
            read(cache.collections), identity.source_id, collection_analysis)
    catch error
        throw(ProjectCacheDataError(
            identity.cache_path,
            "cached collection records are incompatible: $(sprint(showerror, error))",
        ))
    end

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
        if row.collection_key !== nothing && !haskey(collections.records, row.collection_key)
            throw(ProjectCacheDataError(
                identity.cache_path,
                "item '$(row.id)' refers to missing collection key $(row.collection_key)",
            ))
        end
        push!(records, ItemRecord(;
            id=row.id,
            source_item_id=row.source_item_id,
            source_item_path=path,
            source_item_timestamp=timestamp,
            label=row.label,
            kind=Symbol(row.kind),
            collection_key=row.collection_key,
            metadata=get(item_metadata_by_key, row.item_key, MetadataDict()),
        ))
        append_item!(collections, row.id, row.collection_key)
    end

    sort!(records; by=record -> (record.source_item_id, record.id))
    # Scan-wide summaries are derived, not stored: a source item with no items was skipped or failed.
    source_ids_with_items = Set(record.source_item_id for record in records)
    skipped_count = count(id -> !(id in source_ids_with_items), all_source_ids)
    collections.skipped_count = skipped_count

    return SourceScan(
        identity.source_id,
        identity.source_label,
        collections,
        records,
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
            item_id = get(key_to_id, key, nothing)
            item_id === nothing && continue
            entries = get(entries_by_key, key, MetadataDict())
            computed = MetadataDict()
            for (name, value) in analyzed
                get(entries, name, nothing) == value || (computed[name] = value)
            end
            isempty(computed) || (item_metadata[item_id] = computed)
        end
        base = ProjectCacheIndex(
            identity, _load_source_scan(cachedb, entries_by_key, collection_analysis))
        result_states = Dict{CacheResultKey,AnyCachedResultState}()
        for state in values(read(cachedb.result_states))
            result_states[CacheResultKey(CacheResultKind(state.kind), state.entity)] = state
        end
        for state in values(read(cachedb.collection_result_states))
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
        loaded_items=length(index.source.items),
        skipped_source_items=index.source.collections.skipped_count,
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
