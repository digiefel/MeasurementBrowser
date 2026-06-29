const PROJECT_CACHE_SCHEMA_VERSION = 6
const ITEM_DATA_VIEW = "__measurementbrowser_item_data"

"""
DuckDB buffer-pool limit (MiB) for cache connections.

DuckDB otherwise defaults to most of system RAM and retains committed table blocks for the workspace's
lifetime; cache writes are bounded batches, so a smaller buffer keeps large caches out of swap without
reducing scan parallelism. Settable at runtime — e.g. from a GUI slider — via
[`set_cache_memory_limit!`](@ref). Changes apply to caches opened afterward; pass an open
[`CacheDB`](@ref) to also resize it live.
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

"""A cacheable item whose concrete item/data types have no native cache implementation."""
struct UnsupportedItemDataCacheError <: Exception
    item_id::String
    item_type::DataType
    data_type::DataType
end

"""Explain which item-data combination needs a cache implementation or explicit opt-out."""
function Base.showerror(io::IO, err::UnsupportedItemDataCacheError)::Nothing
    print(
        io,
        "Cannot persist item data for '$(err.item_id)': item type $(err.item_type) carries " *
        "$(err.data_type), but the native cache currently supports DataItem values carrying " *
        "AbstractDataFrame data. Define cacheable(::$(err.item_type)) = false to opt this item " *
        "type out.",
    )
end

"""The project, source, and DuckDB file belonging to one cache."""
struct ProjectCacheIdentity
    project_name::String
    source_id::String
    source_label::String
    cache_path::String
end

"""Counts describing the differences between a source scan and its cache."""
struct ProjectCacheStatus
    total_source_items::Int
    cached_source_items::Int
    fresh_source_items::Int
    stale_source_items::Int
    new_source_items::Int
    deleted_source_items::Int
    error_source_items::Int
end

"""All data-less content needed to restore and compare a project cache."""
struct ProjectCacheIndex
    identity::ProjectCacheIdentity
    source::SourceScan
    analysis_errors::Dict{String,String}
end

# --------------------------------------------------------------------------------------------------
# Identity
# --------------------------------------------------------------------------------------------------

"""
Bind one project name to a source and its package-owned cache path.

The project name is used exactly as one directory component. Invalid names fail instead of being
silently rewritten into a different cache identity.
"""
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

# --------------------------------------------------------------------------------------------------
# Connection + schema
# --------------------------------------------------------------------------------------------------

"""
One open DuckDB cache file shared by all of a workspace's cache work.

A single [`DuckDB.DB`](@ref) holds the file for the workspace's lifetime. One locked writer keeps
surrogate-key order aligned with physical append order. Reads use short-lived connections so
independent interactive and background reads can run concurrently against committed snapshots.
"""
mutable struct CacheDB
    identity::ProjectCacheIdentity
    db::DuckDB.DB
    # The per-table buffer is the single write/read door; ProjectCache owns and drives it.
    buffer::CacheBuffer
    writer::DuckDB.Connection
    writer_lock::ReentrantLock
    profiler::Profiling.ProfileSession
    # Time threads spend waiting for a free writer vs. holding one doing work.
    writer_wait_ns::Base.Threads.Atomic{Int64}
    writer_busy_ns::Base.Threads.Atomic{Int64}
end

"""
Open (creating if needed) the cache file for one workspace and ensure its schema.

A cache that cannot be opened or whose schema cannot be created — a truncated file, a stale
write-ahead log from a hard kill, or an unreadable older layout — is discarded and rebuilt from
scratch rather than bricking the workspace. The data it held is recoverable by rescanning.
"""
function open_cache_db(
    identity::ProjectCacheIdentity,
    profiler::Profiling.ProfileSession,
    metrics::BuildMetrics=BuildMetrics(),
)::CacheDB
    mkpath(dirname(identity.cache_path))
    try
        return _connect_cache_db(identity, profiler, metrics)
    catch err
        err isa InterruptException && rethrow()
        @warn(
            "Discarding an incompatible or unreadable project cache and rebuilding it",
            cache_path=identity.cache_path,
            exception=err,
        )
        _remove_cache_files(identity.cache_path)
        return _connect_cache_db(identity, profiler, metrics)
    end
end

"""Open a cache with internal profiling disabled for direct cache operations and tests."""
function open_cache_db(identity::ProjectCacheIdentity)::CacheDB
    return open_cache_db(identity, Profiling.ProfileSession(false, false, nothing, nothing))
end

"""Connect to the cache file, ensure its schema, and wrap it in a [`CacheDB`](@ref)."""
function _connect_cache_db(
    identity::ProjectCacheIdentity,
    profiler::Profiling.ProfileSession,
    metrics::BuildMetrics,
)::CacheDB
    db = DBInterface.connect(DuckDB.DB, identity.cache_path)
    try
        DBInterface.execute(db, _memory_limit_sql(CACHE_MEMORY_LIMIT_MIB[]))   # see CACHE_MEMORY_LIMIT_MIB
        ensure_schema!(db)
        item_data_columns = Set(String(row.name) for row in
            DBInterface.execute(db, "PRAGMA table_info('item_data')"))
        "stage" in item_data_columns ||
            throw(ProjectCacheError(identity.cache_path, "cache schema is out of date"))
        schema_versions = String[
            String(row.value)
            for row in DBInterface.execute(
                db,
                "SELECT value FROM meta WHERE key = 'schema_version'",
            )
        ]
        if !isempty(schema_versions) &&
           only(schema_versions) != string(PROJECT_CACHE_SCHEMA_VERSION)
            throw(ProjectCacheError(identity.cache_path, "cache schema is out of date"))
        end
        return CacheDB(identity, db, CacheBuffer(db, metrics), DBInterface.connect(db),
            ReentrantLock(), profiler,
            Base.Threads.Atomic{Int64}(0), Base.Threads.Atomic{Int64}(0))
    catch
        DBInterface.close!(db)
        rethrow()
    end
end

"""Delete a cache file and its DuckDB sidecars (write-ahead log, temp dir)."""
function _remove_cache_files(cache_path::AbstractString)::Nothing
    for path in (cache_path, cache_path * ".wal", cache_path * ".tmp")
        rm(path; force=true, recursive=true)
    end
    return nothing
end

# --------------------------------------------------------------------------------------------------
# The domain write/read surface over the per-table buffer.
#
# These map the pipeline stages (interpreted/processed/failure) and `ItemRecord`s onto the buffer's
# generic rows. The CacheBuffer below knows none of this; it just buffers typed rows per table.
# --------------------------------------------------------------------------------------------------

"""
Store one source item's interpreted result: data-less records to their disk-backed buffers, and the
interpreted payloads to the memory-only buffer.

Item statistics are not written here — they are computed during processing and stored by
[`store_processed!`](@ref), so each metadata key is written exactly once per build (pure append).
"""
function store_interpreted!(
    cache::CacheDB,
    records::Vector{ItemRecord},
    data::Vector{Any},
)::Nothing
    buffer = cache.buffer
    started = time_ns()
    for record in records
        hex, hash = _encode_fingerprint(record.source_item_fingerprint)
        _store!(buffer.source_items, record.source_item_id,
            SourceItemRow(record.source_item_id, hex, hash,
                record.source_item_path, record.source_item_timestamp))
        _store!(buffer.items, record.id,
            ItemRow(record.id, record.source_item_id, record.item_label, String(record.kind),
                record.collection, _serialize_hex(record.item_fingerprint)))
        for (key, value) in record.parameters
            _store!(buffer.metadata, (Int8(SCOPE_ITEM_PARAMETERS), record.id, String(key)),
                MetaRow(Int8(SCOPE_ITEM_PARAMETERS), record.id, String(key), value))
        end
    end
    for (record, item) in zip(records, data)
        item isa AbstractDataItem || continue
        _store!(buffer.interpreted, record.id, item)
    end
    record_cache_phase!(buffer.metrics.interpreted_write_ns, buffer.metrics.interpreted_writes, started)
    return nothing
end

"""
Store one processed result: its statistics to the metadata buffer, and the processed item itself to
either the disk-backed read-through payload buffer (when `cacheable` and the payload is a non-empty
DataFrame) or the memory-only `processed_memory` buffer otherwise.

A non-cacheable (or non-tabular) processed item is therefore held exactly like an interpreted one —
bounded, never written, recomputed from source on a later miss — so a re-selection still skips
reprocessing while the item is resident (docs/cache.md "Two backing modes").
"""
function store_processed!(
    cache::CacheDB,
    record::ItemRecord,
    item::Union{Nothing,AbstractDataItem},
    stats::Union{Nothing,MetadataDict},
    cacheable::Bool,
)::Nothing
    buffer = cache.buffer
    if stats !== nothing
        started = time_ns()
        for (key, value) in stats
            _store!(buffer.metadata, (Int8(SCOPE_ITEM_STATS), record.id, String(key)),
                MetaRow(Int8(SCOPE_ITEM_STATS), record.id, String(key), value))
        end
        record_cache_phase!(buffer.metrics.stats_write_ns, buffer.metrics.stats_writes, started)
    end
    item === nothing && return nothing
    payload = item_data(item)
    disk = cacheable && record.source_item_fingerprint !== nothing &&
        payload isa AbstractDataFrame && !isempty(names(payload))
    if disk
        started = time_ns()
        _store!(buffer.payload, record.id,
            PayloadEntry(record, payload, _dataframe_storage_id(payload, :processed),
                names(payload)))
        record_cache_phase!(
            buffer.metrics.processed_write_ns, buffer.metrics.processed_writes, started)
    else
        _store!(buffer.processed_memory, record.id, item)
    end
    return nothing
end

"""Store one item's processing/statistics failure."""
function store_failure!(cache::CacheDB, record::ItemRecord, message::String)::Nothing
    _store!(cache.buffer.failures, record.id,
        FailureRow(record.id, record.source_item_id, message))
    return nothing
end

"""
Store item and collection statistics through the metadata buffer.

Item statistics are scoped `SCOPE_ITEM_STATS` keyed by item id; collection statistics are
`SCOPE_NODE_STATS` keyed by the node's collection path. Like every other write this is a pure append
that coalesces per `(scope, entity, key)`; re-running analysis is picked up by a full Rebuild, not an
in-place replace (docs/cache.md "Mutations, deferred").
"""
function store_stats!(
    cache::CacheDB,
    item_stats::AbstractDict{String,<:AbstractDict},
    node_stats::AbstractDict{<:Tuple,<:AbstractDict},
)::Nothing
    (isempty(item_stats) && isempty(node_stats)) && return nothing
    metadata = cache.buffer.metadata
    for (entity, stats) in item_stats
        id = String(entity)
        for (key, value) in metadata_dict(stats)
            _store!(metadata, (Int8(SCOPE_ITEM_STATS), id, String(key)),
                MetaRow(Int8(SCOPE_ITEM_STATS), id, String(key), value))
        end
    end
    for (path, stats) in node_stats
        entity = collection_path_key(collect(String, path))
        for (key, value) in metadata_dict(stats)
            _store!(metadata, (Int8(SCOPE_NODE_STATS), entity, String(key)),
                MetaRow(Int8(SCOPE_NODE_STATS), entity, String(key), value))
        end
    end
    return nothing
end

"""
Read item data back through the buffer: memory first, then (for processed) the database.

Returns a vector aligned to `records`, each element a loaded item or `nothing` on a miss; the caller
falls back to the source. Interpreted payloads are memory-only, so a memory miss is a final miss.
"""
function read_item_data(
    cache::CacheDB,
    records::Vector{ItemRecord};
    stage::Symbol,
)::Vector{Any}
    stage in (:interpreted, :processed) ||
        throw(ArgumentError("unknown item-data cache stage '$stage'"))
    isempty(records) && return Any[]
    buffer = cache.buffer

    if stage === :interpreted
        return lock(buffer.interpreted.cond) do
            Any[get(buffer.interpreted.entries, record.id, nothing) for record in records]
        end
    end

    loaded = Vector{Any}(undef, length(records))
    misses = ItemRecord[]
    miss_positions = Int[]
    lock(buffer.payload.cond) do
        for (position, record) in pairs(records)
            entry = get(buffer.payload.entries, record.id, nothing)
            if entry !== nothing && entry.record.item_fingerprint == record.item_fingerprint
                loaded[position] = DataItem(record, entry.data)
            else
                loaded[position] = nothing
                push!(misses, record)
                push!(miss_positions, position)
            end
        end
    end
    # Non-cacheable processed items never reach disk; serve any still resident from memory before the
    # disk read narrows to the records that could actually have a pointer.
    if !isempty(misses)
        disk_misses = ItemRecord[]
        disk_positions = Int[]
        lock(buffer.processed_memory.cond) do
            for (record, position) in zip(misses, miss_positions)
                held = get(buffer.processed_memory.entries, record.id, nothing)
                if held !== nothing
                    loaded[position] = held
                else
                    push!(disk_misses, record)
                    push!(disk_positions, position)
                end
            end
        end
        misses, miss_positions = disk_misses, disk_positions
    end
    if !isempty(misses)
        fetched = _read_payloads_from_disk(buffer.payload, misses)
        for (index, position) in pairs(miss_positions)
            loaded[position] = fetched[index]
        end
    end
    return loaded
end

"""Close a workspace cache's connections and the underlying database file."""
function close_cache_db!(cachedb::CacheDB)::Nothing
    # Fold the write-ahead log into the database file so the next workspace opened on this path within
    # the same process reads the committed state instead of a pre-commit snapshot.
    try
        DBInterface.execute(cachedb.writer, "CHECKPOINT")
    finally
        DBInterface.close!(cachedb.writer)
        DBInterface.close!(cachedb.db)
    end
    return nothing
end

"""Set the cache buffer-pool limit (MiB) and apply it immediately to one open cache."""
function set_cache_memory_limit!(cachedb::CacheDB, mib::Integer)::Int
    set_cache_memory_limit!(mib)
    DBInterface.execute(cachedb.db, _memory_limit_sql(mib))
    return Int(mib)
end

"""Run `work` on a fresh read connection that snapshots the latest committed cache state."""
function with_reader(work::Function, cachedb::CacheDB)::Any
    return @profile_span cachedb.profiler :cache :reader ProfileAttributes() begin
        connection = @profile_span cachedb.profiler :cache :connect_reader ProfileAttributes() begin
            DBInterface.connect(cachedb.db)
        end
        try
            work(connection)
        finally
            DBInterface.close!(connection)
        end
    end
end

"""Run a multi-statement cache read against one committed snapshot."""
function with_reader_snapshot(work::Function, cachedb::CacheDB)::Any
    return with_reader(cachedb) do connection
        DBInterface.execute(connection, "BEGIN TRANSACTION")
        try
            result = work(connection)
            DBInterface.execute(connection, "COMMIT")
            result
        catch error
            try
                DBInterface.execute(connection, "ROLLBACK")
            catch rollback_error
                throw(CompositeException(error, rollback_error))
            end
            rethrow()
        end
    end
end

"""Run `work` on the cache's serialized writer connection."""
function with_writer(work::Function, cachedb::CacheDB)::Any
    profiler = cachedb.profiler
    token = Profiling.should_trace(profiler) ?
        Profiling.start_span!(profiler, :cache, :writer) : nothing
    requested = time_ns()
    lock(cachedb.writer_lock)
    acquired = time_ns()
    status = :ok
    Base.Threads.atomic_add!(cachedb.writer_wait_ns, Int64(acquired - requested))
    try
        if token === nothing
            return work(cachedb.writer)
        end
        return Base.ScopedValues.with(
            Profiling.CURRENT_SPAN => token.id,
            Profiling.CURRENT_SESSION => profiler,
        ) do
            work(cachedb.writer)
        end
    catch
        status = :error
        rethrow()
    finally
        released = time_ns()
        service_ns = Int64(released - acquired)
        Base.Threads.atomic_add!(cachedb.writer_busy_ns, service_ns)
        unlock(cachedb.writer_lock)
        token === nothing || Profiling.finish_span!(token; status, attributes=ProfileAttributes(
            wait_ns=Int64(acquired - requested),
            service_ns=service_ns,
        ))
    end
end

"""
Run one cache mutation as a writer transaction.

Large scans call this once per bounded source-item batch so DuckDB can release transaction-owned
write memory between batches.
"""
function with_writer_transaction(work::Function, cachedb::CacheDB)::Any
    return @profile_span cachedb.profiler :cache :transaction ProfileAttributes() begin
        with_writer(cachedb) do connection
            @profile_span cachedb.profiler :cache :begin_transaction ProfileAttributes() begin
                DBInterface.execute(connection, "BEGIN TRANSACTION")
            end
            try
                result = work(connection)
                @profile_span cachedb.profiler :cache :commit ProfileAttributes() begin
                    DBInterface.execute(connection, "COMMIT")
                end
                result
            catch transaction_error
                try
                    DBInterface.execute(connection, "ROLLBACK")
                catch rollback_error
                    throw(CompositeException(transaction_error, rollback_error))
                end
                throw(transaction_error)
            end
        end
    end
end

"""Create the cache tables if they do not already exist."""
function ensure_schema!(connection)::Nothing
    DBInterface.execute(connection, """
        CREATE TABLE IF NOT EXISTS meta(
            key TEXT PRIMARY KEY,
            value TEXT)
    """)
    DBInterface.execute(connection, """
        CREATE TABLE IF NOT EXISTS source_items(
            source_item_id TEXT PRIMARY KEY,
            fingerprint TEXT,
            fingerprint_hash TEXT,
            path TEXT,
            timestamp TIMESTAMP,
            error TEXT)
    """)
    DBInterface.execute(connection, """
        CREATE TABLE IF NOT EXISTS items(
            id TEXT PRIMARY KEY,
            source_item_id TEXT,
            item_label TEXT,
            kind TEXT,
            collection VARCHAR[],
            item_fingerprint TEXT)
    """)
    DBInterface.execute(connection, """
        CREATE TABLE IF NOT EXISTS metadata(
            scope TINYINT,
            entity TEXT,
            key TEXT,
            vtype TINYINT,
            b BOOLEAN, i BIGINT, d DOUBLE, s TEXT, ts TIMESTAMP,
            lb BOOLEAN[], li BIGINT[], ld DOUBLE[], ls VARCHAR[],
            PRIMARY KEY (scope, entity, key))
    """)
    # Payload rows use this compact surrogate instead of repeating long item ids. Fresh values on every
    # write also let replacement deletes distinguish current rows from stale ones.
    DBInterface.execute(connection, "CREATE SEQUENCE IF NOT EXISTS item_data_seq START 1")
    DBInterface.execute(connection, """
        CREATE TABLE IF NOT EXISTS item_data(
            item_id TEXT,
            stage TEXT,
            source_item_id TEXT,
            sif_hash TEXT,
            if_hash TEXT,
            storage_id TEXT,
            row_count BIGINT,
            seq BIGINT,
            PRIMARY KEY (item_id, stage))
    """)
    DBInterface.execute(connection, """
        CREATE TABLE IF NOT EXISTS dataframe_schemas(
            storage_id TEXT PRIMARY KEY,
            column_names VARCHAR[])
    """)
    DBInterface.execute(connection, """
        CREATE TABLE IF NOT EXISTS item_failures(
            item_id TEXT PRIMARY KEY,
            source_item_id TEXT,
            message TEXT)
    """)
    return nothing
end

# --------------------------------------------------------------------------------------------------
# Encoding helpers
# --------------------------------------------------------------------------------------------------

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

"""Content hash of a fingerprint, or `nothing` when the fingerprint is absent."""
_fingerprint_hash(fingerprint) =
    fingerprint === nothing ? nothing : bytes2hex(sha1(_serialize_bytes(fingerprint)))

"""Hex-text storage and content hash for one fingerprint, as `(hex, hash)` (`nothing` when absent)."""
function _encode_fingerprint(fingerprint)::Tuple{Union{Nothing,String},Union{Nothing,String}}
    fingerprint === nothing && return (nothing, nothing)
    bytes = _serialize_bytes(fingerprint)
    return (bytes2hex(bytes), bytes2hex(sha1(bytes)))
end

"""Bulk-append positional value tuples to `table` through one appender."""
function _append_rows!(connection, table::AbstractString, rows)::Nothing
    isempty(rows) && return nothing
    appender = DuckDB.Appender(connection, table)
    try
        for row in rows
            for value in row
                DuckDB.append(appender, value)
            end
            DuckDB.end_row(appender)
        end
        DuckDB.flush(appender)
    finally
        DuckDB.close(appender)
    end
    return nothing
end

"""
How one [`MetadataValue`](@ref) variant maps to the EAV `metadata` table: its discriminator, the
value column that stores it, and the converters between the Julia value and that column's storage.
"""
struct MetaCodec
    vtype::MetaVType
    column::Symbol
    to_db::Any          # Julia value -> stored column value (a callable; may be a type like `String`)
    from_db::Any        # stored column value -> Julia value
end

# Single source of truth for the EAV mapping. Adding a `MetadataValue` variant means adding one row
# here (and a value column below if it needs a new storage type); encode, decode, and the row layout
# are all derived from this table.
const META_CODECS = (
    MetaCodec(VT_BOOL,     :b,  identity, identity),
    MetaCodec(VT_INT,      :i,  identity, identity),
    MetaCodec(VT_FLOAT,    :d,  identity, identity),
    MetaCodec(VT_STRING,   :s,  identity, identity),
    MetaCodec(VT_SYMBOL,   :s,  String,            Symbol),
    MetaCodec(VT_DATE,     :ts, DateTime,          Date),
    MetaCodec(VT_DATETIME, :ts, identity,          identity),
    MetaCodec(VT_MISSING,  :s,  _ -> missing,      _ -> missing),
    MetaCodec(VT_VBOOL,    :lb, identity, value -> Vector{Bool}(value)),
    MetaCodec(VT_VINT,     :li, identity, value -> Vector{Int64}(value)),
    MetaCodec(VT_VFLOAT,   :ld, identity, value -> Vector{Float64}(value)),
    MetaCodec(VT_VSTRING,  :ls, identity, value -> Vector{String}(value)),
)

# The value columns of the `metadata` table, in the order the appender writes them.
const META_VALUE_COLUMNS = (:b, :i, :d, :s, :ts, :lb, :li, :ld, :ls)

# `MetadataValue` is a small closed union; map each concrete member type to its codec once.
const _META_CODEC_BY_TYPE = Dict{DataType,MetaCodec}(
    T => codec
    for (T, codec) in zip(
        (Bool, Int64, Float64, String, Symbol, Date, DateTime, Missing,
         Vector{Bool}, Vector{Int64}, Vector{Float64}, Vector{String}),
        META_CODECS,
    )
)
const _META_CODEC_BY_VTYPE = Dict(codec.vtype => codec for codec in META_CODECS)

_meta_codec(value)::MetaCodec = _META_CODEC_BY_TYPE[typeof(value)]
_meta_codec(vtype::MetaVType)::MetaCodec = _META_CODEC_BY_VTYPE[vtype]

"""Reconstruct a [`MetadataValue`](@ref) from one `metadata` row given its `vtype` discriminator."""
function _decode_metadata_value(vtype_raw, row)::MetadataValue
    codec = _meta_codec(MetaVType(Int8(vtype_raw)))
    return codec.from_db(getproperty(row, codec.column))
end

# --------------------------------------------------------------------------------------------------
# Index construction + comparison (in memory; unchanged from the previous cache)
# --------------------------------------------------------------------------------------------------

"""
Build the complete cache index represented by one finished source scan.

Analysis failures remain attached to their source items while successful items stay available.
"""
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
    )
end

# --------------------------------------------------------------------------------------------------
# Writing the index
# --------------------------------------------------------------------------------------------------

"""Append every key of one metadata dict as EAV rows under `scope`/`entity` via the appender."""
function _append_metadata!(appender, scope::MetaScope, entity::AbstractString, dict)::Nothing
    for (key, value) in dict
        codec = _meta_codec(value)
        stored = codec.to_db(value)
        DuckDB.append(appender, Int8(scope))
        DuckDB.append(appender, String(entity))
        DuckDB.append(appender, String(key))
        DuckDB.append(appender, Int8(codec.vtype))
        for column in META_VALUE_COLUMNS
            DuckDB.append(appender, column === codec.column ? stored : missing)
        end
        DuckDB.end_row(appender)
    end
    return nothing
end

"""Write the identity/version `meta` rows describing one source scan."""
function _write_meta_rows!(connection, identity::ProjectCacheIdentity, source::SourceScan)::Nothing
    statement = DBInterface.prepare(connection,
        "INSERT INTO meta VALUES (?, ?) ON CONFLICT (key) DO UPDATE SET value = excluded.value")
    for (key, value) in (
        "schema_version" => string(PROJECT_CACHE_SCHEMA_VERSION),
        "project_name" => identity.project_name,
        "source_id" => identity.source_id,
        "source_label" => identity.source_label,
        "skipped_count" => string(source.hierarchy.skipped_count),
        "has_collection_parameters" => string(source.hierarchy.has_collection_parameters),
    )
        DBInterface.execute(statement, (key, value))
    end
    return nothing
end

"""Map each source-item id to one representative `(path, timestamp)` taken from its records."""
function _source_item_locations(source::SourceScan)::Dict{String,Tuple{Union{Nothing,String},Union{Nothing,DateTime}}}
    locations = Dict{String,Tuple{Union{Nothing,String},Union{Nothing,DateTime}}}()
    for record in source.hierarchy.all_items
        haskey(locations, record.source_item_id) && continue
        locations[record.source_item_id] =
            (record.source_item_path, record.source_item_timestamp)
    end
    return locations
end

# --------------------------------------------------------------------------------------------------
# Incremental writing (per source item, as the scan streams)
# --------------------------------------------------------------------------------------------------

"""Delete every index and item-data row, leaving an empty (but schema-valid) cache for a rebuild."""
function clear_cache_index!(cachedb::CacheDB)::Nothing
    with_writer(cachedb) do connection
        DBInterface.execute(connection, "BEGIN TRANSACTION")
        try
            _delete_all_cached_item_data!(connection)
            for table in ("item_failures", "metadata", "items", "source_items", "meta")
                DBInterface.execute(connection, "DELETE FROM $table")
            end
            DBInterface.execute(connection, "COMMIT")
        catch
            DBInterface.execute(connection, "ROLLBACK")
            rethrow()
        end
    end
    return nothing
end

"""
Stamp the identity `meta` rows so a scan's incremental writes are loadable before it finishes.

Without these rows [`load_cache_index`](@ref) treats the cache as unbuilt, so an interrupted scan
would discard the per-item progress already written. The scan-wide `skipped_count` and
`has_collection_parameters` are filled in later by [`finalize_scan!`](@ref).
"""
function write_scan_identity!(cachedb::CacheDB)::Nothing
    identity = cachedb.identity
    with_writer(cachedb) do connection
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
    end
    return nothing
end

"""Delete the index rows and, when invalidated, cached data belonging to one source item."""
function _delete_source_item_rows!(
    connection,
    source_item_id::AbstractString;
    drop_item_data::Bool,
)::Nothing
    drop_item_data && _delete_cached_source_item_data!(connection, source_item_id)
    DBInterface.execute(
        DBInterface.prepare(connection, "DELETE FROM item_failures WHERE source_item_id = ?"),
        (source_item_id,))
    DBInterface.execute(
        DBInterface.prepare(connection, """
            DELETE FROM metadata WHERE scope IN (?, ?)
              AND entity IN (SELECT id FROM items WHERE source_item_id = ?)"""),
        (Int8(SCOPE_ITEM_PARAMETERS), Int8(SCOPE_ITEM_STATS), source_item_id))
    DBInterface.execute(
        DBInterface.prepare(connection, "DELETE FROM items WHERE source_item_id = ?"),
        (source_item_id,))
    DBInterface.execute(
        DBInterface.prepare(connection, "DELETE FROM source_items WHERE source_item_id = ?"),
        (source_item_id,))
    return nothing
end

"""
Finish an incremental scan: record bookkeeping the per-item writes can't carry.

`written_ids` are the source items already persisted by the per-table buffers. This writes the `meta`
rows and collection-node parameters, records bare rows for skipped/failed source items, stamps
per-source-item analysis errors, and deletes source items (cascading their items, metadata, and cached
item data) that the current scan no longer contains.
"""
function finalize_scan!(
    connection,
    identity::ProjectCacheIdentity,
    source::SourceScan,
    written_ids::AbstractSet{<:AbstractString},
)::Nothing
    new_index = ProjectCacheIndex(identity, source)
    locations = _source_item_locations(source)
    current_ids = keys(source.source_item_fingerprints)

    begin
        _write_meta_rows!(connection, identity, source)

        # Bare rows for source items the scan produced no records for (skipped or failed).
        bare_rows = NTuple{6,Any}[]
        for (id, fingerprint) in source.source_item_fingerprints
            id in written_ids && continue
            path, timestamp = get(locations, id, (nothing, nothing))
            hex, hash = _encode_fingerprint(fingerprint)
            push!(bare_rows,
                (id, hex, hash, path, timestamp, get(new_index.analysis_errors, id, nothing)))
        end
        _append_rows!(connection, "source_items", bare_rows)

        # Stamp errors onto source items that were written but later failed analysis interpretation.
        error_stmt = DBInterface.prepare(connection,
            "UPDATE source_items SET error = ? WHERE source_item_id = ?")
        for (id, message) in new_index.analysis_errors
            id in written_ids || continue
            DBInterface.execute(error_stmt, (message, id))
        end

        # Collection-node parameters (node stats arrive later from analysis): one delete, then a
        # single appender pass over every node's metadata.
        node_entities = String[
            collection_path_key(collect(String, path)) for path in keys(source.hierarchy.index)]
        if !isempty(node_entities)
            placeholders = join(fill("?", length(node_entities)), ", ")
            DBInterface.execute(
                DBInterface.prepare(connection,
                    "DELETE FROM metadata WHERE scope = ? AND entity IN ($placeholders)"),
                (Int8(SCOPE_NODE_PARAMETERS), node_entities...))
        end
        node_appender = DuckDB.Appender(connection, "metadata")
        try
            for (path, node) in source.hierarchy.index
                isempty(node.parameters) && continue
                _append_metadata!(node_appender, SCOPE_NODE_PARAMETERS,
                    collection_path_key(collect(String, path)), node.parameters)
            end
            DuckDB.flush(node_appender)
        finally
            DuckDB.close(node_appender)
        end

        # Drop source items (and everything they own) that the scan no longer contains.
        stored_ids = String[
            row.source_item_id
            for row in DBInterface.execute(connection, "SELECT source_item_id FROM source_items")
        ]
        for id in stored_ids
            id in current_ids && continue
            _delete_source_item_rows!(connection, id; drop_item_data=true)
        end
    end
    return nothing
end

"""Finish one scan in its own writer transaction."""
function finalize_scan!(
    cachedb::CacheDB,
    source::SourceScan,
    written_ids::AbstractSet{<:AbstractString},
)::Nothing
    return @profile_span cachedb.profiler :cache :finalize_scan ProfileAttributes(
        source_id=cachedb.identity.source_id,
        items=Int64(length(source.hierarchy.all_items)),
    ) begin
        with_writer_transaction(cachedb) do connection
            finalize_scan!(connection, cachedb.identity, source, written_ids)
        end
    end
end

"""Replace or clear one processing/statistics failure on an open writer connection."""
function _persist_item_failure!(
    connection,
    record::ItemRecord,
    message::Union{Nothing,String},
)::Nothing
    statement = if message === nothing
        DBInterface.prepare(connection, "DELETE FROM item_failures WHERE item_id = ?")
    else
        DBInterface.prepare(connection, """
            INSERT INTO item_failures VALUES (?, ?, ?)
            ON CONFLICT (item_id) DO UPDATE SET
                source_item_id = excluded.source_item_id,
                message = excluded.message
        """)
    end
    values = message === nothing ? (record.id,) :
        (record.id, record.source_item_id, message)
    DBInterface.execute(statement, values)
    return nothing
end


# --------------------------------------------------------------------------------------------------
# Loading the index
# --------------------------------------------------------------------------------------------------

"""Read the `meta` table into a `Dict{String,String}`, validating identity and schema."""
function _load_meta(connection, identity::ProjectCacheIdentity)::Dict{String,String}
    meta = Dict{String,String}()
    for row in DBInterface.execute(connection, "SELECT key, value FROM meta")
        meta[row.key] = row.value
    end
    isempty(meta) &&
        throw(ProjectCacheError(identity.cache_path, "cache has not been built yet"))
    get(meta, "schema_version", "") == string(PROJECT_CACHE_SCHEMA_VERSION) ||
        throw(ProjectCacheError(identity.cache_path, "cache schema is out of date"))
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

"""Read every EAV row, grouped into `(scope, entity) => MetadataDict`."""
function _load_metadata(connection)::Dict{Tuple{MetaScope,String},MetadataDict}
    grouped = Dict{Tuple{MetaScope,String},MetadataDict}()
    for row in DBInterface.execute(connection, """
        SELECT scope, entity, key, vtype, b, i, d, s, ts, lb, li, ld, ls FROM metadata
    """)
        scope = MetaScope(Int8(row.scope))
        dict = get!(() -> MetadataDict(), grouped, (scope, row.entity))
        dict[Symbol(row.key)] = _decode_metadata_value(row.vtype, row)
    end
    return grouped
end

"""Reconstruct the full `SourceScan` stored in one DuckDB cache."""
function _load_source_scan(connection, identity::ProjectCacheIdentity)::SourceScan
    meta = _load_meta(connection, identity)
    skipped_count = parse(Int, get(meta, "skipped_count", "0"))
    has_collection_parameters = get(meta, "has_collection_parameters", "false") == "true"

    fingerprints = Dict{String,Any}()
    locations = Dict{String,Tuple{Union{Nothing,String},Union{Nothing,DateTime}}}()
    failures = ItemFailure[]
    for row in DBInterface.execute(connection, """
        SELECT source_item_id, fingerprint, path, timestamp, error FROM source_items
    """)
        id = row.source_item_id
        fingerprints[id] = _deserialize_hex(row.fingerprint)
        locations[id] = (_null_to_nothing(row.path), _null_to_nothing(row.timestamp))
        error_message = _null_to_nothing(row.error)
        error_message === nothing ||
            push!(failures, ItemFailure(id, "", String(error_message)))
    end
    for row in DBInterface.execute(connection, """
        SELECT item_id, source_item_id, message FROM item_failures
    """)
        push!(failures, ItemFailure(row.source_item_id, row.item_id, row.message))
    end

    metadata = _load_metadata(connection)

    records = ItemRecord[]
    for row in DBInterface.execute(connection, """
        SELECT id, source_item_id, item_label, kind, collection, item_fingerprint FROM items
    """)
        path, timestamp = get(locations, row.source_item_id, (nothing, nothing))
        push!(records, ItemRecord(;
            id=row.id,
            source_item_id=row.source_item_id,
            source_item_fingerprint=get(fingerprints, row.source_item_id, nothing),
            source_item_path=path,
            source_item_timestamp=timestamp,
            item_label=row.item_label,
            kind=Symbol(row.kind),
            collection=Vector{String}(row.collection),
            parameters=get(metadata, (SCOPE_ITEM_PARAMETERS, row.id), MetadataDict()),
            stats=get(metadata, (SCOPE_ITEM_STATS, row.id), MetadataDict()),
            item_fingerprint=_deserialize_hex(row.item_fingerprint),
        ))
    end

    sort!(records; by=record -> (record.source_item_id, record.id))
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

"""
Load the complete browser index from a workspace's open cache, snapshotting committed state.

Reuses the workspace's already-open [`CacheDB`](@ref) instead of opening a second handle to the file.
"""
function load_cache_index(
    cachedb::CacheDB;
    on_progress::Union{Nothing,Function}=nothing,
)::ProjectCacheIndex
    identity = cachedb.identity
    index = @profile_span cachedb.profiler :cache :load_index ProfileAttributes(
        source_id=identity.source_id,
    ) begin
        with_reader(cachedb) do connection
            ProjectCacheIndex(identity, _load_source_scan(connection, identity))
        end
    end
    _emit_cache_load_progress(on_progress, index)
    return index
end

"""Emit the standard cache-load progress event for a freshly read index."""
function _emit_cache_load_progress(
    on_progress::Union{Nothing,Function},
    index::ProjectCacheIndex,
)::Nothing
    emit_progress(
        on_progress;
        phase=:cache_load,
        total_source_items=1,
        processed_source_items=1,
        loaded_items=length(index.source.hierarchy.all_items),
        skipped_source_items=index.source.hierarchy.skipped_count,
        current_source_item=index.identity.cache_path,
    )
    return nothing
end
