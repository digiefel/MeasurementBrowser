using DuckDB
using DBInterface
using SHA
using Serialization
using Dates

const PROJECT_CACHE_SCHEMA_VERSION = 1
const PROJECT_CACHE_LOCK = ReentrantLock()

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

"""The source and DuckDB file belonging to one cache."""
struct ProjectCacheIdentity
    cache_id::String
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

"""All non-payload content needed to restore and compare a project cache."""
struct ProjectCacheIndex
    identity::ProjectCacheIdentity
    source::SourceScan
    analysis_errors::Dict{String,String}
end

# --------------------------------------------------------------------------------------------------
# Identity
# --------------------------------------------------------------------------------------------------

"""
Return the deterministic cache id for a source.

The folder name keeps ids recognizable; the digest distinguishes equal labels at different origins.
"""
function project_cache_id(source)::String
    id = source_id(source)
    name = replace(lowercase(source_label(source)), r"[^a-z0-9_.-]+" => "_")
    isempty(name) && (name = "source")
    return "$name-$(bytes2hex(sha1(id))[1:12])"
end

"""
Bind one cache id to a source.

The returned identity also contains the package-owned DuckDB path.
"""
function project_cache_identity(
    cache_id::AbstractString,
    source,
)::ProjectCacheIdentity
    depot = isempty(DEPOT_PATH) ? homedir() : DEPOT_PATH[1]
    source_dir = replace(lowercase(source_label(source)), r"[^a-z0-9_.-]+" => "_")
    return ProjectCacheIdentity(
        String(cache_id),
        source_id(source),
        source_label(source),
        joinpath(
            depot,
            "measurementbrowser",
            "cache",
            source_dir,
            "$(String(cache_id)).duckdb",
        ),
    )
end

# --------------------------------------------------------------------------------------------------
# Connection + schema
# --------------------------------------------------------------------------------------------------

"""
Open one DuckDB cache connection while holding the process-wide cache lock, then close it.

Serializing all access keeps background cache updates from overlapping plot reads. A persistent
writer/reader connection split is introduced when scanning writes incrementally; for now every
operation opens and closes its own connection.
"""
function with_cache_db(work::Function, path::AbstractString)::Any
    return lock(PROJECT_CACHE_LOCK) do
        connection = DBInterface.connect(DuckDB.DB, String(path))
        try
            return work(connection)
        finally
            DBInterface.close!(connection)
        end
    end
end

"""
One open DuckDB cache file shared by all of a workspace's cache work.

A single [`DuckDB.DB`](@ref) holds the file (and its lock) for the workspace's lifetime. The
persistent `writer` connection serializes every mutation through `writer_lock`; the persistent
`reader` (guarded by `reader_lock`) serves the interactive payload reads without per-call connection
setup. DuckDB's MVCC lets the reader snapshot committed state while a write is in flight, so plot
reads never block the streaming scan. One-shot bulk reads (index load) still use a fresh connection.
"""
mutable struct CacheDB
    identity::ProjectCacheIdentity
    db::DuckDB.DB
    writer::DuckDB.Connection
    writer_lock::ReentrantLock
    reader::DuckDB.Connection
    reader_lock::ReentrantLock
end

"""
Open (creating if needed) the cache file for one workspace and ensure its schema.

A cache that cannot be opened or whose schema cannot be created — a truncated file, a stale
write-ahead log from a hard kill, or an unreadable older layout — is discarded and rebuilt from
scratch rather than bricking the workspace. The data it held is recoverable by rescanning.
"""
function open_cache_db(identity::ProjectCacheIdentity)::CacheDB
    mkpath(dirname(identity.cache_path))
    try
        return _connect_cache_db(identity)
    catch err
        err isa InterruptException && rethrow()
        @warn(
            "Discarding an unreadable project cache and rebuilding it",
            cache_path=identity.cache_path,
            exception=err,
        )
        _remove_cache_files(identity.cache_path)
        return _connect_cache_db(identity)
    end
end

"""Connect to the cache file, ensure its schema, and wrap it in a [`CacheDB`](@ref)."""
function _connect_cache_db(identity::ProjectCacheIdentity)::CacheDB
    db = DBInterface.connect(DuckDB.DB, identity.cache_path)
    try
        ensure_schema!(db)
        writer = DBInterface.connect(db)
        reader = DBInterface.connect(db)
        return CacheDB(identity, db, writer, ReentrantLock(), reader, ReentrantLock())
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

"""Close a workspace cache's connections and the underlying database file."""
function close_cache_db!(cachedb::CacheDB)::Nothing
    DBInterface.close!(cachedb.writer)
    DBInterface.close!(cachedb.reader)
    DBInterface.close!(cachedb.db)
    return nothing
end

"""Run `work` on the persistent reader connection, serialized against other readers."""
function with_persistent_reader(work::Function, cachedb::CacheDB)::Any
    return lock(cachedb.reader_lock) do
        work(cachedb.reader)
    end
end

"""Run `work` on a fresh read connection that snapshots the latest committed cache state."""
function with_reader(work::Function, cachedb::CacheDB)::Any
    connection = DBInterface.connect(cachedb.db)
    try
        return work(connection)
    finally
        DBInterface.close!(connection)
    end
end

"""Run `work` on the persistent writer connection, serialized against other writers."""
function with_writer(work::Function, cachedb::CacheDB)::Any
    return lock(cachedb.writer_lock) do
        work(cachedb.writer)
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
            fingerprint BLOB,
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
            item_fingerprint BLOB)
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
    DBInterface.execute(connection, """
        CREATE TABLE IF NOT EXISTS payloads(
            item_id TEXT PRIMARY KEY,
            source_item_id TEXT,
            sif_hash TEXT,
            if_hash TEXT,
            format TEXT,
            blob BLOB)
    """)
    return nothing
end

# --------------------------------------------------------------------------------------------------
# Encoding helpers
# --------------------------------------------------------------------------------------------------

"""Serialize any value to a byte vector for a BLOB column."""
function _serialize_bytes(value)::Vector{UInt8}
    io = IOBuffer()
    serialize(io, value)
    return take!(io)
end

"""Deserialize a BLOB column value (or `nothing` for SQL NULL)."""
function _deserialize_blob(value)
    value === nothing && return nothing
    ismissing(value) && return nothing
    return deserialize(IOBuffer(Vector{UInt8}(value)))
end

"""Map a SQL NULL (`missing`) to `nothing`, passing other values through."""
_null_to_nothing(value) = ismissing(value) ? nothing : value

"""Content hash of a fingerprint, or `nothing` when the fingerprint is absent."""
_fingerprint_hash(fingerprint) =
    fingerprint === nothing ? nothing : bytes2hex(sha1(_serialize_bytes(fingerprint)))

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
        push!(get!(errors, failure.source_item_id, String[]), failure.message)
    end
    return ProjectCacheIndex(
        identity,
        source,
        Dict(id => join(messages, "\n") for (id, messages) in errors),
    )
end

"""
Compare a cached index with the index produced by the current source scan.

The result counts matching, changed, new, deleted, and failed source items.
"""
function cache_status(
    cached::ProjectCacheIndex,
    source::SourceScan,
)::ProjectCacheStatus
    current = ProjectCacheIndex(cached.identity, source)
    stale = 0
    fresh = 0
    cached_parameters = source_parameter_state(cached.source)
    current_parameters = source_parameter_state(source)
    for (id, fingerprint) in source.source_item_fingerprints
        previous = get(cached.source.source_item_fingerprints, id, nothing)
        previous === nothing && continue
        changed = previous != fingerprint ||
            get(cached_parameters, id, nothing) != get(current_parameters, id, nothing) ||
            get(cached.analysis_errors, id, "") != get(current.analysis_errors, id, "")
        stale += changed
        fresh += !changed && !haskey(current.analysis_errors, id)
    end
    current_ids = keys(source.source_item_fingerprints)
    cached_ids = keys(cached.source.source_item_fingerprints)
    return ProjectCacheStatus(
        length(current_ids),
        length(cached_ids),
        fresh,
        stale,
        count(id -> !haskey(cached.source.source_item_fingerprints, id), current_ids),
        count(id -> !haskey(source.source_item_fingerprints, id), cached_ids),
        length(current.analysis_errors),
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

# A positional INSERT covering every metadata column: scope, entity, key, vtype, then the nine value
# slots in `META_VALUE_COLUMNS` order. Used for the incremental (small-batch) writes; the bulk
# `write_project_cache!` path keeps the faster appender.
const _METADATA_INSERT_SQL = "INSERT INTO metadata VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)"

"""Insert every key of one metadata dict as EAV rows under `scope`/`entity` via a prepared insert."""
function _insert_metadata!(stmt, scope::MetaScope, entity::AbstractString, dict)::Nothing
    for (key, value) in dict
        codec = _meta_codec(value)
        stored = codec.to_db(value)
        DBInterface.execute(stmt, (
            Int8(scope), String(entity), String(key), Int8(codec.vtype),
            (column === codec.column ? stored : nothing for column in META_VALUE_COLUMNS)...,
        ))
    end
    return nothing
end

"""Write the identity/version `meta` rows describing one source scan."""
function _write_meta_rows!(connection, identity::ProjectCacheIdentity, source::SourceScan)::Nothing
    statement = DBInterface.prepare(connection,
        "INSERT INTO meta VALUES (?, ?) ON CONFLICT (key) DO UPDATE SET value = excluded.value")
    for (key, value) in (
        "schema_version" => string(PROJECT_CACHE_SCHEMA_VERSION),
        "cache_id" => identity.cache_id,
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

"""
Create or update a cache from a finished source scan. `replace=true` discards all existing payloads;
a normal update preserves payloads for unchanged source items and removes only those invalidated by a
changed, deleted, or newly failing source item. The whole index (source items, items, metadata) is
rewritten from the scan.
"""
function write_project_cache!(
    identity::ProjectCacheIdentity,
    source::SourceScan;
    replace::Bool=false,
    on_progress::Union{Nothing,Function}=nothing,
)::ProjectCacheIndex
    identity.source_id == source.source_id ||
        error("Cache source does not match source scan")

    new_index = ProjectCacheIndex(identity, source)
    locations = _source_item_locations(source)
    new_hashes = Dict(
        id => _fingerprint_hash(fingerprint)
        for (id, fingerprint) in source.source_item_fingerprints
    )

    mkpath(dirname(identity.cache_path))
    with_cache_db(identity.cache_path) do connection
        ensure_schema!(connection)
        check_cancel()

        # Decide which payloads survive before the index rows are replaced.
        if replace
            DBInterface.execute(connection, "DELETE FROM payloads")
        else
            invalid = String[]
            for row in DBInterface.execute(
                connection, "SELECT source_item_id, fingerprint_hash FROM source_items")
                id = row.source_item_id
                old_hash = _null_to_nothing(row.fingerprint_hash)
                if !haskey(new_hashes, id) || new_hashes[id] != old_hash
                    push!(invalid, id)
                end
            end
            if !isempty(invalid)
                delete_stmt =
                    DBInterface.prepare(connection, "DELETE FROM payloads WHERE source_item_id = ?")
                for id in invalid
                    DBInterface.execute(delete_stmt, (id,))
                end
            end
        end

        DBInterface.execute(connection, "BEGIN TRANSACTION")
        try
            DBInterface.execute(connection, "DELETE FROM source_items")
            DBInterface.execute(connection, "DELETE FROM items")
            DBInterface.execute(connection, "DELETE FROM metadata")
            DBInterface.execute(connection, "DELETE FROM meta")

            _write_meta_rows!(connection, identity, source)

            source_item_stmt = DBInterface.prepare(connection,
                "INSERT INTO source_items VALUES (?, ?, ?, ?, ?, ?)")
            for (id, fingerprint) in source.source_item_fingerprints
                path, timestamp = get(locations, id, (nothing, nothing))
                DBInterface.execute(source_item_stmt, (
                    id,
                    fingerprint === nothing ? nothing : _serialize_bytes(fingerprint),
                    new_hashes[id],
                    path,
                    timestamp,
                    get(new_index.analysis_errors, id, nothing),
                ))
            end

            # `items`/`source_items` carry BLOB fingerprints, which the DuckDB.jl appender mishandles,
            # so they use prepared statements. The far higher-volume `metadata` table is blob-free and
            # written through the (much faster) appender.
            item_stmt = DBInterface.prepare(connection,
                "INSERT INTO items VALUES (?, ?, ?, ?, ?, ?)")
            metadata_appender = DuckDB.Appender(connection, "metadata")
            try
                for record in source.hierarchy.all_items
                    DBInterface.execute(item_stmt, (
                        record.id,
                        record.source_item_id,
                        record.item_label,
                        String(record.kind),
                        record.collection,
                        record.item_fingerprint === nothing ? nothing :
                            _serialize_bytes(record.item_fingerprint),
                    ))
                    _append_metadata!(
                        metadata_appender, SCOPE_ITEM_PARAMETERS, record.id, record.parameters)
                    _append_metadata!(
                        metadata_appender, SCOPE_ITEM_STATS, record.id, record.stats)
                end

                for (path, node) in source.hierarchy.index
                    entity = collection_path_key(collect(String, path))
                    isempty(node.parameters) || _append_metadata!(
                        metadata_appender, SCOPE_NODE_PARAMETERS, entity, node.parameters)
                    isempty(node.stats) || _append_metadata!(
                        metadata_appender, SCOPE_NODE_STATS, entity, node.stats)
                end
                DuckDB.flush(metadata_appender)   # flush buffered rows into the transaction
            finally
                DuckDB.close(metadata_appender)
            end

            DBInterface.execute(connection, "COMMIT")
        catch
            DBInterface.execute(connection, "ROLLBACK")
            rethrow()
        end

        emit_progress(
            on_progress;
            phase=:cache_finalize,
            total_source_items=1,
            processed_source_items=1,
            loaded_items=length(source.hierarchy.all_items),
            skipped_source_items=0,
            current_source_item=identity.cache_path,
        )
    end

    return new_index
end

# --------------------------------------------------------------------------------------------------
# Incremental writing (per source item, as the scan streams)
# --------------------------------------------------------------------------------------------------

"""Delete every index and payload row, leaving an empty (but schema-valid) cache for a rebuild."""
function clear_cache_index!(cachedb::CacheDB)::Nothing
    with_writer(cachedb) do connection
        DBInterface.execute(connection, "BEGIN TRANSACTION")
        try
            for table in ("payloads", "metadata", "items", "source_items", "meta")
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
            "cache_id" => identity.cache_id,
            "source_id" => identity.source_id,
            "source_label" => identity.source_label,
        )
            DBInterface.execute(statement, (key, value))
        end
    end
    return nothing
end

"""Read the current `fingerprint_hash` stored for one source item, or `nothing` if absent."""
function _stored_fingerprint_hash(connection, source_item_id::AbstractString)
    statement = DBInterface.prepare(connection,
        "SELECT fingerprint_hash FROM source_items WHERE source_item_id = ?")
    for row in DBInterface.execute(statement, (source_item_id,))
        return _null_to_nothing(row.fingerprint_hash)
    end
    return missing   # row absent (distinct from a present row whose hash is NULL)
end

"""Delete the index rows (and, when invalidated, payloads) belonging to one source item."""
function _delete_source_item_rows!(
    connection,
    source_item_id::AbstractString;
    drop_payloads::Bool,
)::Nothing
    drop_payloads && DBInterface.execute(
        DBInterface.prepare(connection, "DELETE FROM payloads WHERE source_item_id = ?"),
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
Durably write one source item's records as the scan produces them.

All records belong to the same source item. Existing rows for that source item are replaced in a
single transaction; its payloads are dropped only when the source-item fingerprint changed, so
parameter-only edits keep cached payloads. Returns the source-item id that was written.
"""
function reconcile_source_item!(cachedb::CacheDB, records::Vector{ItemRecord})::String
    isempty(records) && throw(ArgumentError("reconcile_source_item! needs at least one record"))
    source_item_id = first(records).source_item_id
    all(record -> record.source_item_id == source_item_id, records) ||
        throw(ArgumentError("reconcile_source_item! records must share one source item"))
    fingerprint = first(records).source_item_fingerprint
    new_hash = _fingerprint_hash(fingerprint)
    path = first(records).source_item_path
    timestamp = first(records).source_item_timestamp

    with_writer(cachedb) do connection
        DBInterface.execute(connection, "BEGIN TRANSACTION")
        try
            previous_hash = _stored_fingerprint_hash(connection, source_item_id)
            fingerprint_changed = previous_hash === missing || previous_hash != new_hash
            _delete_source_item_rows!(
                connection, source_item_id; drop_payloads=fingerprint_changed)

            DBInterface.execute(
                DBInterface.prepare(connection, "INSERT INTO source_items VALUES (?, ?, ?, ?, ?, ?)"),
                (
                    source_item_id,
                    fingerprint === nothing ? nothing : _serialize_bytes(fingerprint),
                    new_hash,
                    path,
                    timestamp,
                    nothing,
                ))

            item_stmt = DBInterface.prepare(connection, "INSERT INTO items VALUES (?, ?, ?, ?, ?, ?)")
            metadata_stmt = DBInterface.prepare(connection, _METADATA_INSERT_SQL)
            for record in records
                DBInterface.execute(item_stmt, (
                    record.id,
                    record.source_item_id,
                    record.item_label,
                    String(record.kind),
                    record.collection,
                    record.item_fingerprint === nothing ? nothing :
                        _serialize_bytes(record.item_fingerprint),
                ))
                _insert_metadata!(metadata_stmt, SCOPE_ITEM_PARAMETERS, record.id, record.parameters)
                _insert_metadata!(metadata_stmt, SCOPE_ITEM_STATS, record.id, record.stats)
            end

            DBInterface.execute(connection, "COMMIT")
        catch
            DBInterface.execute(connection, "ROLLBACK")
            rethrow()
        end
    end
    return source_item_id
end

"""
Finish an incremental scan: record bookkeeping the per-item writes can't carry.

`written_ids` are the source items already persisted by [`reconcile_source_item!`](@ref). This writes
the `meta` rows and collection-node parameters, records bare rows for skipped/failed source items,
stamps per-source-item analysis errors, and deletes source items (cascading their items, metadata,
and payloads) that the current scan no longer contains.
"""
function finalize_scan!(
    cachedb::CacheDB,
    source::SourceScan,
    written_ids::AbstractSet{<:AbstractString},
)::Nothing
    identity = cachedb.identity
    new_index = ProjectCacheIndex(identity, source)
    locations = _source_item_locations(source)
    current_ids = keys(source.source_item_fingerprints)

    with_writer(cachedb) do connection
        DBInterface.execute(connection, "BEGIN TRANSACTION")
        try
            _write_meta_rows!(connection, identity, source)

            # Bare rows for source items the scan produced no records for (skipped or failed).
            bare_stmt = DBInterface.prepare(connection, "INSERT INTO source_items VALUES (?, ?, ?, ?, ?, ?)")
            for (id, fingerprint) in source.source_item_fingerprints
                id in written_ids && continue
                path, timestamp = get(locations, id, (nothing, nothing))
                DBInterface.execute(bare_stmt, (
                    id,
                    fingerprint === nothing ? nothing : _serialize_bytes(fingerprint),
                    _fingerprint_hash(fingerprint),
                    path,
                    timestamp,
                    get(new_index.analysis_errors, id, nothing),
                ))
            end

            # Stamp errors onto source items that were written but later failed analysis interpretation.
            error_stmt = DBInterface.prepare(connection,
                "UPDATE source_items SET error = ? WHERE source_item_id = ?")
            for (id, message) in new_index.analysis_errors
                id in written_ids || continue
                DBInterface.execute(error_stmt, (message, id))
            end

            # Collection-node parameters (node stats arrive later from analysis).
            node_stmt = DBInterface.prepare(connection, _METADATA_INSERT_SQL)
            for (path, node) in source.hierarchy.index
                entity = collection_path_key(collect(String, path))
                DBInterface.execute(
                    DBInterface.prepare(connection,
                        "DELETE FROM metadata WHERE scope = ? AND entity = ?"),
                    (Int8(SCOPE_NODE_PARAMETERS), entity))
                isempty(node.parameters) ||
                    _insert_metadata!(node_stmt, SCOPE_NODE_PARAMETERS, entity, node.parameters)
            end

            # Drop source items (and everything they own) that the scan no longer contains.
            stored_ids = String[
                row.source_item_id
                for row in DBInterface.execute(connection, "SELECT source_item_id FROM source_items")
            ]
            for id in stored_ids
                id in current_ids && continue
                _delete_source_item_rows!(connection, id; drop_payloads=true)
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
Persist analysis stats as they complete: item stats (scope 1) and collection-node stats (scope 3).

Each entity's stats rows are replaced wholesale so re-analysis overwrites cleanly. Keys are coerced
through [`metadata_dict`](@ref) so analysis output that isn't already a [`MetadataValue`](@ref)
(e.g. an `Int32`) is normalized exactly as the in-memory merge does.
"""
function persist_stats!(
    cachedb::CacheDB,
    item_stats::AbstractDict{String,<:AbstractDict},
    node_stats::AbstractDict{<:Tuple,<:AbstractDict},
)::Nothing
    (isempty(item_stats) && isempty(node_stats)) && return nothing
    with_writer(cachedb) do connection
        DBInterface.execute(connection, "BEGIN TRANSACTION")
        try
            delete_stmt = DBInterface.prepare(connection,
                "DELETE FROM metadata WHERE scope = ? AND entity = ?")
            insert_stmt = DBInterface.prepare(connection, _METADATA_INSERT_SQL)
            for (entity, stats) in item_stats
                DBInterface.execute(delete_stmt, (Int8(SCOPE_ITEM_STATS), entity))
                _insert_metadata!(insert_stmt, SCOPE_ITEM_STATS, entity, metadata_dict(stats))
            end
            for (path, stats) in node_stats
                entity = collection_path_key(collect(String, path))
                DBInterface.execute(delete_stmt, (Int8(SCOPE_NODE_STATS), entity))
                _insert_metadata!(insert_stmt, SCOPE_NODE_STATS, entity, metadata_dict(stats))
            end
            DBInterface.execute(connection, "COMMIT")
        catch
            DBInterface.execute(connection, "ROLLBACK")
            rethrow()
        end
    end
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
    (get(meta, "cache_id", "") == identity.cache_id &&
     get(meta, "source_id", "") == identity.source_id &&
     get(meta, "source_label", "") == identity.source_label) ||
        throw(ProjectCacheError(identity.cache_path, "cache identity does not match"))
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
        fingerprints[id] = _deserialize_blob(row.fingerprint)
        locations[id] = (_null_to_nothing(row.path), _null_to_nothing(row.timestamp))
        error_message = _null_to_nothing(row.error)
        error_message === nothing ||
            push!(failures, ItemFailure(id, "", String(error_message)))
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
            item_fingerprint=_deserialize_blob(row.item_fingerprint),
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
Load the complete browser index without scanning or opening any source file.
"""
function load_project_cache(
    identity::ProjectCacheIdentity;
    on_progress::Union{Nothing,Function}=nothing,
)::ProjectCacheIndex
    isfile(identity.cache_path) ||
        throw(ProjectCacheError(identity.cache_path, "file does not exist"))
    index = with_cache_db(identity.cache_path) do connection
        ProjectCacheIndex(identity, _load_source_scan(connection, identity))
    end
    _emit_cache_load_progress(on_progress, index)
    return index
end

"""
Load the complete browser index from a workspace's open cache, snapshotting committed state.

Unlike [`load_project_cache`](@ref) this reuses the workspace's already-open [`CacheDB`](@ref)
instead of opening a second handle to the same file.
"""
function load_cache_index(
    cachedb::CacheDB;
    on_progress::Union{Nothing,Function}=nothing,
)::ProjectCacheIndex
    identity = cachedb.identity
    index = with_reader(cachedb) do connection
        ProjectCacheIndex(identity, _load_source_scan(connection, identity))
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

# --------------------------------------------------------------------------------------------------
# Payloads
# --------------------------------------------------------------------------------------------------

# A payload is keyed by item id and validated against the current source-item and item fingerprints,
# so a stale blob is never returned. The source-item fingerprint is the primary invalidation token;
# the item fingerprint is optional (recipe items have none) and compared NULL-for-NULL when absent.
# Without a source-item fingerprint there is nothing to validate against, so the item is never cached.
_payload_cacheable(item::ItemRecord)::Bool = item.source_item_fingerprint !== nothing

"""Read valid payloads for `items` on one connection; missing or stale entries stay `nothing`."""
function _read_payloads(connection, items::Vector{ItemRecord})::Vector{Any}
    data = Any[nothing for _ in items]
    statement = DBInterface.prepare(connection,
        "SELECT sif_hash, if_hash, blob FROM payloads WHERE item_id = ?")
    for (position, item) in pairs(items)
        check_cancel()
        _payload_cacheable(item) || continue
        sif = _fingerprint_hash(item.source_item_fingerprint)
        iff = _fingerprint_hash(item.item_fingerprint)   # `nothing` when the item carries no fingerprint
        for row in DBInterface.execute(statement, (item.id,))
            _null_to_nothing(row.sif_hash) == sif || continue
            _null_to_nothing(row.if_hash) == iff || continue
            data[position] = _deserialize_blob(row.blob)
        end
    end
    return data
end

"""Upsert payloads for `items` on one connection; items with absent fingerprints are skipped."""
function _write_payloads!(connection, items::Vector{ItemRecord}, data::Vector)::Nothing
    statement = DBInterface.prepare(connection, """
        INSERT INTO payloads VALUES (?, ?, ?, ?, ?, ?)
        ON CONFLICT (item_id) DO UPDATE SET
            source_item_id = excluded.source_item_id,
            sif_hash = excluded.sif_hash,
            if_hash = excluded.if_hash,
            format = excluded.format,
            blob = excluded.blob
    """)
    for (item, value) in zip(items, data)
        check_cancel()
        (_payload_cacheable(item) && value !== nothing) || continue
        DBInterface.execute(statement, (
            item.id,
            item.source_item_id,
            _fingerprint_hash(item.source_item_fingerprint),
            _fingerprint_hash(item.item_fingerprint),
            "serialize",
            _serialize_bytes(value),
        ))
    end
    return nothing
end

"""
Read every valid requested payload in one DuckDB operation. Missing or stale entries return
`nothing`. A payload is valid only when both the source-item and item fingerprints still match.
"""
function cached_item_data(
    index::Union{Nothing,ProjectCacheIndex},
    items::Vector{ItemRecord},
)::Vector{Any}
    (index === nothing || isempty(items) || !isfile(index.identity.cache_path)) &&
        return Any[nothing for _ in items]
    return with_cache_db(index.identity.cache_path) do connection
        _read_payloads(connection, items)
    end
end

"""
Write valid item payloads into the cache in one DuckDB operation. Items with a `nothing`
source-item or item fingerprint are skipped (never persisted).
"""
function write_item_data_cache!(
    index::Union{Nothing,ProjectCacheIndex},
    items::Vector{ItemRecord},
    data::Vector,
)::Nothing
    length(items) == length(data) ||
        throw(DimensionMismatch("items and data must have equal lengths"))
    (index === nothing || isempty(items) || !isfile(index.identity.cache_path)) && return nothing
    with_cache_db(index.identity.cache_path) do connection
        ensure_schema!(connection)
        _write_payloads!(connection, items, data)
    end
    return nothing
end

"""Read valid cached payloads for `items` from a workspace's open cache (no origin access)."""
function read_item_payloads(cachedb::CacheDB, items::Vector{ItemRecord})::Vector{Any}
    isempty(items) && return Any[]
    return with_persistent_reader(cachedb) do connection
        _read_payloads(connection, items)
    end
end

"""Persist `data` for `items` into a workspace's open cache (skipping uncacheable entries)."""
function write_item_payloads!(
    cachedb::CacheDB,
    items::Vector{ItemRecord},
    data::Vector,
)::Nothing
    length(items) == length(data) ||
        throw(DimensionMismatch("items and data must have equal lengths"))
    isempty(items) && return nothing
    with_writer(cachedb) do connection
        _write_payloads!(connection, items, data)
    end
    return nothing
end
