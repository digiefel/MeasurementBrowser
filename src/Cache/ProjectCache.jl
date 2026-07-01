const PROJECT_CACHE_SCHEMA_VERSION = 8
const MB_SEQ_COLUMN = "__mb_seq"
const MB_ROW_COLUMN = "__mb_row"

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

struct CachedResultState
    key::CacheResultKey
    status::CacheResultStatus
    source_item_id::String
    message::Union{Nothing,String}
end

"""All data-less content needed to restore and compare a project cache."""
struct ProjectCacheIndex
    identity::ProjectCacheIdentity
    source::SourceScan
    analysis_errors::Dict{String,String}
    result_states::Dict{CacheResultKey,CachedResultState}
end

struct SourceItemRow
    id::String
    fingerprint_hex::Union{Nothing,String}
    fingerprint_hash::Union{Nothing,String}
    path::Union{Nothing,String}
    timestamp::Union{Nothing,DateTime}
end

SourceItemRow(row)::SourceItemRow = SourceItemRow(
    String(row.source_item_id),
    _null_to_nothing(row.fingerprint),
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
    _null_to_nothing(row.item_fingerprint),
)

struct MetaRow
    scope::Int8
    entity::String
    key::String
    value::MetadataValue
end

MetaRow(row)::MetaRow = MetaRow(
    Int8(row.scope),
    String(row.entity),
    String(row.key),
    _decode_metadata_value(row.vtype, row),
)

struct FailureRow
    item_id::String
    source_item_id::String
    message::String
end

FailureRow(row)::FailureRow =
    FailureRow(String(row.item_id), String(row.source_item_id), String(row.message))

CachedResultState(row)::CachedResultState = CachedResultState(
    CacheResultKey(CacheResultKind(Int8(row.kind)), String(row.entity)),
    CacheResultStatus(Int8(row.status)),
    String(row.source_item_id),
    _null_to_nothing(row.message),
)

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

A single [`DuckDB.DB`](@ref) holds the file for the workspace's lifetime. Each domain result store is
owned directly as a [`CacheBuffer`](@ref), whose persistent read and write connections allow reads
alongside its coalesced flush transaction. Only explicit schema, meta-header, and checkpoint
operations open short-lived connections. `CacheDB` owns the processed-item locations, DataFrame
shape IDs, and next sequence assigned before payloads enter their buffer.
"""
mutable struct CacheDB
    identity::ProjectCacheIdentity
    db::DuckDB.DB
    metrics::BuildMetrics
    source_items::CacheBuffer{SourceItemRow}
    items::CacheBuffer{ItemRow}
    metadata::CacheBuffer{MetaRow}
    failures::CacheBuffer{FailureRow}
    result_states::CacheBuffer{CachedResultState}
    payload::CacheBuffer{DataItem}
    item_data_index::Dict{String,Tuple{UInt16,UInt32}}
    dataframe_schemas::Dict{
        Tuple{Tuple{Vararg{String}},Tuple{Vararg{String}}},
        UInt16,
    }
    next_item_data_seq::UInt32
    interpreted::CacheBuffer{AbstractDataItem}
    processed_memory::CacheBuffer{AbstractDataItem}
    profiler::Profiling.ProfileSession
end

"""
Open (creating if needed) the cache file for one workspace and ensure its schema.

Open and schema failures propagate with their original cause. Rebuilding generated state is an
explicit caller decision.
"""
function open_cache_db(
    identity::ProjectCacheIdentity,
    profiler::Profiling.ProfileSession,
    metrics::BuildMetrics=BuildMetrics(),
)::CacheDB
    mkpath(dirname(identity.cache_path))
    db = DBInterface.connect(DuckDB.DB, identity.cache_path)
    connection = DBInterface.connect(db)
    try
        try
            DBInterface.execute(
                connection,
                _memory_limit_sql(CACHE_MEMORY_LIMIT_MIB[]),
            )
            ensure_schema!(connection)
            schema_versions = String[
                String(row.value)
                for row in DBInterface.execute(
                    connection,
                    "SELECT value FROM meta WHERE key = 'schema_version'",
                )
            ]
            if !isempty(schema_versions) &&
               only(schema_versions) != string(PROJECT_CACHE_SCHEMA_VERSION)
                throw(ProjectCacheError(identity.cache_path, "cache schema is out of date"))
            end
        finally
            DBInterface.close!(connection)
        end
        payload = CacheBuffer{DataItem}(
            db,
            "item_data",
            ("storage_id", "seq");
            row_limit=CACHE_BUFFER_ROW_LIMIT,
        )
        payload_state = read(payload)
        return CacheDB(
            identity,
            db,
            metrics,
            CacheBuffer{SourceItemRow}(db, "source_items", ("source_item_id",)),
            CacheBuffer{ItemRow}(db, "items", ("id",)),
            CacheBuffer{MetaRow}(db, "metadata", ("scope", "entity", "key")),
            CacheBuffer{FailureRow}(db, "item_failures", ("item_id", "source_item_id")),
            CacheBuffer{CachedResultState}(db, "result_states", ("kind", "entity")),
            payload,
            payload_state.item_data_index,
            payload_state.dataframe_schemas,
            payload_state.next_item_data_seq,
            CacheBuffer{AbstractDataItem}(
                row_limit=CACHE_BUFFER_ROW_LIMIT,
            ),
            CacheBuffer{AbstractDataItem}(
                row_limit=CACHE_BUFFER_ROW_LIMIT,
            ),
            profiler,
        )
    catch
        DBInterface.close!(db)
        rethrow()
    end
end

"""Open a cache with internal profiling disabled for direct cache operations and tests."""
function open_cache_db(identity::ProjectCacheIdentity)::CacheDB
    return open_cache_db(identity, Profiling.ProfileSession(false, false, nothing, nothing))
end

function _buffer_rows(item::AbstractDataItem)::Int
    data = item_data(item)
    return data isa AbstractDataFrame ? nrow(data) : 1
end

function _transaction(work::Function, connection)::Nothing
    DBInterface.execute(connection, "BEGIN TRANSACTION")
    try
        work()
        DBInterface.execute(connection, "COMMIT")
    catch
        DBInterface.execute(connection, "ROLLBACK")
        rethrow()
    end
    return nothing
end

function _batch_rows(
    batch::Dict{Any,BufferMutation{R}},
    kind::BufferMutationKind,
)::Vector{R} where {R}
    return R[
        mutation.row
        for mutation in values(batch)
        if mutation.kind === kind
    ]
end

function _batch_keys(
    batch::Dict,
    kind::BufferMutationKind,
)::Vector{Any}
    return Any[
        key for (key, mutation) in batch if mutation.kind === kind
    ]
end

"""Execute one parameterized multi-row insert/upsert statement."""
function _execute_batched_values!(
    connection,
    statement_head::String,
    rows::Vector{<:Tuple},
    statement_tail::String="",
)::Nothing
    isempty(rows) && return nothing
    values_sql = join(
        fill("(" * join(fill("?", length(first(rows))), ", ") * ")", length(rows)),
        ", ",
    )
    DBInterface.execute(
        DBInterface.prepare(connection, "$statement_head VALUES $values_sql $statement_tail"),
        Tuple(Iterators.flatten(rows)),
    )
    return nothing
end

function _delete_single_keys!(
    connection,
    table::String,
    column::String,
    keys::Vector,
)::Nothing
    isempty(keys) && return nothing
    placeholders = join(fill("?", length(keys)), ", ")
    DBInterface.execute(
        DBInterface.prepare(
            connection,
            "DELETE FROM $table WHERE $column IN ($placeholders)",
        ),
        Tuple(keys),
    )
    return nothing
end

function _append_table!(emit, connection, table::AbstractString)::Nothing
    appender = DuckDB.Appender(connection, table)
    try
        emit(appender)
        DuckDB.flush(appender)
    finally
        DuckDB.close(appender)
    end
    return nothing
end

function _flush_to_db!(
    buffer::CacheBuffer{SourceItemRow},
    batch::Dict{Any,BufferMutation{SourceItemRow}},
)::Nothing
    connection = buffer.write_connection
    _transaction(connection) do
        _delete_single_keys!(
            connection,
            "source_items",
            "source_item_id",
            _batch_keys(batch, BUFFER_DELETE),
        )
        appends = _batch_rows(batch, BUFFER_APPEND)
        if !isempty(appends)
            _append_table!(connection, "source_items") do appender
                for row in appends
                    for value in (
                        row.id,
                        row.fingerprint_hex,
                        row.fingerprint_hash,
                        row.path,
                        row.timestamp,
                    )
                        DuckDB.append(appender, value)
                    end
                    DuckDB.end_row(appender)
                end
            end
        end
        edit_rows = _batch_rows(batch, BUFFER_EDIT)
        if !isempty(edit_rows)
            upserts = Tuple[
                (
                    row.id,
                    row.fingerprint_hex,
                    row.fingerprint_hash,
                    row.path,
                    row.timestamp,
                )
                for row in edit_rows
            ]
            _execute_batched_values!(
                connection,
                "INSERT INTO source_items",
                upserts,
                """
                ON CONFLICT (source_item_id) DO UPDATE SET
                    fingerprint = excluded.fingerprint,
                    fingerprint_hash = excluded.fingerprint_hash,
                    path = excluded.path,
                    timestamp = excluded.timestamp
                """,
            )
        end
    end
    return nothing
end

function _flush_to_db!(
    buffer::CacheBuffer{ItemRow},
    batch::Dict{Any,BufferMutation{ItemRow}},
)::Nothing
    connection = buffer.write_connection
    _transaction(connection) do
        _delete_single_keys!(
            connection,
            "items",
            "id",
            _batch_keys(batch, BUFFER_DELETE),
        )
        appends = _batch_rows(batch, BUFFER_APPEND)
        if !isempty(appends)
            _append_table!(connection, "items") do appender
                for row in appends
                    for value in (
                        row.id,
                        row.source_item_id,
                        row.item_label,
                        row.kind,
                        row.collection,
                        row.item_fingerprint_hex,
                    )
                        DuckDB.append(appender, value)
                    end
                    DuckDB.end_row(appender)
                end
            end
        end
        edit_rows = _batch_rows(batch, BUFFER_EDIT)
        if !isempty(edit_rows)
            upserts = Tuple[
                (
                    row.id,
                    row.source_item_id,
                    row.item_label,
                    row.kind,
                    row.collection,
                    row.item_fingerprint_hex,
                )
                for row in edit_rows
            ]
            _execute_batched_values!(
                connection,
                "INSERT INTO items",
                upserts,
                """
                ON CONFLICT (id) DO UPDATE SET
                    source_item_id = excluded.source_item_id,
                    item_label = excluded.item_label,
                    kind = excluded.kind,
                    collection = excluded.collection,
                    item_fingerprint = excluded.item_fingerprint
                """,
            )
        end
    end
    return nothing
end

"""Encode one typed metadata row in physical column order."""
function _metadata_db_values(row::MetaRow)::Tuple
    codec = _meta_codec(row.value)
    stored = codec.to_db(row.value)
    values = Any[row.scope, row.entity, row.key, Int8(codec.vtype)]
    append!(
        values,
        (column === codec.column ? stored : missing for column in META_VALUE_COLUMNS),
    )
    return Tuple(values)
end

function _flush_to_db!(
    buffer::CacheBuffer{MetaRow},
    batch::Dict{Any,BufferMutation{MetaRow}},
)::Nothing
    connection = buffer.write_connection
    _transaction(connection) do
        keys = _batch_keys(batch, BUFFER_DELETE)
        if !isempty(keys)
            predicates = join(fill("(scope = ? AND entity = ? AND key = ?)", length(keys)), " OR ")
            DBInterface.execute(
                DBInterface.prepare(connection, "DELETE FROM metadata WHERE $predicates"),
                Tuple(Iterators.flatten(keys)),
            )
        end
        appends = _batch_rows(batch, BUFFER_APPEND)
        if !isempty(appends)
            _append_table!(connection, "metadata") do appender
                for row in appends
                    foreach(value -> DuckDB.append(appender, value), _metadata_db_values(row))
                    DuckDB.end_row(appender)
                end
            end
        end
        edit_rows = _batch_rows(batch, BUFFER_EDIT)
        if !isempty(edit_rows)
            upserts = Tuple[
                _metadata_db_values(row)
                for row in edit_rows
            ]
            _execute_batched_values!(
                connection,
                "INSERT INTO metadata",
                upserts,
                """
                ON CONFLICT (scope, entity, key) DO UPDATE SET
                    vtype = excluded.vtype,
                    b = excluded.b,
                    i = excluded.i,
                    d = excluded.d,
                    s = excluded.s,
                    ts = excluded.ts,
                    lb = excluded.lb,
                    li = excluded.li,
                    ld = excluded.ld,
                    ls = excluded.ls
                """,
            )
        end
    end
    return nothing
end

function _flush_to_db!(
    buffer::CacheBuffer{FailureRow},
    batch::Dict{Any,BufferMutation{FailureRow}},
)::Nothing
    connection = buffer.write_connection
    _transaction(connection) do
        keys = _batch_keys(batch, BUFFER_DELETE)
        if !isempty(keys)
            predicates =
                join(fill("(item_id = ? AND source_item_id = ?)", length(keys)), " OR ")
            DBInterface.execute(
                DBInterface.prepare(connection, "DELETE FROM item_failures WHERE $predicates"),
                Tuple(Iterators.flatten(keys)),
            )
        end
        appends = _batch_rows(batch, BUFFER_APPEND)
        if !isempty(appends)
            _append_table!(connection, "item_failures") do appender
                for row in appends
                    for value in (row.item_id, row.source_item_id, row.message)
                        DuckDB.append(appender, value)
                    end
                    DuckDB.end_row(appender)
                end
            end
        end
        edit_rows = _batch_rows(batch, BUFFER_EDIT)
        if !isempty(edit_rows)
            upserts = Tuple[
                (row.item_id, row.source_item_id, row.message)
                for row in edit_rows
            ]
            _execute_batched_values!(
                connection,
                "INSERT INTO item_failures",
                upserts,
                """
                ON CONFLICT (item_id, source_item_id) DO UPDATE SET
                    message = excluded.message
                """,
            )
        end
    end
    return nothing
end

function _flush_to_db!(
    buffer::CacheBuffer{CachedResultState},
    batch::Dict{Any,BufferMutation{CachedResultState}},
)::Nothing
    connection = buffer.write_connection
    _transaction(connection) do
        keys = _batch_keys(batch, BUFFER_DELETE)
        if !isempty(keys)
            predicates = join(fill("(kind = ? AND entity = ?)", length(keys)), " OR ")
            DBInterface.execute(
                DBInterface.prepare(connection, "DELETE FROM result_states WHERE $predicates"),
                Tuple(Iterators.flatten(keys)),
            )
        end
        appends = _batch_rows(batch, BUFFER_APPEND)
        if !isempty(appends)
            _append_table!(connection, "result_states") do appender
                for row in appends
                    for value in (
                        Int8(row.key.kind),
                        row.key.entity,
                        Int8(row.status),
                        row.source_item_id,
                        row.message,
                    )
                        DuckDB.append(appender, value)
                    end
                    DuckDB.end_row(appender)
                end
            end
        end
        edit_rows = _batch_rows(batch, BUFFER_EDIT)
        if !isempty(edit_rows)
            upserts = Tuple[
                (
                    Int8(row.key.kind),
                    row.key.entity,
                    Int8(row.status),
                    row.source_item_id,
                    row.message,
                )
                for row in edit_rows
            ]
            _execute_batched_values!(
                connection,
                "INSERT INTO result_states",
                upserts,
                """
                ON CONFLICT (kind, entity) DO UPDATE SET
                    status = excluded.status,
                    source_item_id = excluded.source_item_id,
                    message = excluded.message
                """,
            )
        end
    end
    return nothing
end

function _quote_identifier(identifier::AbstractString)::String
    return "\"" * replace(String(identifier), "\"" => "\"\"") * "\""
end

function _dataframe_shape(
    data::AbstractDataFrame,
)::Tuple{Tuple{Vararg{String}},Tuple{Vararg{String}}}
    return (
        Tuple(String(name) for name in names(data)),
        Tuple(_duckdb_sql_type(eltype(data[!, name])) for name in names(data)),
    )
end

_dataframe_table_name(storage_id::UInt16)::String = "dataframe_$(storage_id)"

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

function _delete_payloads!(
    connection,
    locations::Vector{Tuple{UInt16,UInt32}},
)::Set{UInt16}
    isempty(locations) && return Set{UInt16}()
    by_storage = Dict{UInt16,Vector{UInt32}}()
    for (storage_id, seq) in locations
        push!(get!(() -> UInt32[], by_storage, storage_id), seq)
    end
    for (storage_id, seqs) in by_storage
        seq_placeholders = join(fill("?", length(seqs)), ", ")
        table = _quote_identifier(_dataframe_table_name(storage_id))
        DBInterface.execute(
            DBInterface.prepare(
                connection,
                "DELETE FROM $table WHERE $MB_SEQ_COLUMN IN ($seq_placeholders)",
            ),
            Tuple(seqs),
        )
    end
    predicates = join(fill("(storage_id = ? AND seq = ?)", length(locations)), " OR ")
    DBInterface.execute(
        DBInterface.prepare(connection, "DELETE FROM item_data WHERE $predicates"),
        Tuple(Iterators.flatten(locations)),
    )
    return Set(keys(by_storage))
end

function _drop_unreferenced_dataframe_schemas!(
    connection,
    storage_ids::Set{UInt16},
)::Nothing
    isempty(storage_ids) && return nothing
    ids = collect(storage_ids)
    placeholders = join(fill("?", length(ids)), ", ")
    retained = Set(
        UInt16(row.storage_id) for row in DBInterface.execute(
            DBInterface.prepare(connection, """
                SELECT DISTINCT storage_id
                FROM item_data
                WHERE storage_id IN ($placeholders)
            """),
            Tuple(ids),
        )
    )
    unused = setdiff(storage_ids, retained)
    isempty(unused) && return nothing
    for storage_id in unused
        DBInterface.execute(
            connection,
            "DROP TABLE $(_quote_identifier(_dataframe_table_name(storage_id)))",
        )
    end
    unused_ids = collect(unused)
    unused_placeholders = join(fill("?", length(unused_ids)), ", ")
    DBInterface.execute(
        DBInterface.prepare(
            connection,
            "DELETE FROM dataframe_schemas WHERE storage_id IN ($unused_placeholders)",
        ),
        Tuple(unused_ids),
    )
    return nothing
end

function _flush_to_db!(
    buffer::CacheBuffer{DataItem},
    batch::Dict{Any,BufferMutation{DataItem}},
)::Nothing
    connection = buffer.write_connection
    deleted = Tuple{UInt16,UInt32}[]
    appends = Tuple{Tuple{UInt16,UInt32},DataItem}[]
    for (key, mutation) in batch
        key isa Tuple{UInt16,UInt32} ||
            error("Payload buffer key must be (UInt16 storage_id, UInt32 seq)")
        if mutation.kind === BUFFER_DELETE
            push!(deleted, key)
        elseif mutation.kind === BUFFER_APPEND
            push!(appends, (key, something(mutation.row)))
        else
            error("Payload buffers replace physical data with delete plus append, never edit")
        end
    end
    _transaction(connection) do
        old_storage_ids = _delete_payloads!(connection, deleted)
        if !isempty(appends)
            existing_storage_ids = Set(
                UInt16(row.storage_id) for row in DBInterface.execute(
                    connection,
                    "SELECT storage_id FROM dataframe_schemas",
                )
            )
            groups = Dict{UInt16,Vector{Tuple{UInt32,DataItem}}}()
            for (location, item) in appends
                storage_id, seq = location
                push!(get!(() -> Tuple{UInt32,DataItem}[], groups, storage_id), (seq, item))
            end
            for (storage_id, payloads) in groups
                first_data = item_data(first(payloads)[2])
                columns, types = _dataframe_shape(first_data)
                if !(storage_id in existing_storage_ids)
                    definitions = join(
                        ("$(_quote_identifier(name)) $type" for
                         (name, type) in zip(columns, types)),
                        ", ",
                    )
                    table = _quote_identifier(_dataframe_table_name(storage_id))
                    DBInterface.execute(
                        connection,
                        "CREATE TABLE $table (" *
                        "$MB_SEQ_COLUMN UINTEGER, $MB_ROW_COLUMN UBIGINT, $definitions)",
                    )
                    DBInterface.execute(
                        DBInterface.prepare(
                            connection,
                            "INSERT INTO dataframe_schemas VALUES (?, ?, ?)",
                        ),
                        (storage_id, collect(columns), collect(types)),
                    )
                end
                table = _dataframe_table_name(storage_id)
                _append_table!(connection, table) do appender
                    for (seq, item) in payloads
                        data = item_data(item)
                        _dataframe_shape(data) == (columns, types) ||
                            error(
                                "Storage id $storage_id was assigned to incompatible " *
                                "DataFrame shapes",
                            )
                        vectors = [data[!, name] for name in columns]
                        for row_index in 1:nrow(data)
                            DuckDB.append(appender, seq)
                            DuckDB.append(appender, UInt64(row_index))
                            foreach(column -> DuckDB.append(appender, column[row_index]), vectors)
                            DuckDB.end_row(appender)
                        end
                    end
                end
            end
            _append_table!(connection, "item_data") do appender
                for ((storage_id, seq), item) in appends
                    for value in (item.id, storage_id, seq)
                        DuckDB.append(appender, value)
                    end
                    DuckDB.end_row(appender)
                end
            end
        end
        _drop_unreferenced_dataframe_schemas!(connection, old_storage_ids)
    end
    return nothing
end

"""Read one processed payload body through queued, writing, then its in-memory disk location."""
function Base.read(
    buffer::CacheBuffer{DataItem},
    location::Tuple{UInt16,UInt32},
)::Any
    pending = lock(buffer.condition) do
        _require_open(buffer)
        queued = _read_mutation(buffer.queued, location)
        queued[1] ? queued : _read_mutation(buffer.writing, location)
    end
    pending[1] && return pending[2] === nothing ? nothing : item_data(pending[2])
    storage_id, seq = location
    table = _quote_identifier(_dataframe_table_name(storage_id))
    data = DataFrame(DBInterface.execute(
        DBInterface.prepare(buffer.read_connection, """
            SELECT * EXCLUDE ($MB_SEQ_COLUMN, $MB_ROW_COLUMN)
            FROM $table
            WHERE $MB_SEQ_COLUMN = ?
            ORDER BY $MB_ROW_COLUMN
        """),
        (seq,),
    ))
    latest = lock(buffer.condition) do
        _require_open(buffer)
        queued = _read_mutation(buffer.queued, location)
        queued[1] ? queued : _read_mutation(buffer.writing, location)
    end
    return latest[1] ?
        (latest[2] === nothing ? nothing : item_data(latest[2])) :
        data
end

"""Load payload routing and allocation state without materializing payload bodies."""
function Base.read(
    buffer::CacheBuffer{DataItem},
)::NamedTuple{
    (:item_data_index, :dataframe_schemas, :next_item_data_seq),
    Tuple{
        Dict{String,Tuple{UInt16,UInt32}},
        Dict{Tuple{Tuple{Vararg{String}},Tuple{Vararg{String}}},UInt16},
        UInt32,
    },
}
    locations = Dict{String,Tuple{UInt16,UInt32}}()
    for row in DBInterface.execute(
        buffer.read_connection,
        "SELECT item_id, storage_id, seq FROM item_data",
    )
        locations[String(row.item_id)] = (UInt16(row.storage_id), UInt32(row.seq))
    end
    dataframe_schemas =
        Dict{Tuple{Tuple{Vararg{String}},Tuple{Vararg{String}}},UInt16}()
    for row in DBInterface.execute(
        buffer.read_connection,
        "SELECT storage_id, column_names, column_types FROM dataframe_schemas",
    )
        shape = (
            Tuple(String(name) for name in row.column_names),
            Tuple(String(type) for type in row.column_types),
        )
        dataframe_schemas[shape] = UInt16(row.storage_id)
    end
    max_seq = Int128(only(DBInterface.execute(
        buffer.read_connection,
        "SELECT coalesce(max(seq), 0) AS seq FROM item_data",
    )).seq)
    max_seq < typemax(UInt32) ||
        error("Processed-data cache exhausted its UInt32 sequence space")
    lock(buffer.condition) do
        _require_open(buffer)
        isempty(buffer.writing) && isempty(buffer.queued) ||
            error("Complete payload state may only be loaded before producers start")
    end
    return (
        item_data_index=locations,
        dataframe_schemas=dataframe_schemas,
        next_item_data_seq=UInt32(max_seq + 1),
    )
end

"""Clear the complete processed-payload table family in one transaction."""
function _clear_db!(
    buffer::CacheBuffer{DataItem},
)::Nothing
    connection = buffer.write_connection
    _transaction(connection) do
        storage_ids = UInt16[
            UInt16(row.storage_id) for
            row in DBInterface.execute(connection, "SELECT storage_id FROM dataframe_schemas")
        ]
        for storage_id in storage_ids
            DBInterface.execute(
                connection,
                "DROP TABLE $(_quote_identifier(_dataframe_table_name(storage_id)))",
            )
        end
        DBInterface.execute(connection, "DELETE FROM item_data")
        DBInterface.execute(connection, "DELETE FROM dataframe_schemas")
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

Item statistics are not written here; processing publishes them independently.
"""
function store_interpreted!(
    cache::CacheDB,
    records::Vector{ItemRecord},
    data::Vector{<:AbstractDataItem},
)::Nothing
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
                MetaRow(Int8(SCOPE_ITEM_PARAMETERS), record.id, String(key), value))
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

"""
Store one processed result in the disk-backed payload buffer when cacheable, or in the memory-only
processed buffer otherwise. Item statistics are stored independently with [`store_item_stats!`](@ref).

A non-cacheable (or non-tabular) processed item is therefore held exactly like an interpreted one —
bounded, never written, recomputed from source on a later miss — so a re-selection still skips
reprocessing while the item is resident (docs/cache.md "Two backing modes").
"""
function store_processed!(
    cache::CacheDB,
    record::ItemRecord,
    item::AbstractDataItem,
)::Nothing
    payload = item_data(item)
    disk = cacheable(item) && payload isa AbstractDataFrame && !isempty(names(payload))
    if disk
        started = time_ns()
        old_location = get(cache.item_data_index, record.id, nothing)
        any(name -> name in (MB_SEQ_COLUMN, MB_ROW_COLUMN), names(payload)) &&
            error(
                "Processed item '$(record.id)' uses a reserved cache column name " *
                "'$MB_SEQ_COLUMN' or '$MB_ROW_COLUMN'",
            )
        shape = _dataframe_shape(payload)
        storage_id = get(cache.dataframe_schemas, shape, nothing)
        if storage_id === nothing
            next_storage_id = isempty(cache.dataframe_schemas) ?
                0 : maximum(Int, values(cache.dataframe_schemas)) + 1
            next_storage_id <= typemax(UInt16) ||
                error("Processed-data cache exhausted its UInt16 DataFrame shape space")
            storage_id = UInt16(next_storage_id)
            cache.dataframe_schemas[shape] = storage_id
        end
        cache.next_item_data_seq < typemax(UInt32) ||
            error("Processed-data cache exhausted its UInt32 sequence space")
        seq = cache.next_item_data_seq
        cache.next_item_data_seq += UInt32(1)
        new_location = (storage_id, seq)
        old_location === nothing || delete!(cache.payload, old_location)
        cache.item_data_index[record.id] = new_location
        delete!(cache.processed_memory, record.id)
        append!(cache.payload, new_location, DataItem(record, payload))
        record_cache_phase!(
            cache.metrics.processed_write_ns,
            cache.metrics.processed_writes,
            started,
        )
    else
        old_location = pop!(cache.item_data_index, record.id, nothing)
        old_location === nothing || delete!(cache.payload, old_location)
        append!(cache.processed_memory, record.id, item)
    end
    delete!(cache.result_states, (Int8(PROCESSING_RESULT), record.id))
    return nothing
end

"""Persist one independent failed work result."""
function store_result_failure!(
    cache::CacheDB,
    key::CacheResultKey,
    source_item_id::AbstractString,
    message::AbstractString,
)::Nothing
    edit!(
        cache.result_states,
        (Int8(key.kind), key.entity),
        CachedResultState(key, RESULT_FAILED, String(source_item_id), String(message)),
    )
    return nothing
end

"""Persist one item-stat result, including a successful empty result."""
function store_item_stats!(
    cache::CacheDB,
    record::ItemRecord,
    stats::AbstractDict,
)::Nothing
    for (key, value) in metadata_dict(stats)
        edit!(
            cache.metadata,
            (Int8(SCOPE_ITEM_STATS), record.id, String(key)),
            MetaRow(Int8(SCOPE_ITEM_STATS), record.id, String(key), value),
        )
    end
    result_key = CacheResultKey(ITEM_STATS_RESULT, record.id)
    edit!(
        cache.result_states,
        (Int8(result_key.kind), result_key.entity),
        CachedResultState(result_key, RESULT_READY, record.source_item_id, nothing),
    )
    return nothing
end

"""Persist one collection-stat result, including a successful empty result."""
function store_collection_stats!(
    cache::CacheDB,
    collection_key::AbstractString,
    stats::AbstractDict,
)::Nothing
    entity = String(collection_key)
    for (key, value) in metadata_dict(stats)
        edit!(
            cache.metadata,
            (Int8(SCOPE_NODE_STATS), entity, String(key)),
            MetaRow(Int8(SCOPE_NODE_STATS), entity, String(key), value),
        )
    end
    result_key = CacheResultKey(COLLECTION_STATS_RESULT, entity)
    edit!(
        cache.result_states,
        (Int8(result_key.kind), result_key.entity),
        CachedResultState(result_key, RESULT_READY, "", nothing),
    )
    return nothing
end

"""Delete one collection's persisted statistics and completion state."""
function delete_collection_stats!(
    cache::CacheDB,
    collection_key::AbstractString,
)::Nothing
    entity = String(collection_key)
    for key in keys(read(cache.metadata))
        key[1] == Int8(SCOPE_NODE_STATS) && key[2] == entity ||
            continue
        delete!(cache.metadata, key)
    end
    delete!(
        cache.result_states,
        (Int8(COLLECTION_STATS_RESULT), entity),
    )
    return nothing
end

"""Delete every cached descendant of one source item through typed buffer mutations."""
function delete_source_item!(
    cache::CacheDB,
    source_item_id::AbstractString,
)::Nothing
    source_id = String(source_item_id)
    owned_items = String[
        row.id for row in values(read(cache.items)) if row.source_item_id == source_id
    ]
    owned = Set(owned_items)
    delete!(cache.source_items, source_id)
    for item_id in owned_items
        old_location = pop!(cache.item_data_index, item_id, nothing)
        delete!(cache.items, item_id)
        old_location === nothing || delete!(cache.payload, old_location)
        delete!(cache.interpreted, item_id)
        delete!(cache.processed_memory, item_id)
    end
    for (key, row) in read(cache.metadata)
        row.scope in (Int8(SCOPE_ITEM_PARAMETERS), Int8(SCOPE_ITEM_STATS)) &&
            row.entity in owned || continue
        delete!(cache.metadata, key)
    end
    for (key, row) in read(cache.failures)
        row.source_item_id == source_id || continue
        delete!(cache.failures, key)
    end
    for (key, row) in read(cache.result_states)
        row.source_item_id == source_id || continue
        delete!(cache.result_states, key)
    end
    return nothing
end

"""Aggregate staged durable keys and payload rows."""
function cache_pending_counts(
    cache::CacheDB,
)::NamedTuple{(:items, :rows),Tuple{Int,Int}}
    items = 0
    rows = 0
    for buffer in (
        cache.source_items,
        cache.items,
        cache.metadata,
        cache.failures,
        cache.result_states,
        cache.payload,
    )
        lock(buffer.condition) do
            items += length(buffer.queued) + length(buffer.writing)
            rows += buffer.queued_rows + buffer.writing_rows
        end
    end
    return (items=items, rows=rows)
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

    if stage === :interpreted
        return Any[read(cache.interpreted, record.id) for record in records]
    end

    loaded = Any[]
    for record in records
        held = read(cache.processed_memory, record.id)
        if held !== nothing
            push!(loaded, held)
            continue
        end
        location = get(cache.item_data_index, record.id, nothing)
        data = location === nothing ? nothing : read(cache.payload, location)
        push!(loaded, data === nothing ? nothing : DataItem(record, data))
    end
    return loaded
end

"""Close a workspace cache's database file, checkpointing first."""
function close_cache_db!(cachedb::CacheDB)::Nothing
    failure::Union{Nothing,Exception} = nothing
    for buffer in (
        cachedb.source_items,
        cachedb.items,
        cachedb.metadata,
        cachedb.failures,
        cachedb.result_states,
        cachedb.payload,
        cachedb.interpreted,
        cachedb.processed_memory,
    )
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
            timestamp TIMESTAMP)
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
    DBInterface.execute(connection, """
        CREATE TABLE IF NOT EXISTS item_failures(
            item_id TEXT,
            source_item_id TEXT,
            message TEXT,
            PRIMARY KEY (item_id, source_item_id))
    """)
    DBInterface.execute(connection, """
        CREATE TABLE IF NOT EXISTS result_states(
            kind TINYINT,
            entity TEXT,
            status TINYINT,
            source_item_id TEXT,
            message TEXT,
            PRIMARY KEY (kind, entity))
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

"""Hex-text storage and content hash for one fingerprint, as `(hex, hash)` (`nothing` when absent)."""
function _encode_fingerprint(fingerprint)::Tuple{Union{Nothing,String},Union{Nothing,String}}
    fingerprint === nothing && return (nothing, nothing)
    bytes = _serialize_bytes(fingerprint)
    return (bytes2hex(bytes), bytes2hex(sha1(bytes)))
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
        Dict{CacheResultKey,CachedResultState}(),
    )
end

# --------------------------------------------------------------------------------------------------
# Writing the index
# --------------------------------------------------------------------------------------------------

# --------------------------------------------------------------------------------------------------
# Incremental writing (per source item, as the scan streams)
# --------------------------------------------------------------------------------------------------

"""Delete every index and item-data row, leaving an empty (but schema-valid) cache for a rebuild."""
function clear_cache_index!(cachedb::CacheDB)::Nothing
    for buffer in (
        cachedb.source_items,
        cachedb.items,
        cachedb.metadata,
        cachedb.failures,
        cachedb.result_states,
        cachedb.payload,
        cachedb.interpreted,
        cachedb.processed_memory,
    )
        clear!(buffer)
    end
    empty!(cachedb.item_data_index)
    empty!(cachedb.dataframe_schemas)
    cachedb.next_item_data_seq = UInt32(1)
    connection = DBInterface.connect(cachedb.db)
    try
        DBInterface.execute(connection, "DELETE FROM meta")
    finally
        DBInterface.close!(connection)
    end
    return nothing
end

"""
Stamp the identity `meta` rows so a scan's incremental writes are loadable before it finishes.

Without these rows [`load_cache_index`](@ref) treats the cache as unbuilt, so an interrupted scan
would discard the per-item progress already written. The scan-wide `skipped_count` and
`has_collection_parameters` are derived on load, not stored.
"""
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

# --------------------------------------------------------------------------------------------------
# Loading the index
# --------------------------------------------------------------------------------------------------

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
        dict[Symbol(row.key)] = row.value
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
        base = ProjectCacheIndex(identity, _load_source_scan(cachedb))
        records = Dict(record.id => record for record in base.source.hierarchy.all_items)
        result_states = Dict{CacheResultKey,CachedResultState}()
        for item_id in keys(cachedb.item_data_index)
            record = get(records, item_id, nothing)
            record === nothing && continue
            key = CacheResultKey(PROCESSING_RESULT, item_id)
            result_states[key] =
                CachedResultState(key, RESULT_READY, record.source_item_id, nothing)
        end
        for state in values(read(cachedb.result_states))
            result_states[state.key] = state
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
