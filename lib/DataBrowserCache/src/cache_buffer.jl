const CACHE_BUFFER_ROW_LIMIT = 2_000_000
const CACHE_BUFFER_FLUSH_INTERVAL = 2.0
const CACHE_DELETE_KEY_BATCH = 1_000

const MB_SEQ_COLUMN = "__mb_seq"
const MB_ROW_COLUMN = "__mb_row"

@enum BufferMutationKind::UInt8 begin
    BUFFER_APPEND
    BUFFER_EDIT
    BUFFER_DELETE
end

struct BufferMutation{R}
    kind::BufferMutationKind
    row::Union{Nothing,R}
end

abstract type AbstractDiskStore{K,R} end

function _quote_identifier(identifier::AbstractString)::String
    return "\"" * replace(String(identifier), "\"" => "\"\"") * "\""
end

"""The DuckDB column type storing Julia element type `T`, or `nothing` when none does."""
function _duckdb_sql_type_maybe(T::Type)::Union{Nothing,String}
    candidates = Type[m for m in Base.uniontypes(T) if m !== Nothing && m !== Missing]
    length(candidates) == 1 || return nothing
    U = only(candidates)
    U <: Bool && return "BOOLEAN"
    U <: Int8 && return "TINYINT"
    U <: Int16 && return "SMALLINT"
    U <: Int32 && return "INTEGER"
    U <: Int64 && return "BIGINT"
    U <: Int128 && return "HUGEINT"
    U <: UInt8 && return "UTINYINT"
    U <: UInt16 && return "USMALLINT"
    U <: UInt32 && return "UINTEGER"
    U <: UInt64 && return "UBIGINT"
    U <: Float32 && return "FLOAT"
    U <: Float64 && return "DOUBLE"
    U <: AbstractString && return "VARCHAR"
    U <: DateTime && return "TIMESTAMP"
    U <: Date && return "DATE"
    if U <: AbstractArray
        inner = _duckdb_sql_type_maybe(eltype(U))
        return inner === nothing ? nothing : "$(inner)[]"
    end
    return nothing
end

function _duckdb_sql_type(T::Type)::String
    sql = _duckdb_sql_type_maybe(T)
    sql === nothing && error("No DuckDB column type for Julia element type $T")
    return sql
end

function create_table_sql(
    R::Type,
    table::AbstractString,
    key_columns::Tuple{Vararg{Symbol}}=(),
)::String
    definitions = String[
        "$(_quote_identifier(String(name))) $(_duckdb_sql_type(fieldtype(R, name)))"
        for name in fieldnames(R)
    ]
    if !isempty(key_columns)
        key_list = join((_quote_identifier(String(key)) for key in key_columns), ", ")
        push!(definitions, "PRIMARY KEY ($key_list)")
    end
    return "CREATE TABLE IF NOT EXISTS $(_quote_identifier(table)) " *
        "($(join(definitions, ", ")))"
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

function _append_table!(emit::Function, connection, table::AbstractString)::Nothing
    appender = DuckDB.Appender(connection, table)
    try
        emit(appender)
        DuckDB.flush(appender)
    finally
        DuckDB.close(appender)
    end
    return nothing
end

const PayloadShape = Tuple{Tuple{Vararg{String}},Tuple{Vararg{String}}}

"""Column names and DuckDB column types of one tabular payload, via the Tables.jl interface."""
function _payload_shape(data)::PayloadShape
    columns = Tables.columns(data)
    column_names = Tables.columnnames(columns)
    return (
        Tuple(String(name) for name in column_names),
        Tuple(_duckdb_sql_type(eltype(Tables.getcolumn(columns, name)))
              for name in column_names),
    )
end

"""
Whether a payload can land in the columnar payload store: it implements Tables.jl and every
column's element type maps to a DuckDB column type. Anything else stays memory-only.
"""
function _storable_table(payload)::Bool
    Tables.istable(payload) || return false
    columns = Tables.columns(payload)
    column_names = Tables.columnnames(columns)
    isempty(column_names) && return false
    return all(
        _duckdb_sql_type_maybe(eltype(Tables.getcolumn(columns, name))) !== nothing
        for name in column_names)
end

"""Number of rows in one tabular payload."""
_payload_rows(data)::Int = Tables.rowcount(Tables.columns(data))

_payload_table_name(storage_id::UInt16)::String = "payload_$(storage_id)"

function _buffer_rows end

# Disk-backed ordinary tables

mutable struct RowStore{K,R} <: AbstractDiskStore{K,R}
    read_connection::DuckDB.Connection
    # DuckDB connections are not thread-safe; every read_connection use holds this lock.
    read_lock::ReentrantLock
    write_connection::DuckDB.Connection
    table_name::String
    quoted_table::String
    key_columns::Tuple{Vararg{Symbol}}
    quoted_key_columns::Vector{String}
    flush_condition::Base.Threads.Condition
    queued::Dict{K,BufferMutation{R}}
    writing::Dict{K,BufferMutation{R}}
    queued_rows::Int
    writing_rows::Int
    row_limit::Union{Nothing,Int}
    last_flush::Float64
    closing::Bool
    flush_task::Task
end

function RowStore{K,R}(
    database::DuckDB.DB,
    table::AbstractString,
    key_columns::Tuple{Vararg{AbstractString}};
    row_limit::Union{Nothing,Int}=nothing,
) where {K,R}
    row_limit === nothing || row_limit > 0 ||
        throw(ArgumentError("cache buffer row limit must be positive"))
    isempty(key_columns) &&
        throw(ArgumentError("cache buffer requires at least one key column"))
    key_symbols = Tuple(Symbol(String(column)) for column in key_columns)
    quoted_table = _quote_identifier(table)
    quoted_key_columns = String[_quote_identifier(String(column)) for column in key_symbols]
    read_connection = DBInterface.connect(database)
    write_connection = DBInterface.connect(database)
    initial_flush_task = current_task()
    store = RowStore{K,R}(
        read_connection,
        ReentrantLock(),
        write_connection,
        String(table),
        quoted_table,
        key_symbols,
        quoted_key_columns,
        Base.Threads.Condition(),
        Dict{K,BufferMutation{R}}(),
        Dict{K,BufferMutation{R}}(),
        0,
        0,
        row_limit,
        time(),
        false,
        initial_flush_task,
    )
    store.flush_task = Base.Threads.@spawn _flush_loop!(store)
    return store
end

# Disk-backed dynamically-widened metadata table

"""A wide metadata column's logical type: the discriminator restoring the [`MetadataValue`](@ref)."""
@enum MetaVType::Int8 begin
    VT_BOOL = 1
    VT_INT = 2
    VT_FLOAT = 3
    VT_STRING = 4
    VT_SYMBOL = 5
    VT_DATE = 6
    VT_DATETIME = 7
    VT_VBOOL = 9
    VT_VINT = 10
    VT_VFLOAT = 11
    VT_VSTRING = 12
end

"""Classify one [`MetadataValue`](@ref) as its wide-column logical type. `missing` is not stored."""
function _meta_vtype(value)::MetaVType
    value isa Bool && return VT_BOOL
    value isa Int64 && return VT_INT
    value isa Float64 && return VT_FLOAT
    value isa String && return VT_STRING
    value isa Symbol && return VT_SYMBOL
    value isa Date && return VT_DATE
    value isa DateTime && return VT_DATETIME
    value isa Vector{Bool} && return VT_VBOOL
    value isa Vector{Int64} && return VT_VINT
    value isa Vector{Float64} && return VT_VFLOAT
    value isa Vector{String} && return VT_VSTRING
    error("No wide-column vtype for $(typeof(value))")
end

"""The Julia type name a wide-column logical type stands for, for conflict messages."""
function _vtype_julia(vtype::MetaVType)::String
    vtype === VT_BOOL && return "Bool"
    vtype === VT_INT && return "Int64"
    vtype === VT_FLOAT && return "Float64"
    vtype === VT_STRING && return "String"
    vtype === VT_SYMBOL && return "Symbol"
    vtype === VT_DATE && return "Date"
    vtype === VT_DATETIME && return "DateTime"
    vtype === VT_VBOOL && return "Vector{Bool}"
    vtype === VT_VINT && return "Vector{Int64}"
    vtype === VT_VFLOAT && return "Vector{Float64}"
    vtype === VT_VSTRING && return "Vector{String}"
    error("No Julia type name for wide-column vtype $vtype")
end

"""One wide-column type conflict: the name, its registered type, and the rejected value's type."""
const WideConflict = Tuple{Symbol,MetaVType,MetaVType}

"""The DuckDB column type storing one wide-column logical type."""
function _meta_vtype_sql(vtype::MetaVType)::String
    vtype === VT_BOOL && return "BOOLEAN"
    vtype === VT_INT && return "BIGINT"
    vtype === VT_FLOAT && return "DOUBLE"
    vtype === VT_STRING && return "VARCHAR"
    vtype === VT_SYMBOL && return "VARCHAR"
    vtype === VT_DATE && return "TIMESTAMP"
    vtype === VT_DATETIME && return "TIMESTAMP"
    vtype === VT_VBOOL && return "BOOLEAN[]"
    vtype === VT_VINT && return "BIGINT[]"
    vtype === VT_VFLOAT && return "DOUBLE[]"
    vtype === VT_VSTRING && return "VARCHAR[]"
    error("No SQL type for wide-column vtype $vtype")
end

"""The value stored in a wide column for one [`MetadataValue`](@ref) (Symbol/Date normalized to SQL)."""
function _encode_wide_value(value)
    value isa Symbol && return String(value)
    value isa Date && return DateTime(value)
    return value
end

"""Restore one [`MetadataValue`](@ref) from a wide-column cell using the registry's logical type."""
function _decode_wide_value(vtype::MetaVType, raw)
    vtype === VT_BOOL && return Bool(raw)
    vtype === VT_INT && return Int64(raw)
    vtype === VT_FLOAT && return Float64(raw)
    vtype === VT_STRING && return String(raw)
    vtype === VT_SYMBOL && return Symbol(String(raw))
    vtype === VT_DATE && return Date(raw)
    vtype === VT_DATETIME && return DateTime(raw)
    vtype === VT_VBOOL && return Vector{Bool}(raw)
    vtype === VT_VINT && return Vector{Int64}(raw)
    vtype === VT_VFLOAT && return Vector{Float64}(raw)
    vtype === VT_VSTRING && return Vector{String}(raw)
    error("No decoder for wide-column vtype $vtype")
end

"""
A dynamically-widened metadata table: one row per entity, one bare column per metadata name.

Reuses the generic disk-store queue (`append!`/`edit!`/`delete!`/`clear!`/`close!`) with whole-row
replace semantics. `edit!` type-checks each value against the first type registered for its name,
registers unseen names in memory immediately, and returns the keys it dropped for a type conflict;
the domain layer surfaces those as failures. The flush additively `ALTER TABLE ADD COLUMN`s new names
in one transaction with the keyed deletes and the Appender pass.
"""
mutable struct WideRowStore{K} <: AbstractDiskStore{K,MetadataDict}
    read_connection::DuckDB.Connection
    read_lock::ReentrantLock
    write_connection::DuckDB.Connection
    table_name::String
    quoted_table::String
    key_column::Symbol
    quoted_key_column::String
    flush_condition::Base.Threads.Condition
    queued::Dict{K,BufferMutation{MetadataDict}}
    writing::Dict{K,BufferMutation{MetadataDict}}
    queued_rows::Int
    writing_rows::Int
    row_limit::Union{Nothing,Int}
    last_flush::Float64
    closing::Bool
    flush_task::Task
    # Registered columns: their logical type and their physical order (matching the Appender).
    column_vtypes::Dict{Symbol,MetaVType}
    column_order::Vector{Symbol}
    # Names first seen in memory this session, flushed as ALTER TABLE ADD COLUMN + registry insert.
    pending_columns::Vector{Symbol}
end

_buffer_rows(::MetadataDict)::Int = 1

function WideRowStore{K}(
    database::DuckDB.DB,
    table::AbstractString,
    key_column::AbstractString,
) where {K}
    key_symbol = Symbol(String(key_column))
    quoted_table = _quote_identifier(table)
    read_connection = DBInterface.connect(database)
    write_connection = DBInterface.connect(database)
    column_vtypes = Dict{Symbol,MetaVType}()
    for row in DBInterface.execute(
        DBInterface.prepare(
            read_connection,
            "SELECT column_name, vtype FROM wide_columns WHERE table_name = ?"),
        (String(table),),
    )
        column_vtypes[Symbol(String(row.column_name))] = MetaVType(Int8(row.vtype))
    end
    column_order = Symbol[]
    for row in DBInterface.execute(
        DBInterface.prepare(read_connection, """
            SELECT column_name FROM information_schema.columns
            WHERE table_schema = 'main' AND table_name = ?
            ORDER BY ordinal_position
        """),
        (String(table),),
    )
        name = Symbol(String(row.column_name))
        name === key_symbol || push!(column_order, name)
    end
    store = WideRowStore{K}(
        read_connection,
        ReentrantLock(),
        write_connection,
        String(table),
        quoted_table,
        key_symbol,
        _quote_identifier(String(key_column)),
        Base.Threads.Condition(),
        Dict{K,BufferMutation{MetadataDict}}(),
        Dict{K,BufferMutation{MetadataDict}}(),
        0,
        0,
        nothing,
        time(),
        false,
        current_task(),
        column_vtypes,
        column_order,
        Symbol[],
    )
    store.flush_task = Base.Threads.@spawn _flush_loop!(store)
    return store
end

"""
Buffer one whole-row metadata replacement, type-checked against registered columns.

Returns the conflicts dropped from the row: each name whose registered logical type differs from
the incoming value's, with both types. Unseen names register their type immediately, in memory.
The surviving keys are stored; conflicting keys are omitted from the row.
"""
function edit!(store::WideRowStore{K}, key::K, row::MetadataDict)::Vector{WideConflict} where {K}
    dropped = WideConflict[]
    accepted = MetadataDict()
    lock(store.flush_condition) do
        _require_writable(store)
        for (name, value) in row
            value === missing && continue
            incoming = _meta_vtype(value)
            registered = get(store.column_vtypes, name, nothing)
            if registered === nothing
                store.column_vtypes[name] = incoming
                push!(store.column_order, name)
                push!(store.pending_columns, name)
            elseif registered !== incoming
                push!(dropped, (name, registered, incoming))
                continue
            end
            accepted[name] = value
        end
        _accept_rows!(store, key, accepted)
        previous = get(store.queued, key, nothing)
        kind = previous !== nothing && previous.kind === BUFFER_APPEND ?
            BUFFER_APPEND : BUFFER_EDIT
        _set_queued!(store, key, BufferMutation{MetadataDict}(kind, accepted))
    end
    return dropped
end

function Base.read(store::WideRowStore{K})::Dict{K,MetadataDict} where {K}
    while true
        writing, queued, vtypes, disk_columns = lock(store.flush_condition) do
            _require_open(store)
            # Only decode flushed columns from disk; pending names live in queued/writing.
            (
                store.writing,
                store.queued,
                copy(store.column_vtypes),
                setdiff(store.column_order, store.pending_columns),
            )
        end
        rows = Dict{K,MetadataDict}()
        lock(store.read_lock) do
            for database_row in DBInterface.execute(
                store.read_connection, "SELECT * FROM $(store.quoted_table)")
                raw_key = getproperty(database_row, store.key_column)
                dict = MetadataDict()
                for name in disk_columns
                    vtype = vtypes[name]
                    raw = getproperty(database_row, name)
                    ismissing(raw) && continue
                    dict[name] = _decode_wide_value(vtype, raw)
                end
                rows[convert(K, raw_key)] = dict
            end
        end
        stable = lock(store.flush_condition) do
            _require_open(store)
            store.writing === writing && store.queued === queued || return false
            _overlay!(rows, store.writing)
            _overlay!(rows, store.queued)
            true
        end
        stable && return rows
    end
end

function _flush_to_db!(
    store::WideRowStore{K},
    batch::Dict{K,BufferMutation{MetadataDict}},
)::Nothing where {K}
    deleted_keys = K[]
    rows = Pair{K,MetadataDict}[]
    for (key, mutation) in batch
        if mutation.kind === BUFFER_DELETE
            push!(deleted_keys, key)
        else
            mutation.kind === BUFFER_EDIT && push!(deleted_keys, key)
            push!(rows, key => something(mutation.row))
        end
    end
    new_columns, column_order, vtypes = lock(store.flush_condition) do
        (copy(store.pending_columns), copy(store.column_order), copy(store.column_vtypes))
    end
    _transaction(store.write_connection) do
        for name in new_columns
            DBInterface.execute(
                store.write_connection,
                "ALTER TABLE $(store.quoted_table) ADD COLUMN " *
                "$(_quote_identifier(String(name))) $(_meta_vtype_sql(vtypes[name]))",
            )
            DBInterface.execute(
                DBInterface.prepare(
                    store.write_connection, "INSERT INTO wide_columns VALUES (?, ?, ?)"),
                (store.table_name, String(name), Int8(vtypes[name])),
            )
        end
        _delete_keys!(
            store.write_connection,
            store.quoted_table,
            String[store.quoted_key_column],
            deleted_keys,
        )
        isempty(rows) || _append_table!(store.write_connection, store.table_name) do appender
            for (key, dict) in rows
                DuckDB.append(appender, key)
                for name in column_order
                    value = get(dict, name, missing)
                    DuckDB.append(
                        appender, value === missing ? missing : _encode_wide_value(value))
                end
                DuckDB.end_row(appender)
            end
        end
    end
    lock(store.flush_condition) do
        filter!(name -> !(name in new_columns), store.pending_columns)
    end
    return nothing
end

_flush_operation(store::WideRowStore)::Symbol = Symbol("flush_", store.table_name)

_clear_db!(store::WideRowStore)::Nothing =
    (_transaction(() -> DBInterface.execute(
        store.write_connection, "DELETE FROM $(store.quoted_table)"), store.write_connection);
     nothing)

_reset_memory!(::WideRowStore)::Nothing = nothing

# Memory-only FIFO

mutable struct MemoryStore{K,R}
    lock::ReentrantLock
    entries::Dict{K,R}
    order::Vector{K}
    rows::Int
    row_limit::Int
    closed::Bool
end

function MemoryStore{K,R}(; row_limit::Int)::MemoryStore{K,R} where {K,R}
    row_limit > 0 || throw(ArgumentError("memory store row limit must be positive"))
    return MemoryStore{K,R}(ReentrantLock(), Dict{K,R}(), K[], 0, row_limit, false)
end

function Base.append!(store::MemoryStore{K,R}, key::K, row::R)::Bool where {K,R}
    lock(store.lock) do
        store.closed && error("memory store is closed")
        previous = pop!(store.entries, key, nothing)
        if previous !== nothing
            store.rows -= _buffer_rows(previous)
            deleteat!(store.order, something(findfirst(isequal(key), store.order)))
        end
        store.entries[key] = row
        push!(store.order, key)
        store.rows += _buffer_rows(row)
        while store.rows > store.row_limit && length(store.order) > 1
            evicted_key = popfirst!(store.order)
            store.rows -= _buffer_rows(pop!(store.entries, evicted_key))
        end
    end
    return true
end

function Base.delete!(store::MemoryStore{K}, key::K)::Nothing where {K}
    lock(store.lock) do
        store.closed && error("memory store is closed")
        previous = pop!(store.entries, key, nothing)
        previous === nothing && return
        store.rows -= _buffer_rows(previous)
        deleteat!(store.order, something(findfirst(isequal(key), store.order)))
    end
    return nothing
end

function Base.read(store::MemoryStore{K,R}, key::K)::Union{Nothing,R} where {K,R}
    return lock(store.lock) do
        store.closed && error("memory store is closed")
        get(store.entries, key, nothing)
    end
end

function clear!(store::MemoryStore)::Nothing
    lock(store.lock) do
        store.closed && error("memory store is closed")
        empty!(store.entries)
        empty!(store.order)
        store.rows = 0
    end
    return nothing
end

function close!(store::MemoryStore)::Nothing
    lock(store.lock) do
        store.closed && error("memory store is already closed")
        store.closed = true
        empty!(store.entries)
        empty!(store.order)
        store.rows = 0
    end
    return nothing
end

# Disk-backed tabular table family

const PAYLOAD_STAGE_PROCESSED = Int8(0)
const PAYLOAD_STAGE_COLLECTION_PROCESSED = Int8(1)

_payload_stage_code(stage::Symbol)::Int8 =
    stage === :processed ? PAYLOAD_STAGE_PROCESSED :
    stage === :collection_processed ? PAYLOAD_STAGE_COLLECTION_PROCESSED :
    throw(ArgumentError("unknown payload stage '$stage'"))

"""One item's payload key: its integer surrogate and the payload stage it belongs to."""
const PayloadKey = Tuple{Int64,Int8}

# `body` is the caller's original Tables.jl container; pending reads hand it back untouched, and
# the flush loop reads it only through the interface. `container` is the hex-serialized container
# type, carried in the batch because the store's `containers` map can drop a superseded location
# while its append is still being flushed.
struct TabularBodyBatch
    item_key::Int64
    stage::Int8
    body::Any
    container::String
end

_buffer_rows(batch::TabularBodyBatch)::Int = _payload_rows(batch.body)

mutable struct TabularFamilyStore <:
               AbstractDiskStore{Tuple{UInt16,UInt32},TabularBodyBatch}
    read_connection::DuckDB.Connection
    # DuckDB connections are not thread-safe; every read_connection use holds this lock.
    read_lock::ReentrantLock
    write_connection::DuckDB.Connection
    flush_condition::Base.Threads.Condition
    queued::Dict{Tuple{UInt16,UInt32},BufferMutation{TabularBodyBatch}}
    writing::Dict{Tuple{UInt16,UInt32},BufferMutation{TabularBodyBatch}}
    queued_rows::Int
    writing_rows::Int
    row_limit::Union{Nothing,Int}
    last_flush::Float64
    closing::Bool
    flush_task::Task
    item_locations::Dict{PayloadKey,Tuple{UInt16,UInt32}}
    payload_schemas::Dict{PayloadShape,UInt16}
    # Hex-serialized container type per stored payload, rebuilt on read via Tables.materializer.
    containers::Dict{Tuple{UInt16,UInt32},String}
    next_seq::UInt32
end

function TabularFamilyStore(
    database::DuckDB.DB;
    row_limit::Union{Nothing,Int}=nothing,
)::TabularFamilyStore
    row_limit === nothing || row_limit > 0 ||
        throw(ArgumentError("cache buffer row limit must be positive"))
    read_connection = DBInterface.connect(database)
    write_connection = DBInterface.connect(database)
    locations = Dict{PayloadKey,Tuple{UInt16,UInt32}}()
    containers = Dict{Tuple{UInt16,UInt32},String}()
    for row in DBInterface.execute(
        read_connection, "SELECT item_key, stage, storage_id, seq, container FROM item_data")
        location = (UInt16(row.storage_id), UInt32(row.seq))
        locations[(Int64(row.item_key), Int8(row.stage))] = location
        containers[location] = String(row.container)
    end
    schemas = Dict{PayloadShape,UInt16}()
    for row in DBInterface.execute(
        read_connection,
        "SELECT storage_id, column_names, column_types FROM payload_schemas",
    )
        shape = (
            Tuple(String(name) for name in row.column_names),
            Tuple(String(type) for type in row.column_types),
        )
        schemas[shape] = UInt16(row.storage_id)
    end
    max_seq = Int128(only(DBInterface.execute(
        read_connection,
        "SELECT coalesce(max(seq), 0) AS seq FROM item_data",
    )).seq)
    max_seq < typemax(UInt32) ||
        error("Processed-data cache exhausted its UInt32 sequence space")
    store = TabularFamilyStore(
        read_connection,
        ReentrantLock(),
        write_connection,
        Base.Threads.Condition(),
        Dict{Tuple{UInt16,UInt32},BufferMutation{TabularBodyBatch}}(),
        Dict{Tuple{UInt16,UInt32},BufferMutation{TabularBodyBatch}}(),
        0,
        0,
        row_limit,
        time(),
        false,
        current_task(),
        locations,
        schemas,
        containers,
        UInt32(max_seq + 1),
    )
    store.flush_task = Base.Threads.@spawn _flush_loop!(store)
    return store
end

# Shared disk queue

_require_open(store::AbstractDiskStore)::Nothing =
    store.closing ? error("cache store is closed") : nothing

_require_writable(store::AbstractDiskStore)::Nothing =
    (_require_open(store); istaskfailed(store.flush_task) && fetch(store.flush_task); nothing)

function _mutation_rows(
    store::AbstractDiskStore{K,R},
    mutation::Union{Nothing,BufferMutation{R}},
)::Int where {K,R}
    (store.row_limit === nothing || mutation === nothing || mutation.row === nothing) && return 0
    return _buffer_rows(mutation.row)
end

function _set_queued!(
    store::AbstractDiskStore{K,R},
    key::K,
    mutation::Union{Nothing,BufferMutation{R}},
)::Nothing where {K,R}
    previous = get(store.queued, key, nothing)
    previous === nothing || (store.queued_rows -= _mutation_rows(store, previous))
    if mutation === nothing
        delete!(store.queued, key)
    else
        store.queued[key] = mutation
        store.queued_rows += _mutation_rows(store, mutation)
    end
    notify(store.flush_condition)
    return nothing
end

function _accept_rows!(
    store::AbstractDiskStore{K,R},
    key::K,
    row::R,
)::Nothing where {K,R}
    store.row_limit === nothing && return nothing
    incoming = _buffer_rows(row)
    previous = get(store.queued, key, nothing)
    previous_rows = _mutation_rows(store, previous)
    while store.queued_rows - previous_rows > 0 &&
          store.queued_rows - previous_rows + incoming > store.row_limit
        wait(store.flush_condition)
        _require_writable(store)
        previous = get(store.queued, key, nothing)
        previous_rows = _mutation_rows(store, previous)
    end
    return nothing
end

function Base.append!(
    store::AbstractDiskStore{K,R},
    key::K,
    row::R,
)::Bool where {K,R}
    lock(store.flush_condition)
    try
        _require_writable(store)
        _accept_rows!(store, key, row)
        previous = get(store.queued, key, nothing)
        if previous === nothing
            writing = get(store.writing, key, nothing)
            kind = writing === nothing || writing.kind === BUFFER_DELETE ?
                BUFFER_APPEND : BUFFER_EDIT
        else
            kind = previous.kind === BUFFER_APPEND ? BUFFER_APPEND : BUFFER_EDIT
        end
        _set_queued!(store, key, BufferMutation{R}(kind, row))
        return true
    finally
        unlock(store.flush_condition)
    end
end

function edit!(
    store::AbstractDiskStore{K,R},
    key::K,
    row::R,
)::Bool where {K,R}
    lock(store.flush_condition)
    try
        _require_writable(store)
        _accept_rows!(store, key, row)
        previous = get(store.queued, key, nothing)
        kind = previous !== nothing && previous.kind === BUFFER_APPEND ?
            BUFFER_APPEND : BUFFER_EDIT
        _set_queued!(store, key, BufferMutation{R}(kind, row))
        return true
    finally
        unlock(store.flush_condition)
    end
end

function Base.delete!(store::AbstractDiskStore{K,R}, key::K)::Nothing where {K,R}
    lock(store.flush_condition) do
        _require_writable(store)
        previous = get(store.queued, key, nothing)
        if previous !== nothing && previous.kind === BUFFER_APPEND
            _set_queued!(store, key, nothing)
        else
            _set_queued!(store, key, BufferMutation{R}(BUFFER_DELETE, nothing))
        end
    end
    return nothing
end

function _queue_delete_location!(
    store::TabularFamilyStore,
    location::Tuple{UInt16,UInt32},
)::Nothing
    delete!(store.containers, location)
    previous = get(store.queued, location, nothing)
    if previous !== nothing && previous.kind === BUFFER_APPEND
        _set_queued!(store, location, nothing)
    else
        _set_queued!(
            store,
            location,
            BufferMutation{TabularBodyBatch}(BUFFER_DELETE, nothing),
        )
    end
    return nothing
end

function Base.append!(
    store::TabularFamilyStore,
    payload_key::PayloadKey,
    body,
)::Bool
    any(
        name -> String(name) in (MB_SEQ_COLUMN, MB_ROW_COLUMN),
        Tables.columnnames(Tables.columns(body)),
    ) && error(
        "Tabular item '$(payload_key[1])' uses a reserved cache column name " *
        "'$MB_SEQ_COLUMN' or '$MB_ROW_COLUMN'",
    )
    batch = TabularBodyBatch(
        payload_key[1], payload_key[2], body, something(_serialize_hex(typeof(body))))
    lock(store.flush_condition)
    try
        _require_writable(store)
        incoming = _buffer_rows(batch)
        old_location = get(store.item_locations, payload_key, nothing)
        old_mutation = old_location === nothing ?
            nothing : get(store.queued, old_location, nothing)
        old_rows = _mutation_rows(store, old_mutation)
        while store.row_limit !== nothing &&
              store.queued_rows - old_rows > 0 &&
              store.queued_rows - old_rows + incoming > store.row_limit
            wait(store.flush_condition)
            _require_writable(store)
            old_location = get(store.item_locations, payload_key, nothing)
            old_mutation = old_location === nothing ?
                nothing : get(store.queued, old_location, nothing)
            old_rows = _mutation_rows(store, old_mutation)
        end
        shape = _payload_shape(body)
        store.next_seq < typemax(UInt32) ||
            error("Processed-data cache exhausted its UInt32 sequence space")
        storage_id = get(store.payload_schemas, shape, nothing)
        if storage_id === nothing
            next_storage_id = isempty(store.payload_schemas) ?
                0 : maximum(Int, values(store.payload_schemas)) + 1
            next_storage_id <= typemax(UInt16) ||
                error("Processed-data cache exhausted its UInt16 payload shape space")
            storage_id = UInt16(next_storage_id)
            store.payload_schemas[shape] = storage_id
        end
        location = (storage_id, store.next_seq)
        store.next_seq += UInt32(1)
        old_location = pop!(store.item_locations, payload_key, nothing)
        store.item_locations[payload_key] = location
        store.containers[location] = batch.container
        old_location === nothing || _queue_delete_location!(store, old_location)
        _set_queued!(
            store,
            location,
            BufferMutation{TabularBodyBatch}(BUFFER_APPEND, batch),
        )
        return true
    finally
        unlock(store.flush_condition)
    end
end

function Base.delete!(store::TabularFamilyStore, payload_key::PayloadKey)::Nothing
    lock(store.flush_condition) do
        _require_writable(store)
        location = pop!(store.item_locations, payload_key, nothing)
        location === nothing || _queue_delete_location!(store, location)
    end
    return nothing
end

"""Wait on a held condition, waking no later than `deadline`; whether the deadline is still ahead."""
function wait_condition_deadline(condition::Base.Threads.Condition, deadline::Float64)::Bool
    remaining = deadline - time()
    remaining <= 0 && return false
    timer = Timer(remaining) do _
        lock(() -> notify(condition; all=true), condition)
    end
    try
        wait(condition)
    finally
        close(timer)
    end
    return time() < deadline
end

_queue_empty(store::AbstractDiskStore)::Bool = isempty(store.queued)

function _flush_due(store::AbstractDiskStore, now::Float64)::Bool
    _queue_empty(store) && return false
    at_capacity = store.row_limit !== nothing && store.queued_rows >= store.row_limit
    return store.closing || at_capacity ||
        now - store.last_flush >= CACHE_BUFFER_FLUSH_INTERVAL
end

_flush_operation(store::RowStore)::Symbol = Symbol("flush_", store.table_name)
_flush_operation(::TabularFamilyStore)::Symbol = :flush_payload

function _flush_rows(batch::Dict{K,BufferMutation{R}})::Int64 where {K,R}
    rows = 0
    for mutation in values(batch)
        mutation.kind === BUFFER_DELETE && continue
        rows += 1
    end
    return Int64(rows)
end

function _flush_rows(
    batch::Dict{
        Tuple{UInt16,UInt32},
        BufferMutation{TabularBodyBatch},
    },
)::Int64
    rows = 0
    for mutation in values(batch)
        mutation.kind === BUFFER_DELETE && continue
        row = something(mutation.row)
        rows += _payload_rows(row.body)
    end
    return Int64(rows)
end

function _flush_loop!(store::AbstractDiskStore{K,R})::Nothing where {K,R}
    while true
        writing = lock(store.flush_condition) do
            while !_flush_due(store, time())
                store.closing && _queue_empty(store) && return nothing
                if _queue_empty(store)
                    # No queued rows means no flush deadline; sleep until an append notifies.
                    wait(store.flush_condition)
                else
                    wait_condition_deadline(
                        store.flush_condition,
                        store.last_flush + CACHE_BUFFER_FLUSH_INTERVAL,
                    )
                end
            end
            store.writing = store.queued
            store.writing_rows = store.queued_rows
            store.queued = Dict{K,BufferMutation{R}}()
            store.queued_rows = 0
            store.last_flush = time()
            notify(store.flush_condition)
            store.writing
        end
        writing === nothing && return nothing
        try
            @time_dbg _flush_to_db!(store, writing)
        catch error
            lock(store.flush_condition) do
                notify(store.flush_condition, error; all=true, error=true)
            end
            @error "Cache store flush failed" exception=error
            rethrow()
        end
        lock(store.flush_condition) do
            store.writing = Dict{K,BufferMutation{R}}()
            store.writing_rows = 0
            notify(store.flush_condition)
        end
    end
end

function _read_mutation(
    mutations::Dict{K,BufferMutation{R}},
    key::K,
)::Tuple{Bool,Union{Nothing,R}} where {K,R}
    mutation = get(mutations, key, nothing)
    mutation === nothing && return (false, nothing)
    return mutation.kind === BUFFER_DELETE ? (true, nothing) : (true, mutation.row)
end

# Reads

function _overlay!(
    rows::Dict{K,R},
    mutations::Dict{K,BufferMutation{R}},
)::Nothing where {K,R}
    for (key, mutation) in mutations
        if mutation.kind === BUFFER_DELETE
            delete!(rows, key)
        else
            rows[key] = something(mutation.row)
        end
    end
    return nothing
end

function Base.read(store::RowStore{K,R})::Dict{K,R} where {K,R}
    while true
        writing, queued = lock(store.flush_condition) do
            _require_open(store)
            (store.writing, store.queued)
        end
        rows = Dict{K,R}()
        lock(store.read_lock) do
            for database_row in DBInterface.execute(
                store.read_connection,
                "SELECT * FROM $(store.quoted_table)",
            )
                row = R(database_row)
                values = Tuple(
                    getproperty(database_row, column) for column in store.key_columns)
                raw_key = length(values) == 1 ? only(values) : values
                rows[convert(K, raw_key)] = row
            end
        end
        stable = lock(store.flush_condition) do
            _require_open(store)
            store.writing === writing && store.queued === queued || return false
            _overlay!(rows, store.writing)
            _overlay!(rows, store.queued)
            true
        end
        stable && return rows
    end
end

function _pending_tabular_body(
    store::TabularFamilyStore,
    location::Tuple{UInt16,UInt32},
)::Tuple{Bool,Any}
    queued = _read_mutation(store.queued, location)
    queued[1] && return (true, queued[2] === nothing ? nothing : queued[2].body)
    writing = _read_mutation(store.writing, location)
    return writing[1] ?
        (true, writing[2] === nothing ? nothing : writing[2].body) :
        (false, nothing)
end

function _read_tabular_locations(
    store::TabularFamilyStore,
    locations::Vector{Tuple{UInt16,UInt32}},
)::Dict{Tuple{UInt16,UInt32},Any}
    by_storage = Dict{UInt16,Vector{UInt32}}()
    for (storage_id, seq) in locations
        push!(get!(() -> UInt32[], by_storage, storage_id), seq)
    end
    rows_by_location = Dict{Tuple{UInt16,UInt32},Any}()
    for (storage_id, seqs) in by_storage
        unique_seqs = unique(seqs)
        table = _quote_identifier(_payload_table_name(storage_id))
        placeholders = join(fill("?", length(unique_seqs)), ", ")
        rows = lock(store.read_lock) do
            Tables.columntable(DBInterface.execute(
                DBInterface.prepare(store.read_connection, """
                    SELECT *
                    FROM $table
                    WHERE $MB_SEQ_COLUMN IN ($placeholders)
                    ORDER BY $MB_SEQ_COLUMN, $MB_ROW_COLUMN
                """),
                Tuple(unique_seqs),
            ))
        end
        data_names = Tuple(
            name for name in keys(rows) if String(name) ∉ (MB_SEQ_COLUMN, MB_ROW_COLUMN)
        )
        # The query orders by seq, so one payload is one contiguous run of the seq column.
        seq_column = rows[Symbol(MB_SEQ_COLUMN)]
        run_start = 1
        while run_start <= length(seq_column)
            run_end = run_start
            while run_end < length(seq_column) &&
                  seq_column[run_end + 1] == seq_column[run_start]
                run_end += 1
            end
            rows_by_location[(storage_id, UInt32(seq_column[run_start]))] =
                NamedTuple{data_names}(
                    Tuple(rows[name][run_start:run_end] for name in data_names))
            run_start = run_end + 1
        end
        # A zero-row payload has a pointer and schema but no physical rows; it reconstructs as
        # an empty value with the correct columns and types.
        for seq in unique_seqs
            location = (storage_id, seq)
            haskey(rows_by_location, location) && continue
            rows_by_location[location] =
                NamedTuple{data_names}(Tuple(rows[name][1:0] for name in data_names))
        end
    end
    return rows_by_location
end

"""
Rebuild one disk payload as the container type it was stored with, via `Tables.materializer`.
A container type that no longer deserializes (changed project code) is a cache miss — the cache
may always be rebuilt — so the caller recomputes and replaces the entry.
"""
function _materialize_payload(store::TabularFamilyStore, location::Tuple{UInt16,UInt32}, columns)
    hex = get(store.containers, location, nothing)
    hex === nothing && return columns
    container = try
        _deserialize_hex(hex)
    catch
        return nothing
    end
    return Tables.materializer(container)(columns)
end

function Base.read(
    store::TabularFamilyStore,
    payload_keys::Vector{PayloadKey},
)::Dict{PayloadKey,Any}
    remaining = unique(payload_keys)
    results = Dict{PayloadKey,Any}()
    while !isempty(remaining)
        locations = Dict{PayloadKey,Tuple{UInt16,UInt32}}()
        lock(store.flush_condition) do
            _require_open(store)
            for key in remaining
                location = get(store.item_locations, key, nothing)
                if location === nothing
                    results[key] = nothing
                    continue
                end
                handled, data = _pending_tabular_body(store, location)
                if handled
                    results[key] = data
                else
                    locations[key] = location
                end
            end
        end
        isempty(locations) && break
        disk_rows = _read_tabular_locations(
            store,
            collect(Set(values(locations))),
        )
        retry = PayloadKey[]
        lock(store.flush_condition) do
            _require_open(store)
            for (key, location) in locations
                current_location = get(store.item_locations, key, nothing)
                if current_location === nothing
                    results[key] = nothing
                    continue
                end
                handled, data = _pending_tabular_body(store, current_location)
                if handled
                    results[key] = data
                elseif current_location == location
                    results[key] = _materialize_payload(store, location, disk_rows[location])
                else
                    push!(retry, key)
                end
            end
        end
        remaining = retry
    end
    return results
end

# Flush

function _delete_keys_with_temp_table!(
    connection,
    quoted_table::String,
    quoted_key_columns::Vector{String},
    keys::Vector{K},
)::Nothing where {K}
    key_types = Base.fieldtypes(K)
    table_name = "mb_delete_keys_$(time_ns())"
    quoted_temp = _quote_identifier(table_name)
    definitions = join(
        (
            "$(quoted_key_columns[index]) $(_duckdb_sql_type(key_types[index]))"
            for index in eachindex(quoted_key_columns)
        ),
        ", ",
    )
    DBInterface.execute(connection, "CREATE TEMPORARY TABLE $quoted_temp ($definitions)")
    try
        _append_table!(connection, table_name) do appender
            for key in keys
                for value in key
                    DuckDB.append(appender, value)
                end
                DuckDB.end_row(appender)
            end
        end
        predicate = join(
            (
                "$quoted_table.$column = $quoted_temp.$column"
                for column in quoted_key_columns
            ),
            " AND ",
        )
        DBInterface.execute(
            connection,
            "DELETE FROM $quoted_table USING $quoted_temp WHERE $predicate",
        )
    finally
        DBInterface.execute(connection, "DROP TABLE $quoted_temp")
    end
    return nothing
end

function _delete_keys!(
    connection,
    quoted_table::String,
    quoted_key_columns::Vector{String},
    keys::Vector{K},
)::Nothing where {K}
    isempty(keys) && return nothing
    arity = length(quoted_key_columns)
    if arity > 1
        _delete_keys_with_temp_table!(connection, quoted_table, quoted_key_columns, keys)
        return nothing
    end
    for first_index in 1:CACHE_DELETE_KEY_BATCH:length(keys)
        chunk = @view keys[first_index:min(first_index + CACHE_DELETE_KEY_BATCH - 1, end)]
        placeholders = join(fill("?", length(chunk)), ", ")
        sql = "DELETE FROM $quoted_table WHERE $(only(quoted_key_columns)) " *
            "IN ($placeholders)"
        DBInterface.execute(DBInterface.prepare(connection, sql), Tuple(chunk))
    end
    return nothing
end

function _flush_to_db!(
    store::RowStore{K,R},
    batch::Dict{K,BufferMutation{R}},
)::Nothing where {K,R}
    deleted_keys = K[]
    rows = R[]
    for (key, mutation) in batch
        if mutation.kind === BUFFER_DELETE
            push!(deleted_keys, key)
        else
            mutation.kind === BUFFER_EDIT && push!(deleted_keys, key)
            push!(rows, something(mutation.row))
        end
    end
    _transaction(store.write_connection) do
        _delete_keys!(
            store.write_connection,
            store.quoted_table,
            store.quoted_key_columns,
            deleted_keys,
        )
        isempty(rows) || _append_table!(store.write_connection, store.table_name) do appender
            for row in rows
                for column in fieldnames(R)
                    DuckDB.append(appender, getproperty(row, column))
                end
                DuckDB.end_row(appender)
            end
        end
    end
    return nothing
end

function _delete_payloads!(
    connection,
    locations::Vector{Tuple{UInt16,UInt32}},
)::Nothing
    # TODO: add a cold processed-data maintenance pass that drops unreferenced
    # payload tables, prunes stale payload_schemas rows, checkpoints/compacts
    # DuckDB, and optionally rewrites rows in WorkspaceIndex item order for locality.
    # Keep it out of flush: deletes here should stay blind and cheap.
    isempty(locations) && return nothing
    by_storage = Dict{UInt16,Vector{UInt32}}()
    for (storage_id, seq) in locations
        push!(get!(() -> UInt32[], by_storage, storage_id), seq)
    end
    for (storage_id, seqs) in by_storage
        table = _quote_identifier(_payload_table_name(storage_id))
        for first_index in 1:CACHE_DELETE_KEY_BATCH:length(seqs)
            chunk = @view seqs[first_index:min(first_index + CACHE_DELETE_KEY_BATCH - 1, end)]
            placeholders = join(fill("?", length(chunk)), ", ")
            DBInterface.execute(
                DBInterface.prepare(
                    connection,
                    "DELETE FROM $table WHERE $MB_SEQ_COLUMN IN ($placeholders)",
                ),
                Tuple(chunk),
            )
        end
    end
    _delete_keys!(
        connection,
        _quote_identifier("item_data"),
        String[_quote_identifier("storage_id"), _quote_identifier("seq")],
        locations,
    )
    return nothing
end

function _flush_to_db!(
    store::TabularFamilyStore,
    batch::Dict{
        Tuple{UInt16,UInt32},
        BufferMutation{TabularBodyBatch},
    },
)::Nothing
    deleted = Tuple{UInt16,UInt32}[]
    appends = Pair{Tuple{UInt16,UInt32},TabularBodyBatch}[]
    for (location, mutation) in batch
        if mutation.kind === BUFFER_DELETE
            push!(deleted, location)
        elseif mutation.kind === BUFFER_APPEND
            push!(appends, location => something(mutation.row))
        else
            error("TabularFamilyStore replaces physical data with delete plus append, never edit")
        end
    end
    _transaction(store.write_connection) do
        _delete_payloads!(store.write_connection, deleted)
        if !isempty(appends)
            groups =
                Dict{UInt16,Vector{Pair{Tuple{UInt16,UInt32},TabularBodyBatch}}}()
            for append in appends
                push!(get!(
                    () -> Pair{Tuple{UInt16,UInt32},TabularBodyBatch}[],
                    groups,
                    append.first[1],
                ), append)
            end
            for (storage_id, storage_appends) in groups
                columns, types = _payload_shape(first(storage_appends).second.body)
                definitions = join(
                    ("$(_quote_identifier(name)) $type" for
                     (name, type) in zip(columns, types)),
                    ", ",
                )
                table = _quote_identifier(_payload_table_name(storage_id))
                DBInterface.execute(
                    store.write_connection,
                    "CREATE TABLE IF NOT EXISTS $table (" *
                    "$MB_SEQ_COLUMN UINTEGER, $MB_ROW_COLUMN UBIGINT, $definitions)",
                )
                DBInterface.execute(
                    DBInterface.prepare(
                        store.write_connection,
                        """
                        INSERT INTO payload_schemas VALUES (?, ?, ?)
                        ON CONFLICT (storage_id) DO NOTHING
                        """,
                    ),
                    (storage_id, collect(columns), collect(types)),
                )
                _append_table!(
                    store.write_connection,
                    _payload_table_name(storage_id),
                ) do appender
                    for (location, entry) in storage_appends
                        body_columns = Tables.columns(entry.body)
                        _payload_shape(entry.body) == (columns, types) ||
                            error("Storage id $storage_id has incompatible payload shapes")
                        vectors = [
                            Tables.getcolumn(body_columns, Symbol(name)) for name in columns
                        ]
                        for row_index in 1:_payload_rows(entry.body)
                            DuckDB.append(appender, location[2])
                            DuckDB.append(appender, UInt64(row_index))
                            for column in vectors
                                DuckDB.append(appender, column[row_index])
                            end
                            DuckDB.end_row(appender)
                        end
                    end
                end
            end
            _append_table!(store.write_connection, "item_data") do appender
                for (location, entry) in appends
                    DuckDB.append(appender, entry.item_key)
                    DuckDB.append(appender, entry.stage)
                    DuckDB.append(appender, location[1])
                    DuckDB.append(appender, location[2])
                    DuckDB.append(appender, entry.container)
                    DuckDB.end_row(appender)
                end
            end
        end
    end
    return nothing
end

# Lifecycle

function _clear_db!(store::RowStore)::Nothing
    _transaction(store.write_connection) do
        DBInterface.execute(store.write_connection, "DELETE FROM $(store.quoted_table)")
    end
    return nothing
end

function _clear_db!(store::TabularFamilyStore)::Nothing
    _transaction(store.write_connection) do
        for row in DBInterface.execute(
            store.write_connection,
            "SELECT storage_id FROM payload_schemas",
        )
            table = _quote_identifier(_payload_table_name(UInt16(row.storage_id)))
            DBInterface.execute(store.write_connection, "DROP TABLE $table")
        end
        DBInterface.execute(store.write_connection, "DELETE FROM item_data")
        DBInterface.execute(store.write_connection, "DELETE FROM payload_schemas")
    end
    return nothing
end

_reset_memory!(store::RowStore)::Nothing = nothing

function _reset_memory!(store::TabularFamilyStore)::Nothing
    empty!(store.item_locations)
    empty!(store.payload_schemas)
    empty!(store.containers)
    store.next_seq = UInt32(1)
    return nothing
end

function clear!(store::AbstractDiskStore)::Nothing
    lock(store.flush_condition)
    try
        _require_writable(store)
        empty!(store.queued)
        store.queued_rows = 0
        while !isempty(store.writing)
            wait(store.flush_condition)
            _require_writable(store)
        end
        _clear_db!(store)
        _reset_memory!(store)
        notify(store.flush_condition)
    finally
        unlock(store.flush_condition)
    end
    return nothing
end

function close!(store::AbstractDiskStore)::Nothing
    lock(store.flush_condition) do
        store.closing && error("cache store is already closed")
        store.closing = true
        notify(store.flush_condition)
    end
    task_error::Union{Nothing,Exception} = nothing
    try
        wait(store.flush_task)
    catch error
        task_error = error
    end
    try
        DBInterface.close!(store.read_connection)
    finally
        DBInterface.close!(store.write_connection)
    end
    task_error === nothing || throw(task_error)
    return nothing
end
