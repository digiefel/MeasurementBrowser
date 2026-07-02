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

function _duckdb_sql_type(T::Type)::String
    candidates = Type[m for m in Base.uniontypes(T) if m !== Nothing && m !== Missing]
    isempty(candidates) && error("No DuckDB column type for Julia element type $T")
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
    U <: AbstractArray && return "$(_duckdb_sql_type(eltype(U)))[]"
    error("No DuckDB column type for Julia element type $U")
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

const DataFrameShape = Tuple{Tuple{Vararg{String}},Tuple{Vararg{String}}}

function _dataframe_shape(data::AbstractDataFrame)::DataFrameShape
    return (
        Tuple(String(name) for name in names(data)),
        Tuple(_duckdb_sql_type(eltype(data[!, name])) for name in names(data)),
    )
end

_dataframe_table_name(storage_id::UInt16)::String = "dataframe_$(storage_id)"

function _buffer_rows end

# Disk-backed ordinary tables

mutable struct RowStore{K,R} <: AbstractDiskStore{K,R}
    read_connection::DuckDB.Connection
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

struct TabularBodyBatch
    item_id::String
    body::AbstractDataFrame
end

_buffer_rows(batch::TabularBodyBatch)::Int = nrow(batch.body)

mutable struct TabularFamilyStore <:
               AbstractDiskStore{Tuple{UInt16,UInt32},TabularBodyBatch}
    read_connection::DuckDB.Connection
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
    item_locations::Dict{String,Tuple{UInt16,UInt32}}
    dataframe_schemas::Dict{DataFrameShape,UInt16}
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
    locations = Dict{String,Tuple{UInt16,UInt32}}()
    for row in DBInterface.execute(read_connection, "SELECT item_id, storage_id, seq FROM item_data")
        locations[String(row.item_id)] = (UInt16(row.storage_id), UInt32(row.seq))
    end
    schemas = Dict{DataFrameShape,UInt16}()
    for row in DBInterface.execute(
        read_connection,
        "SELECT storage_id, column_names, column_types FROM dataframe_schemas",
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
    item_id::AbstractString,
    body::AbstractDataFrame,
)::Bool
    id = String(item_id)
    any(name -> name in (MB_SEQ_COLUMN, MB_ROW_COLUMN), names(body)) &&
        error(
            "Tabular item '$id' uses a reserved cache column name " *
            "'$MB_SEQ_COLUMN' or '$MB_ROW_COLUMN'",
        )
    batch = TabularBodyBatch(id, body)
    lock(store.flush_condition)
    try
        _require_writable(store)
        incoming = _buffer_rows(batch)
        old_location = get(store.item_locations, id, nothing)
        old_mutation = old_location === nothing ?
            nothing : get(store.queued, old_location, nothing)
        old_rows = _mutation_rows(store, old_mutation)
        while store.row_limit !== nothing &&
              store.queued_rows - old_rows > 0 &&
              store.queued_rows - old_rows + incoming > store.row_limit
            wait(store.flush_condition)
            _require_writable(store)
            old_location = get(store.item_locations, id, nothing)
            old_mutation = old_location === nothing ?
                nothing : get(store.queued, old_location, nothing)
            old_rows = _mutation_rows(store, old_mutation)
        end
        shape = _dataframe_shape(body)
        store.next_seq < typemax(UInt32) ||
            error("Processed-data cache exhausted its UInt32 sequence space")
        storage_id = get(store.dataframe_schemas, shape, nothing)
        if storage_id === nothing
            next_storage_id = isempty(store.dataframe_schemas) ?
                0 : maximum(Int, values(store.dataframe_schemas)) + 1
            next_storage_id <= typemax(UInt16) ||
                error("Processed-data cache exhausted its UInt16 DataFrame shape space")
            storage_id = UInt16(next_storage_id)
            store.dataframe_schemas[shape] = storage_id
        end
        location = (storage_id, store.next_seq)
        store.next_seq += UInt32(1)
        old_location = pop!(store.item_locations, id, nothing)
        store.item_locations[id] = location
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

function Base.delete!(store::TabularFamilyStore, item_id::AbstractString)::Nothing
    id = String(item_id)
    lock(store.flush_condition) do
        _require_writable(store)
        location = pop!(store.item_locations, id, nothing)
        location === nothing || _queue_delete_location!(store, location)
    end
    return nothing
end

cached_item_ids(store::TabularFamilyStore)::Vector{String} =
    lock(() -> (_require_open(store); collect(keys(store.item_locations))), store.flush_condition)

function _flush_due(store::AbstractDiskStore, now::Float64)::Bool
    isempty(store.queued) && return false
    at_capacity = store.row_limit !== nothing && store.queued_rows >= store.row_limit
    return store.closing || at_capacity ||
        now - store.last_flush >= CACHE_BUFFER_FLUSH_INTERVAL
end

function _flush_loop!(store::AbstractDiskStore{K,R})::Nothing where {K,R}
    while true
        writing = lock(store.flush_condition) do
            while !_flush_due(store, time())
                store.closing && isempty(store.queued) && return nothing
                remaining = max(
                    eps(Float64),
                    CACHE_BUFFER_FLUSH_INTERVAL - (time() - store.last_flush),
                )
                timer = Timer(remaining) do _
                    lock(() -> notify(store.flush_condition), store.flush_condition)
                end
                try
                    wait(store.flush_condition)
                finally
                    close(timer)
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
            _flush_to_db!(store, writing)
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
        for database_row in DBInterface.execute(
            store.read_connection,
            "SELECT * FROM $(store.quoted_table)",
        )
            row = R(database_row)
            values = Tuple(getproperty(database_row, column) for column in store.key_columns)
            raw_key = length(values) == 1 ? only(values) : values
            rows[convert(K, raw_key)] = row
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
)::Tuple{Bool,Union{Nothing,AbstractDataFrame}}
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
)::Dict{Tuple{UInt16,UInt32},AbstractDataFrame}
    by_storage = Dict{UInt16,Vector{UInt32}}()
    for (storage_id, seq) in locations
        push!(get!(() -> UInt32[], by_storage, storage_id), seq)
    end
    rows_by_location = Dict{Tuple{UInt16,UInt32},AbstractDataFrame}()
    for (storage_id, seqs) in by_storage
        unique_seqs = unique(seqs)
        table = _quote_identifier(_dataframe_table_name(storage_id))
        placeholders = join(fill("?", length(unique_seqs)), ", ")
        rows = DataFrame(DBInterface.execute(
            DBInterface.prepare(store.read_connection, """
                SELECT *
                FROM $table
                WHERE $MB_SEQ_COLUMN IN ($placeholders)
                ORDER BY $MB_SEQ_COLUMN, $MB_ROW_COLUMN
            """),
            Tuple(unique_seqs),
        ))
        data_columns = String[
            name for name in names(rows) if name ∉ (MB_SEQ_COLUMN, MB_ROW_COLUMN)
        ]
        for group in groupby(rows, MB_SEQ_COLUMN)
            rows_by_location[(storage_id, UInt32(group[1, MB_SEQ_COLUMN]))] =
                DataFrame(group[:, data_columns])
        end
    end
    return rows_by_location
end

function Base.read(
    store::TabularFamilyStore,
    item_ids::Vector{String},
)::Dict{String,Union{Nothing,AbstractDataFrame}}
    remaining = unique(item_ids)
    results = Dict{String,Union{Nothing,AbstractDataFrame}}()
    while !isempty(remaining)
        locations = Dict{String,Tuple{UInt16,UInt32}}()
        lock(store.flush_condition) do
            _require_open(store)
            for id in remaining
                location = get(store.item_locations, id, nothing)
                if location === nothing
                    results[id] = nothing
                    continue
                end
                handled, data = _pending_tabular_body(store, location)
                if handled
                    results[id] = data
                else
                    locations[id] = location
                end
            end
        end
        isempty(locations) && break
        disk_rows = _read_tabular_locations(
            store,
            collect(Set(values(locations))),
        )
        retry = String[]
        lock(store.flush_condition) do
            _require_open(store)
            for (id, location) in locations
                current_location = get(store.item_locations, id, nothing)
                if current_location === nothing
                    results[id] = nothing
                    continue
                end
                handled, data = _pending_tabular_body(store, current_location)
                if handled
                    results[id] = data
                elseif current_location == location
                    results[id] = disk_rows[location]
                else
                    push!(retry, id)
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
    # DataFrame tables, prunes stale dataframe_schemas rows, checkpoints/compacts
    # DuckDB, and optionally rewrites rows in WorkspaceIndex item order for locality.
    # Keep it out of flush: deletes here should stay blind and cheap.
    isempty(locations) && return nothing
    by_storage = Dict{UInt16,Vector{UInt32}}()
    for (storage_id, seq) in locations
        push!(get!(() -> UInt32[], by_storage, storage_id), seq)
    end
    for (storage_id, seqs) in by_storage
        table = _quote_identifier(_dataframe_table_name(storage_id))
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
                columns, types = _dataframe_shape(first(storage_appends).second.body)
                definitions = join(
                    ("$(_quote_identifier(name)) $type" for
                     (name, type) in zip(columns, types)),
                    ", ",
                )
                table = _quote_identifier(_dataframe_table_name(storage_id))
                DBInterface.execute(
                    store.write_connection,
                    "CREATE TABLE IF NOT EXISTS $table (" *
                    "$MB_SEQ_COLUMN UINTEGER, $MB_ROW_COLUMN UBIGINT, $definitions)",
                )
                DBInterface.execute(
                    DBInterface.prepare(
                        store.write_connection,
                        """
                        INSERT INTO dataframe_schemas VALUES (?, ?, ?)
                        ON CONFLICT (storage_id) DO NOTHING
                        """,
                    ),
                    (storage_id, collect(columns), collect(types)),
                )
                _append_table!(
                    store.write_connection,
                    _dataframe_table_name(storage_id),
                ) do appender
                    for (location, entry) in storage_appends
                        body = entry.body
                        _dataframe_shape(body) == (columns, types) ||
                            error("Storage id $storage_id has incompatible DataFrame shapes")
                        vectors = [body[!, name] for name in columns]
                        for row_index in 1:nrow(body)
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
                    DuckDB.append(appender, entry.item_id)
                    DuckDB.append(appender, location[1])
                    DuckDB.append(appender, location[2])
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
            "SELECT storage_id FROM dataframe_schemas",
        )
            table = _quote_identifier(_dataframe_table_name(UInt16(row.storage_id)))
            DBInterface.execute(store.write_connection, "DROP TABLE $table")
        end
        DBInterface.execute(store.write_connection, "DELETE FROM item_data")
        DBInterface.execute(store.write_connection, "DELETE FROM dataframe_schemas")
    end
    return nothing
end

_reset_memory!(store::RowStore)::Nothing = nothing

function _reset_memory!(store::TabularFamilyStore)::Nothing
    empty!(store.item_locations)
    empty!(store.dataframe_schemas)
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
