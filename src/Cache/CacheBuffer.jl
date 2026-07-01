const CACHE_BUFFER_ROW_LIMIT = 2_000_000
const CACHE_BUFFER_FLUSH_INTERVAL = 2.0

@enum BufferMutationKind::UInt8 begin
    BUFFER_APPEND
    BUFFER_EDIT
    BUFFER_DELETE
end

"""One final coalesced mutation for a key."""
struct BufferMutation{R}
    kind::BufferMutationKind
    row::Union{Nothing,R}
end

"""
Live bounded memory in front of one result store.

Disk-backed buffers own separate persistent read and write connections. Ordinary stores also own
their quoted table name, key columns, and prepared keyed read.
"""
mutable struct CacheBuffer{R,S}
    read_connection::Union{Nothing,DuckDB.Connection}
    write_connection::Union{Nothing,DuckDB.Connection}
    table::Union{Nothing,String}
    key_columns::Tuple{Vararg{Symbol}}
    keyed_read::Union{Nothing,S}
    condition::Base.Threads.Condition
    queued::Dict{Any,BufferMutation{R}}
    writing::Dict{Any,BufferMutation{R}}
    queued_rows::Int
    writing_rows::Int
    row_limit::Union{Nothing,Int}
    last_flush::Float64
    closing::Bool
    flush_task::Union{Nothing,Task}
end

"""Validate and quote one configured SQL identifier."""
function _buffer_identifier(identifier::AbstractString)::String
    value = String(identifier)
    occursin(r"^[A-Za-z_][A-Za-z0-9_]*$", value) ||
        throw(ArgumentError("invalid cache buffer SQL identifier '$value'"))
    return "\"$value\""
end

"""
    CacheBuffer{R}(database, table, key_columns; row_limit=nothing)

Construct a live disk-backed buffer with generic ordinary-table read configuration.
"""
function CacheBuffer{R}(
    database::DuckDB.DB,
    table::AbstractString,
    key_columns::Tuple{Vararg{AbstractString}};
    row_limit::Union{Nothing,Int}=nothing,
) where {R}
    row_limit === nothing || row_limit > 0 ||
        throw(ArgumentError("cache buffer row limit must be positive"))
    isempty(key_columns) &&
        throw(ArgumentError("cache buffer requires at least one key column"))
    quoted_table = _buffer_identifier(table)
    column_symbols = Tuple(Symbol(String(column)) for column in key_columns)
    quoted_columns = String[_buffer_identifier(column) for column in key_columns]
    read_connection = DBInterface.connect(database)
    write_connection = try
        DBInterface.connect(database)
    catch
        DBInterface.close!(read_connection)
        rethrow()
    end
    keyed_read = try
        predicate = join(("$column = ?" for column in quoted_columns), " AND ")
        DBInterface.prepare(
            read_connection,
            "SELECT * FROM $quoted_table WHERE $predicate",
        )
    catch
        DBInterface.close!(read_connection)
        DBInterface.close!(write_connection)
        rethrow()
    end
    buffer = CacheBuffer{R,typeof(keyed_read)}(
        read_connection,
        write_connection,
        quoted_table,
        column_symbols,
        keyed_read,
        Base.Threads.Condition(),
        Dict{Any,BufferMutation{R}}(),
        Dict{Any,BufferMutation{R}}(),
        0,
        0,
        row_limit,
        time(),
        false,
        nothing,
    )
    buffer.flush_task = Base.Threads.@spawn _flush_loop!(buffer)
    return buffer
end

"""Construct a live memory-only buffer."""
function CacheBuffer{R}(; row_limit::Union{Nothing,Int}=nothing) where {R}
    row_limit === nothing || row_limit > 0 ||
        throw(ArgumentError("cache buffer row limit must be positive"))
    return CacheBuffer{R,Nothing}(
        nothing,
        nothing,
        nothing,
        (),
        nothing,
        Base.Threads.Condition(),
        Dict{Any,BufferMutation{R}}(),
        Dict{Any,BufferMutation{R}}(),
        0,
        0,
        row_limit,
        time(),
        false,
        nothing,
    )
end

"""Fail an operation after the buffer entered its permanent closing state. Caller holds `condition`."""
function _require_open(buffer::CacheBuffer)::Nothing
    buffer.closing && error("CacheBuffer is closed and cannot accept further operations")
    return nothing
end

"""Reject mutation or clearing after a terminal flush failure. Caller holds `condition`."""
function _require_writable(buffer::CacheBuffer)::Nothing
    _require_open(buffer)
    buffer.flush_task !== nothing && istaskfailed(buffer.flush_task) &&
        fetch(buffer.flush_task)
    return nothing
end

"""Rows retained by one payload. Defined only for payload types configured with a row limit."""
function _buffer_rows end

"""Return one mutation's row count without consulting payload methods for unlimited stores."""
function _mutation_rows(buffer::CacheBuffer, mutation)::Int
    (buffer.row_limit === nothing || mutation === nothing || mutation.row === nothing) && return 0
    return _buffer_rows(mutation.row)
end

"""Replace one queued mutation and update queued row capacity. Caller holds `condition`."""
function _set_queued!(
    buffer::CacheBuffer{R},
    key,
    mutation::Union{Nothing,BufferMutation{R}},
)::Nothing where {R}
    previous = get(buffer.queued, key, nothing)
    previous === nothing || (buffer.queued_rows -= _mutation_rows(buffer, previous))
    if mutation === nothing
        delete!(buffer.queued, key)
    else
        buffer.queued[key] = mutation
        buffer.queued_rows += _mutation_rows(buffer, mutation)
    end
    notify(buffer.condition)
    return nothing
end

"""Wait for queued disk capacity or reject a memory-only payload that would cross its row limit."""
function _accept_rows!(buffer::CacheBuffer, key, row)::Bool
    buffer.row_limit === nothing && return true
    incoming = _buffer_rows(row)
    previous = get(buffer.queued, key, nothing)
    previous_rows = _mutation_rows(buffer, previous)
    if buffer.write_connection === nothing
        return buffer.queued_rows - previous_rows + incoming <= buffer.row_limit
    end
    while length(buffer.queued) > (previous === nothing ? 0 : 1) &&
          buffer.queued_rows - previous_rows + incoming > buffer.row_limit
        wait(buffer.condition)
        _require_writable(buffer)
        previous = get(buffer.queued, key, nothing)
        previous_rows = _mutation_rows(buffer, previous)
    end
    return true
end

"""Record append intent in the active queue."""
function Base.append!(buffer::CacheBuffer{R}, key, row::R)::Bool where {R}
    lock(buffer.condition)
    try
        _require_writable(buffer)
        _accept_rows!(buffer, key, row) || return false
        previous = get(buffer.queued, key, nothing)
        kind = previous === nothing || previous.kind === BUFFER_APPEND ?
            BUFFER_APPEND : BUFFER_EDIT
        _set_queued!(buffer, key, BufferMutation{R}(kind, row))
        return true
    finally
        unlock(buffer.condition)
    end
end

"""Record edit intent in the active queue. A queued append remains an append."""
function edit!(buffer::CacheBuffer{R}, key, row::R)::Bool where {R}
    lock(buffer.condition)
    try
        _require_writable(buffer)
        _accept_rows!(buffer, key, row) || return false
        previous = get(buffer.queued, key, nothing)
        kind = previous !== nothing && previous.kind === BUFFER_APPEND ?
            BUFFER_APPEND : BUFFER_EDIT
        _set_queued!(buffer, key, BufferMutation{R}(kind, row))
        return true
    finally
        unlock(buffer.condition)
    end
end

"""Record delete intent in the active queue. A queued append followed by delete cancels."""
function Base.delete!(buffer::CacheBuffer{R}, key)::Nothing where {R}
    lock(buffer.condition) do
        _require_writable(buffer)
        previous = get(buffer.queued, key, nothing)
        if buffer.write_connection === nothing ||
           (previous !== nothing && previous.kind === BUFFER_APPEND)
            _set_queued!(buffer, key, nothing)
        else
            _set_queued!(buffer, key, BufferMutation{R}(BUFFER_DELETE, nothing))
        end
    end
    return nothing
end

"""Whether queued work is due. Caller holds `condition`; this performs no I/O."""
function _flush_due(buffer::CacheBuffer, now::Float64)::Bool
    isempty(buffer.queued) && return false
    at_capacity = buffer.row_limit !== nothing && buffer.queued_rows >= buffer.row_limit
    return buffer.closing || at_capacity ||
        now - buffer.last_flush >= CACHE_BUFFER_FLUSH_INTERVAL
end

"""
Wait for work, detach the active queue, write it, then clear the completed batch.

The condition is held only for waiting, queue transfer, and final state publication. Database I/O
runs without it.
"""
function _flush_loop!(buffer::CacheBuffer{R})::Nothing where {R}
    while true
        writing = lock(buffer.condition) do
            while !_flush_due(buffer, time())
                buffer.closing && isempty(buffer.queued) && return nothing
                remaining = max(
                    eps(Float64),
                    CACHE_BUFFER_FLUSH_INTERVAL - (time() - buffer.last_flush),
                )
                timer = Timer(remaining) do _
                    lock(() -> notify(buffer.condition), buffer.condition)
                end
                try
                    wait(buffer.condition)
                finally
                    close(timer)
                end
            end
            buffer.writing = buffer.queued
            buffer.writing_rows = buffer.queued_rows
            buffer.queued = Dict{Any,BufferMutation{R}}()
            buffer.queued_rows = 0
            buffer.last_flush = time()
            notify(buffer.condition)
            buffer.writing
        end
        writing === nothing && return nothing
        try
            _flush_to_db!(buffer, writing)
        catch error
            lock(buffer.condition) do
                notify(buffer.condition, error; all=true, error=true)
            end
            @error "Cache buffer flush failed" exception=error
            rethrow()
        end
        lock(buffer.condition) do
            buffer.writing = Dict{Any,BufferMutation{R}}()
            buffer.writing_rows = 0
            notify(buffer.condition)
        end
    end
end

"""Resolve one key from a mutation dictionary. Caller holds `condition`."""
function _read_mutation(mutations::Dict, key)
    mutation = get(mutations, key, nothing)
    mutation === nothing && return (false, nothing)
    return mutation.kind === BUFFER_DELETE ? (true, nothing) : (true, mutation.row)
end

"""Normalize one configured key to the prepared statement's positional parameters."""
function _key_parameters(buffer::CacheBuffer, key)::Tuple
    length(buffer.key_columns) == 1 && return (key,)
    key isa Tuple && length(key) == length(buffer.key_columns) ||
        throw(ArgumentError(
            "cache buffer key must be a $(length(buffer.key_columns))-value tuple",
        ))
    return key
end

"""Read one ordinary-table key through queued, writing, then committed state."""
function Base.read(buffer::CacheBuffer{R}, key)::Union{Nothing,R} where {R}
    pending = lock(buffer.condition) do
        _require_open(buffer)
        queued = _read_mutation(buffer.queued, key)
        queued[1] ? queued : _read_mutation(buffer.writing, key)
    end
    pending[1] && return pending[2]
    buffer.read_connection === nothing && return nothing
    rows = collect(DBInterface.execute(buffer.keyed_read, _key_parameters(buffer, key)))
    length(rows) <= 1 ||
        error("CacheBuffer keyed read returned multiple rows for key $(repr(key))")
    committed = isempty(rows) ? nothing : R(only(rows))
    latest = lock(buffer.condition) do
        _require_open(buffer)
        queued = _read_mutation(buffer.queued, key)
        queued[1] ? queued : _read_mutation(buffer.writing, key)
    end
    return latest[1] ? latest[2] : committed
end

"""Apply one mutation dictionary to a complete logical ordinary-table result."""
function _overlay!(rows::Dict, mutations::Dict)::Nothing
    for (key, mutation) in mutations
        if mutation.kind === BUFFER_DELETE
            delete!(rows, key)
        else
            rows[key] = mutation.row
        end
    end
    return nothing
end

"""Read one complete ordinary table and overlay writing followed by queued state."""
function Base.read(buffer::CacheBuffer{R})::Dict{Any,R} where {R}
    buffer.read_connection === nothing && return lock(buffer.condition) do
        _require_open(buffer)
        rows = Dict{Any,R}()
        _overlay!(rows, buffer.writing)
        _overlay!(rows, buffer.queued)
        rows
    end
    while true
        writing, queued = lock(buffer.condition) do
            _require_open(buffer)
            (buffer.writing, buffer.queued)
        end
        rows = Dict{Any,R}()
        for database_row in DBInterface.execute(
            buffer.read_connection,
            "SELECT * FROM $(buffer.table)",
        )
            row = R(database_row)
            key_values = Tuple(getproperty(database_row, column) for column in buffer.key_columns)
            key = length(key_values) == 1 ? only(key_values) : key_values
            rows[key] = row
        end
        stable = lock(buffer.condition) do
            _require_open(buffer)
            buffer.writing === writing && buffer.queued === queued || return false
            _overlay!(rows, buffer.writing)
            _overlay!(rows, buffer.queued)
            true
        end
        stable && return rows
    end
end

"""Clear one ordinary owned table in one transaction. Caller holds `condition`."""
function _clear_db!(buffer::CacheBuffer)::Nothing
    connection = buffer.write_connection
    DBInterface.execute(connection, "BEGIN TRANSACTION")
    try
        DBInterface.execute(connection, "DELETE FROM $(buffer.table)")
        DBInterface.execute(connection, "COMMIT")
    catch
        DBInterface.execute(connection, "ROLLBACK")
        rethrow()
    end
    return nothing
end

"""
Synchronously discard queued work, wait for the current write, and clear the complete owned store.

The condition remains held during the cold database clear so no operation can enter.
"""
function clear!(buffer::CacheBuffer)::Nothing
    lock(buffer.condition)
    try
        _require_writable(buffer)
        empty!(buffer.queued)
        buffer.queued_rows = 0
        while !isempty(buffer.writing)
            wait(buffer.condition)
            _require_writable(buffer)
        end
        if buffer.write_connection !== nothing
            _clear_db!(buffer)
        end
        notify(buffer.condition)
        return nothing
    finally
        unlock(buffer.condition)
    end
end

"""Permanently close a live buffer after its final queued batch is written."""
function close!(buffer::CacheBuffer)::Nothing
    lock(buffer.condition) do
        buffer.closing && error("CacheBuffer is already closed")
        buffer.closing = true
        notify(buffer.condition)
    end
    task_error::Union{Nothing,Exception} = nothing
    if buffer.flush_task !== nothing
        try
            wait(buffer.flush_task)
        catch error
            task_error = error
        end
    end
    try
        buffer.read_connection === nothing || DBInterface.close!(buffer.read_connection)
    finally
        buffer.write_connection === nothing || DBInterface.close!(buffer.write_connection)
    end
    task_error === nothing || throw(task_error)
    return nothing
end
