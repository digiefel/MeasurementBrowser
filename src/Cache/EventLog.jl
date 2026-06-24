# Crash-durable operation log. Opt-in per workspace via open_workspace(...; event_log=path). Each instrumented cache operation
# writes a `begin` row before it runs and an `end` row after, each flushed to disk immediately under a
# lock. The point is survivability: if the process segfaults mid-operation, the file on disk still ends
# with that operation's `begin` and no matching `end` — so a crash with no error names the exact call
# (and the cache size at the time) that killed it. It also doubles as timing data: plot dur_ms vs
# prior_items per op to see whether a cost grows with the cache (an O(N^2) scan) or stays flat.
#
# This is a diagnostic mode, not a normal-run path: the per-line flush under a global lock serializes
# instrumented ops across threads. That is the cost of durability and is acceptable when hunting a bug.

const EVENT_LOG = Ref{Union{Nothing,IO}}(nothing)
const EVENT_LOCK = ReentrantLock()
const EVENT_T0 = Ref{UInt64}(0)
# Cumulative items written to payload tables so far, used as the x-axis ("how full is the cache").
const ITEMS_WRITTEN = Base.Threads.Atomic{Int}(0)
# The file the live/last build wrote to. Recording is requested per workspace via open_workspace's
# `event_log=` argument; the build lifecycle opens and closes the file here.
const EVENT_LOG_PATH = Ref{String}("")

"""Open the crash-durable event log at `path`, truncating any previous file."""
function start_event_log!(path::AbstractString)::Nothing
    stop_event_log!()
    io = open(String(path), "w")
    println(io, "t_s,thread,phase,op,prior_items,rows,dur_ms")
    flush(io)
    EVENT_T0[] = time_ns()
    ITEMS_WRITTEN[] = 0
    EVENT_LOG_PATH[] = String(path)
    EVENT_LOG[] = io
    return nothing
end

"""Close the event log if open. Idempotent and safe to call on every idle poll."""
function stop_event_log!()::Nothing
    io = EVENT_LOG[]
    EVENT_LOG[] = nothing
    io === nothing && return nothing
    try
        close(io)
    catch
    end
    return nothing
end

"""The path the active (or most recent) build's event log was written to, or empty if none yet."""
current_event_log_path()::String = EVENT_LOG_PATH[]

"""Write one event row and flush it to the OS so it survives a crash of this process."""
function _event_row(phase::Symbol, op::Symbol, rows::Integer, dur_ms::Union{Nothing,Float64})::Nothing
    io = EVENT_LOG[]
    io === nothing && return nothing
    t = round((time_ns() - EVENT_T0[]) / 1e9; digits=6)
    prior = ITEMS_WRITTEN[]
    tid = Base.Threads.threadid()
    duration = dur_ms === nothing ? "" : string(round(dur_ms; digits=4))
    lock(EVENT_LOCK) do
        println(io, "$t,$tid,$phase,$op,$prior,$rows,$duration")
        flush(io)
    end
    return nothing
end

"""
Wrap one cache operation so it logs `begin` before and `end` (with duration) after.

`rows` is this operation's own size; `prior_items` (logged automatically) is the cache size when it
ran. A `begin` with no matching `end` at the tail of the file is the operation that crashed.
"""
macro cache_event(op, rows, expr)
    quote
        local _op = $(esc(op))
        local _rows = $(esc(rows))
        if EVENT_LOG[] === nothing
            $(esc(expr))                       # zero overhead when not recording
        else
            _event_row(:begin, _op, _rows, nothing)
            local _t = time_ns()
            local _result = $(esc(expr))
            _event_row(:end, _op, _rows, (time_ns() - _t) / 1e6)
            _result
        end
    end
end
