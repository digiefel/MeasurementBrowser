"""
Dev-only DataBrowser instrumentation.

Loading this package turns on `@timed_dbg` timing (the macro itself lives in
DataBrowserAPI and is dormant until enabled): on load we install the timing hooks
into DataBrowserAPI and raise its active profiling level, then accumulate
per-task `TimerOutput` segments into one shared master timer.

The collector is task-safe by construction: each Julia task records into its own
`TimerOutput`, and a completed segment is merged into the master under a lock when
that task's outermost `@timed_dbg` section exits. All master access goes through
the lock, so a snapshot can never race a merge.

Public entry points return native `TimerOutputs.TimerOutput` objects; display and
analysis are delegated to TimerOutputs (show, flatten, merge, Tables.jl, ...):

    snapshot_debug_timings()   # independent copy of accumulated timings so far
    take_debug_timings!()      # like snapshot, but also resets the master
    finish_debug_timings!()    # stop recording and return the final timings
    reset_debug_timings!()     # clear the master and (re)enable recording
"""
module DataBrowserProfiling

using TimerOutputs
import DataBrowserAPI

# ===========================================================================
# Instrumentation collector
# ===========================================================================

const _TLS_KEY = gensym(:databrowser_timed_dbg)
const DEFAULT_LEVEL = 1

# Per-task timing state. One task owns one TimerOutput segment and never shares
# it, so it records without a lock. `depth` tracks nesting so we know when the
# outermost section closes.
mutable struct TaskTimingContext
    timer::TimerOutput
    depth::Int
end

# The single accumulated timer. Every access — merge on submit, snapshot, take,
# reset — is serialized by MASTER_LOCK, making this the sole mutation point and
# guaranteeing no reader races a merge.
const MASTER = Base.RefValue{TimerOutput}(TimerOutput("DataBrowser debug timings"))
const MASTER_LOCK = ReentrantLock()

# Whether new outermost sections are admitted. Set on load; cleared by
# finish_debug_timings!. Sections already in flight always finish normally.
const RECORDING = Base.Threads.Atomic{Bool}(false)

function _submit!(timer::TimerOutput)
    lock(MASTER_LOCK) do
        merge!(MASTER[], timer)
    end
    return nothing
end

# Hook implementations installed into DataBrowserAPI on load. `_begin` returns an
# inert `nothing` token when not recording so the paired `_end` is a no-op.
function _begin(label)
    RECORDING[] || return nothing
    tls = task_local_storage()
    ctx = get(tls, _TLS_KEY, nothing)
    if ctx === nothing
        ctx = TaskTimingContext(TimerOutput(), 0)
        tls[_TLS_KEY] = ctx
    end
    ctx = ctx::TaskTimingContext
    ctx.depth += 1
    section = begin_timed_section!(ctx.timer, label)
    return (ctx, section)
end

function _end(token)
    token === nothing && return nothing
    ctx, section = token
    end_timed_section!(ctx.timer, section)
    ctx.depth -= 1
    if ctx.depth == 0
        _submit!(ctx.timer)
        delete!(task_local_storage(), _TLS_KEY)
    end
    return nothing
end

# Install hook implementations into DataBrowserAPI and raise its profiling level.
# Redefining these methods invalidates and recompiles the annotated call sites so
# their timing branch goes live. Done at load time, never during precompilation.
function _enable!()
    beginimpl = _begin
    endimpl = _end
    Core.eval(DataBrowserAPI, quote
        _timed_dbg_begin(label) = $(beginimpl)(label)
        _timed_dbg_end(token) = $(endimpl)(token)
        profile_level() = $(DEFAULT_LEVEL)
    end)
    RECORDING[] = true
    return nothing
end

"""
    snapshot_debug_timings() -> TimerOutput

An independent copy of the debug timings accumulated so far. Sections still
running (whose task-owned segment has not yet been submitted) are excluded; they
appear in a later snapshot once they complete.
"""
function snapshot_debug_timings()
    lock(MASTER_LOCK) do
        snap = TimerOutput("DataBrowser debug timings")
        merge!(snap, MASTER[])
        return snap
    end
end

"""
    take_debug_timings!() -> TimerOutput

Return the accumulated debug timings and reset the master to empty so a fresh
interval begins. A completed outermost section belongs to whichever interval it
reaches the master in.
"""
function take_debug_timings!()
    lock(MASTER_LOCK) do
        taken = MASTER[]
        MASTER[] = TimerOutput("DataBrowser debug timings")
        return taken
    end
end

"""
    finish_debug_timings!() -> TimerOutput

Stop admitting new sections and return the final accumulated debug timings.
Sections still in flight are no longer recorded. Intended to be called once at
shutdown, after the GUI and workspace are closed.
"""
function finish_debug_timings!()
    RECORDING[] = false
    return take_debug_timings!()
end

"""
    reset_debug_timings!() -> Nothing

Discard all accumulated debug timings and (re)enable recording. Useful from the
REPL or tests to isolate a fresh interval.
"""
function reset_debug_timings!()
    lock(MASTER_LOCK) do
        MASTER[] = TimerOutput("DataBrowser debug timings")
    end
    RECORDING[] = true
    return nothing
end

function __init__()
    _enable!()
    return nothing
end

export snapshot_debug_timings, take_debug_timings!, finish_debug_timings!,
    reset_debug_timings!

end # module
