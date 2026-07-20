"""
Dev-only DataBrowser instrumentation.

Loading this package turns on `@time_dbg` timing (the macro itself lives in
DataBrowserAPI and is dormant until enabled): on load we install the timing hooks
into DataBrowserAPI and raise its active profiling level, then accumulate
per-task `TimerOutput` segments into one shared master timer.

The collector is task-safe by construction: each Julia task records into its own
`TimerOutput`, and a completed segment is merged into the master under a lock when
that task's outermost `@time_dbg` section exits. All master access goes through
the lock, so a snapshot can never race a merge.

Public entry points return native `TimerOutputs.TimerOutput` objects; display and
analysis are delegated to TimerOutputs (show, flatten, merge, Tables.jl, ...):

    snapshot_debug_timings()   # independent copy of accumulated timings so far
    take_debug_timings!()      # like snapshot, but also resets the master
    finish_debug_timings!()    # stop recording and return the final timings
    reset_debug_timings!()     # clear the master and (re)enable recording

The Julia sampling-profiler helpers and the process-RSS helper also live here but
share no state with the collector above.
"""
module DataBrowserProfiling

using Profile
using TimerOutputs
import DataBrowserAPI

# ===========================================================================
# Instrumentation collector
# ===========================================================================

const _TLS_KEY = gensym(:databrowser_time_dbg)
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
        _time_dbg_begin(label) = $(beginimpl)(label)
        _time_dbg_end(token) = $(endimpl)(token)
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

# ===========================================================================
# Julia sampling profiler (dev-only; shares no state with the collector above)
# ===========================================================================

const SAMPLING_TOTAL_BUFFER_SIZE = 2_000_000
const SAMPLING_DELAY_SECONDS = 0.01

struct MachTaskBasicInfo
    virtual_size::UInt64
    resident_size::UInt64
    resident_size_max::UInt64
    user_seconds::Int32
    user_microseconds::Int32
    system_seconds::Int32
    system_microseconds::Int32
    policy::Int32
    suspend_count::Int32
end

"""One reduced source line from Julia's sampling profiler."""
struct SamplingProfileRow
    samples::Int
    self_samples::Int
    function_name::String
    file::String
    line::Int
end

"""Saved Julia sampling result used by diagnostic tooling."""
struct SamplingProfile
    total_samples::Int
    delay_seconds::Float64
    truncated::Bool
    rows::Vector{SamplingProfileRow}
end

"""Start Julia's bounded all-thread CPU sampling profiler."""
function start_sampling!()::Nothing
    Profile.is_running() && error("Julia's CPU profiler is already running")
    Profile.init(n=max(50_000, cld(SAMPLING_TOTAL_BUFFER_SIZE, Base.Threads.nthreads())),
        delay=SAMPLING_DELAY_SECONDS)
    Profile.clear()
    Profile.start_timer()
    return nothing
end

function stop_sampling!()::SamplingProfile
    running = Profile.is_running()
    truncated = Profile.is_buffer_full()
    (running || truncated) || error("Julia's CPU profiler is not running")
    running && Profile.stop_timer()
    data, info = Profile.retrieve(include_meta=false, limitwarn=false)
    flat, flat_info = Profile.flatten(data, info)
    counts = Dict{Tuple{String,String,Int},Int}()
    self_counts = Dict{Tuple{String,String,Int},Int}()
    frames = Tuple{String,String,Int}[]
    total_samples = 0
    for ip in flat
        if ip == 0
            isempty(frames) && continue
            total_samples += 1
            for frame in unique(frames)
                counts[frame] = get(counts, frame, 0) + 1
            end
            self_counts[frames[1]] = get(self_counts, frames[1], 0) + 1
            empty!(frames)
        else
            frame = flat_info[ip]
            frame.from_c || push!(
                frames,
                (String(frame.func), String(frame.file), frame.line),
            )
        end
    end
    rows = SamplingProfileRow[
        SamplingProfileRow(samples, get(self_counts, frame, 0), frame...)
        for (frame, samples) in counts
    ]
    sort!(rows; by=row -> (row.self_samples, row.samples), rev=true)
    Profile.clear()
    return SamplingProfile(total_samples, SAMPLING_DELAY_SECONDS, truncated, rows)
end

function cancel_sampling!()::Nothing
    Profile.is_running() && Profile.stop_timer()
    Profile.clear()
    return nothing
end

"""Return current process resident bytes."""
function process_rss_bytes()::Int64
    if Sys.isapple()
        info = Ref(MachTaskBasicInfo(0, 0, 0, 0, 0, 0, 0, 0, 0))
        count = Ref{UInt32}(UInt32(div(sizeof(MachTaskBasicInfo), sizeof(Int32))))
        task = ccall(:mach_task_self, UInt32, ())
        result = ccall(:task_info, Int32,
            (UInt32, Int32, Ref{MachTaskBasicInfo}, Ref{UInt32}), task, 20, info, count)
        result == 0 || error("mach task_info failed while sampling RSS with code $result")
        return Int64(info[].resident_size)
    elseif Sys.islinux()
        fields = split(read("/proc/self/statm", String))
        return parse(Int64, fields[2]) * Int64(Sys.page_size())
    end
    error("Current RSS sampling is unsupported on $(Sys.KERNEL)")
end

export snapshot_debug_timings, take_debug_timings!, finish_debug_timings!,
    reset_debug_timings!,
    SamplingProfile, SamplingProfileRow, start_sampling!, stop_sampling!,
    cancel_sampling!, process_rss_bytes

end # module
