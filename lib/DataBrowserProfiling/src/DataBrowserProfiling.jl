"""Debug timing summaries and Julia sampling-profile helpers."""
module DataBrowserProfiling

using Profile
using Printf: @printf
using TimerOutputs

const SAMPLING_TOTAL_BUFFER_SIZE = 2_000_000
const SAMPLING_DELAY_SECONDS = 0.01
const DEBUG_TIMER_TLS_KEY = gensym(:databrowser_debug_timer)

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

"""Saved Julia sampling result used by the plot diagnostic UI."""
struct SamplingProfile
    total_samples::Int
    delay_seconds::Float64
    truncated::Bool
    rows::Vector{SamplingProfileRow}
end

"""One explicit debug-timing measurement owned by a benchmark or diagnostic caller."""
mutable struct DebugTimings
    start_ns::UInt64
    stop_ns::UInt64
    generation::UInt64
    active::Base.Threads.Atomic{Bool}
    timers::Vector{TimerOutput}
    lock::ReentrantLock
    merged::Union{Nothing,TimerOutput}
end

DebugTimings(; start_ns::UInt64=time_ns()) =
    DebugTimings(
        start_ns,
        UInt64(0),
        UInt64(1),
        Base.Threads.Atomic{Bool}(false),
        TimerOutput[],
        ReentrantLock(),
        nothing,
    )

const ACTIVE_DEBUG_TIMINGS =
    Base.ScopedValues.ScopedValue{Union{Nothing,DebugTimings}}(nothing)
const ACTIVE_LOCK = ReentrantLock()
const ACTIVE_MEASUREMENT = Ref{Union{Nothing,DebugTimings}}(nothing)
const DEBUG_TIMINGS_ACTIVE = Base.Threads.Atomic{Bool}(false)

struct DebugTimingRow
    label::String
    calls::Int
    total_elapsed_ns::Int64
    average_elapsed_ns::Float64
    process_allocation_delta_bytes::Int64
    average_process_allocation_delta_bytes::Float64
end

struct TaskDebugTimer
    timings::DebugTimings
    generation::UInt64
    timer::TimerOutput
end

"""Reset one explicit measurement before entering `with_debug_timings`."""
function reset_debug_timings!(timings::DebugTimings; start_ns::UInt64=time_ns())::DebugTimings
    lock(timings.lock) do
        timings.start_ns = start_ns
        timings.stop_ns = 0
        timings.generation += 1
        empty!(timings.timers)
        timings.merged = nothing
    end
    return timings
end

function _task_timer(timings::DebugTimings)::TimerOutput
    tls = task_local_storage()
    cached = get(tls, DEBUG_TIMER_TLS_KEY, nothing)
    if cached isa TaskDebugTimer && cached.timings === timings &&
       cached.generation == timings.generation
        return cached.timer
    end
    timer = TimerOutput()
    tls[DEBUG_TIMER_TLS_KEY] = TaskDebugTimer(timings, timings.generation, timer)
    lock(timings.lock) do
        push!(timings.timers, timer)
        timings.merged = nothing
    end
    return timer
end

"""Return the current task's timer, or `nothing` outside an explicit measurement."""
function debug_timer()::Union{Nothing,TimerOutput}
    cached = get(task_local_storage(), DEBUG_TIMER_TLS_KEY, nothing)
    if cached isa TaskDebugTimer && cached.timings.active[] &&
       cached.generation == cached.timings.generation
        return cached.timer
    end
    DEBUG_TIMINGS_ACTIVE[] || return nothing
    timings = ACTIVE_DEBUG_TIMINGS[]
    timings === nothing && (timings = lock(ACTIVE_LOCK) do
        ACTIVE_MEASUREMENT[]
    end)
    timings === nothing && return nothing
    return _task_timer(timings::DebugTimings)
end

function _time_dbg_label(ex)
    ex isa Expr && ex.head === :call || throw(ArgumentError(
        "@time_dbg expects a call; use @time_dbg \"label\" begin ... end for a compound expression",
    ))
    callee = ex.args[1]
    callee isa Symbol && return String(callee)
    callee isa Expr && callee.head === :. || throw(ArgumentError(
        "@time_dbg needs an explicit label for this callee",
    ))
    name = callee.args[end]
    name isa QuoteNode && (name = name.value)
    name isa Symbol || throw(ArgumentError("@time_dbg needs an explicit label for this callee"))
    return String(name)
end

"""Time one debug-only call with a label derived from its callee."""
macro time_dbg(args...)
    label, ex = if length(args) == 1
        (_time_dbg_label(args[1]), args[1])
    elseif length(args) == 2 && args[1] isa String
        (args[1], args[2])
    else
        throw(ArgumentError("use @time_dbg call(...) or @time_dbg \"label\" expression"))
    end
    timer = gensym(:debug_timer)
    return quote
        local $timer = $(GlobalRef(DataBrowserProfiling, :debug_timer))()
        if $timer === nothing
            $(esc(ex))
        else
            $(Expr(:macrocall,
                GlobalRef(TimerOutputs, Symbol("@timeit")),
                __source__,
                timer,
                label,
                esc(ex),
            ))
        end
    end
end

"""Run `work` while instrumented calls record into the explicit `timings` object."""
function with_debug_timings(work::Function, timings::DebugTimings)
    lock(ACTIVE_LOCK) do
        ACTIVE_MEASUREMENT[] === nothing || error("A DebugTimings measurement is already active")
        ACTIVE_MEASUREMENT[] = timings
        timings.active[] = true
        DEBUG_TIMINGS_ACTIVE[] = true
    end
    try
        return Base.ScopedValues.with(ACTIVE_DEBUG_TIMINGS => timings) do
            work()
        end
    finally
        lock(ACTIVE_LOCK) do
            DEBUG_TIMINGS_ACTIVE[] = false
            timings.active[] = false
            ACTIVE_MEASUREMENT[] = nothing
        end
        lock(timings.lock) do
            timings.stop_ns = time_ns()
        end
    end
end

"""Merge task-owned timers after the caller has reached its intended quiescent point."""
function collect_debug_timings!(timings::DebugTimings)::TimerOutput
    lock(timings.lock) do
        merged = TimerOutput("DataBrowser debug timings")
        for timer in timings.timers
            merge!(merged, timer)
        end
        timings.merged = merged
        return merged
    end
end

function _flatten_rows!(rows::Vector{DebugTimingRow}, timer::TimerOutput)::Nothing
    for (label, child) in timer.inner_timers
        calls = TimerOutputs.ncalls(child)
        elapsed = TimerOutputs.time(child)
        allocated = TimerOutputs.allocated(child)
        push!(rows, DebugTimingRow(
            label,
            calls,
            elapsed,
            calls == 0 ? 0.0 : elapsed / calls,
            allocated,
            calls == 0 ? 0.0 : allocated / calls,
        ))
    end
    return nothing
end

"""Return flattened, label-level rows for a completed measurement."""
function debug_timing_rows(timings::DebugTimings)::Vector{DebugTimingRow}
    merged = timings.merged === nothing ? collect_debug_timings!(timings) : timings.merged
    flattened = TimerOutputs.flatten(merged)
    rows = DebugTimingRow[]
    _flatten_rows!(rows, flattened)
    sort!(rows; by=row -> row.total_elapsed_ns, rev=true)
    return rows
end

"""Write text and CSV debug timing summaries under `outdir`."""
function write_debug_timings(outdir::AbstractString, timings::DebugTimings)::Nothing
    timings.stop_ns >= timings.start_ns || error(
        "Run with_debug_timings before writing its results",
    )
    mkpath(outdir)
    rows = debug_timing_rows(timings)
    open(joinpath(outdir, "debug_timings.txt"), "w") do io
        elapsed = (timings.stop_ns - timings.start_ns) / 1e9
        @printf(io, "Measurement wall time: %.6f s\n\n", elapsed)
        println(io, "Times are inclusive: a marked parent call includes its marked child calls.")
        println(io, "Calls on different tasks can overlap, so row totals can exceed wall time.")
        println(io, "Allocation changes use Julia's process-wide counter and can overlap too.\n")
        label_width = max(5, maximum(length(row.label) for row in rows; init=0))
        println(io,
            rpad("label", label_width), "  ",
            lpad("calls", 8), "  ",
            lpad("total seconds", 14), "  ",
            lpad("average ms", 14), "  ",
            lpad("allocation change (bytes)", 25),
        )
        for row in rows
            @printf(
                io,
                "%s  %8d  %14.6f  %14.6f  %25d\n",
                rpad(row.label, label_width),
                row.calls,
                row.total_elapsed_ns / 1e9,
                row.average_elapsed_ns / 1e6,
                row.process_allocation_delta_bytes,
            )
        end
    end
    open(joinpath(outdir, "debug_timings.csv"), "w") do io
        println(io, "label,calls,total_elapsed_ns,average_elapsed_ns,process_allocation_delta_bytes,average_process_allocation_delta_bytes")
        for row in debug_timing_rows(timings)
            label = replace(row.label, '"' => "\"\"")
            println(io, '"', label, "\",", row.calls, ',', row.total_elapsed_ns,
                ',', row.average_elapsed_ns, ',', row.process_allocation_delta_bytes,
                ',', row.average_process_allocation_delta_bytes)
        end
    end
    return nothing
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

export DebugTimings, DebugTimingRow, @time_dbg, with_debug_timings, reset_debug_timings!,
    collect_debug_timings!, debug_timing_rows, write_debug_timings,
    SamplingProfile, SamplingProfileRow, start_sampling!, stop_sampling!, cancel_sampling!, process_rss_bytes

end # module
