"""
Profiling and timing infrastructure shared by the whole engine.

Hot regions across the pipeline are wrapped with `@timeit_debug TIMER "label" expr`. These compile to
nothing until [`enable!`](@ref) turns on debug timings, so the instrumentation is permanent in the
source yet costs nothing in normal runs. Enable it, exercise a workload, then [`report`](@ref) the
accumulated per-region time, call counts, and allocations.

`@timeit_debug` writes are not thread-safe into a shared timer, so only sequential regions (cache
reads/writes, index load, the analysis loop) are instrumented. The always-on scan profile records
source-item callback phases separately. An explicit profiled rebuild uses Julia's bounded sampling
profiler and reduces its stacks to source-line rows for the Performance window.
"""
module Profiling

using Profile
using TimerOutputs
export TIMER, @timeit_debug

"""One source line in a completed sampling profile."""
struct SamplingProfileRow
    samples::Int
    self_samples::Int
    function_name::String
    file::String
    line::Int
end

"""Bounded sampling result retained by the workspace after a profiled rebuild."""
struct SamplingProfile
    total_samples::Int
    delay_seconds::Float64
    truncated::Bool
    rows::Vector{SamplingProfileRow}
end

const SAMPLING_TOTAL_BUFFER_SIZE = 2_000_000
const SAMPLING_DELAY_SECONDS = 0.01

"""The one timer every instrumented region in the package writes to."""
const TIMER = TimerOutput()

"""Modules that contain `@timeit_debug` regions; [`enable!`](@ref) flips them together."""
const INSTRUMENTED_MODULES = Module[]

"""Register `m` as carrying `@timeit_debug` regions. Called once from each instrumented module."""
function register_instrumented!(m::Module)::Nothing
    m in INSTRUMENTED_MODULES || push!(INSTRUMENTED_MODULES, m)
    return nothing
end

"""Turn on debug timings for every instrumented module and clear the timer."""
function enable!()::Nothing
    reset_timer!(TIMER)
    for m in INSTRUMENTED_MODULES
        TimerOutputs.enable_debug_timings(m)
    end
    return nothing
end

"""Turn off debug timings for every instrumented module."""
function disable!()::Nothing
    for m in INSTRUMENTED_MODULES
        TimerOutputs.disable_debug_timings(m)
    end
    return nothing
end

"""Clear all accumulated timings without changing whether timing is enabled."""
reset!()::Nothing = (reset_timer!(TIMER); nothing)

"""Print the accumulated per-region timings (time, calls, allocations)."""
function report(io::IO=stdout; kwargs...)::Nothing
    print_timer(io, TIMER; kwargs...)
    println(io)
    return nothing
end

"""Run `f()` with debug timings enabled, print the report, then restore the previous state."""
function profile(f; kwargs...)
    enable!()
    try
        result = f()
        report(; kwargs...)
        return result
    finally
        disable!()
    end
end

"""Start a bounded, all-thread CPU sampling profile for one explicit workspace operation."""
function start_sampling!()::Nothing
    Profile.is_running() && error(
        "Cannot profile the rebuild because Julia's CPU profiler is already running",
    )
    buffer_per_thread = max(
        50_000,
        cld(SAMPLING_TOTAL_BUFFER_SIZE, Base.Threads.nthreads()),
    )
    Profile.init(n=buffer_per_thread, delay=SAMPLING_DELAY_SECONDS)
    Profile.clear()
    Profile.start_timer()
    return nothing
end

"""Stop sampling and reduce the raw stacks to flat source-line rows."""
function stop_sampling!()::SamplingProfile
    running = Profile.is_running()
    truncated = Profile.is_buffer_full()
    (running || truncated) || error(
        "Cannot finish the rebuild profile because sampling is not running",
    )
    running && Profile.stop_timer()
    data, line_info = Profile.retrieve(include_meta=false, limitwarn=false)
    flat_data, flat_info = Profile.flatten(data, line_info)

    counts = Dict{Tuple{String,String,Int},Int}()
    self_counts = Dict{Tuple{String,String,Int},Int}()
    frames = Tuple{String,String,Int}[]
    total_samples = 0
    for instruction_pointer in flat_data
        if instruction_pointer == 0
            isempty(frames) && continue
            total_samples += 1
            for frame in unique(frames)
                counts[frame] = get(counts, frame, 0) + 1
            end
            # Julia stores each sampled stack from the active leaf toward the root.
            leaf = frames[1]
            self_counts[leaf] = get(self_counts, leaf, 0) + 1
            empty!(frames)
            continue
        end
        frame = flat_info[instruction_pointer]
        frame.from_c && continue
        push!(frames, (String(frame.func), String(frame.file), frame.line))
    end

    rows = SamplingProfileRow[
        SamplingProfileRow(
            samples,
            get(self_counts, frame, 0),
            frame[1],
            frame[2],
            frame[3],
        )
        for (frame, samples) in counts
    ]
    sort!(rows; by=row -> (row.self_samples, row.samples), rev=true)
    result = SamplingProfile(total_samples, SAMPLING_DELAY_SECONDS, truncated, rows)

    # Retain only the compact report. Raw instruction pointers and line dictionaries can otherwise
    # remain reachable until a later full collection, directly after the cache build's own peak.
    empty!(data)
    empty!(flat_data)
    empty!(line_info)
    empty!(flat_info)
    empty!(frames)
    empty!(counts)
    empty!(self_counts)
    Profile.clear()
    GC.gc(true)
    return result
end

"""Stop an interrupted sampling run without producing a report."""
function cancel_sampling!()::Nothing
    Profile.is_running() && Profile.stop_timer()
    Profile.clear()
    GC.gc(true)
    return nothing
end

end # module Profiling
