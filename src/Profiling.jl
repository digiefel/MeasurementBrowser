"""
Profiling and timing infrastructure shared by the whole engine.

Hot regions across the pipeline are wrapped with `@timeit_debug TIMER "label" expr`. These compile to
nothing until [`enable!`](@ref) turns on debug timings, so the instrumentation is permanent in the
source yet costs nothing in normal runs. Enable it, exercise a workload, then [`report`](@ref) the
accumulated per-region time, call counts, and allocations.

`@timeit_debug` writes are not thread-safe into a shared timer, so only sequential regions (cache
reads/writes, index load, the analysis loop) are instrumented; the parallel scan's per-file parse
cost is captured separately by the per-kind scan profile (`scan_profile_summary`). When profiling,
drive overlapping phases one at a time.
"""
module Profiling

using TimerOutputs
export TIMER, @timeit_debug

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

end # module Profiling
