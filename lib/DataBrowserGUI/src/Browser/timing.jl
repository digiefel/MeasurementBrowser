import TimerOutputs

# Always-on render-loop instrumentation.
#
# `MAIN_TIMER` is the single TimerOutput that `@timed` records into. It is written
# ONLY from the main (GUI) task — the render loop and everything it calls
# synchronously. That invariant is what makes it lock-free: TimerOutputs keeps an
# active-section stack per timer with no synchronization, so a `@timed` section
# entered off the main task would corrupt it. Do not use `@timed` from a spawned
# task; the multi-task engine pipeline uses the dev-only `@timed_dbg` instead.

const MAIN_TIMER = TimerOutputs.TimerOutput("main task")

"""
    @timed [label] expr

Time `expr` in the always-on main-task timer (`MAIN_TIMER`), recording call count,
time, and allocations, and nesting automatically inside any enclosing `@timed`
section. Like `@timed_dbg` but always on and with no `level`.

Forms:

    @timed some_call(...)          # label is the expression text
    @timed "label" begin ... end   # explicit label

Main-task-only (see `MAIN_TIMER`). Internal; pulled in with `using ...: @timed`.
"""
macro timed(args...)
    rest = collect(args)
    label = nothing
    if length(rest) >= 2 && rest[1] isa AbstractString
        label = String(rest[1])
        rest = rest[2:end]
    end
    length(rest) == 1 || error(
        "@timed expects a single expression, optionally preceded by a string label",
    )
    expr = rest[1]
    lbl = label === nothing ? string(expr) : label
    return :(TimerOutputs.@timeit MAIN_TIMER $lbl $(esc(expr)))
end
