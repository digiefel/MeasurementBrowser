# Internal debug-timing instrumentation.
#
# When the dev-only `DataBrowserProfiling` package is not loaded, `profile_level()`
# returns 0, so every `@timed_dbg level=N ...` (with N >= 1) has a statically-false
# guard and the compiler folds it away to just the wrapped expression — the
# `_timed_dbg_begin`/`_timed_dbg_end` stubs below are never called and add nothing.
#
# `DataBrowserProfiling` activates instrumentation at runtime by (a) redefining
# `profile_level()` to a positive value, which recompiles annotated methods with
# the timing branch live, and (b) overriding the two hook methods to record into
# per-task TimerOutputs. None of that lives here.
#
# `@timed_dbg` is INTERNAL: it is deliberately not exported. Consuming modules pull
# it in with `using DataBrowserAPI: @timed_dbg`; it never reaches DataBrowser's
# public API.

"""
    profile_level() -> Int

Active instrumentation level. `0` (the default) disables every `@timed_dbg`
section and lets the compiler remove it entirely. `DataBrowserProfiling` raises
this at runtime to switch instrumentation on. A `@timed_dbg level=N` section
records only when `profile_level() >= N`, so higher levels can gate chattier
annotations behind an explicit opt-in.
"""
profile_level() = 0

# Dormant timing hooks. `DataBrowserProfiling` overrides these when it loads.
# While `profile_level()` is 0 they are unreachable (dead-code eliminated at the
# call site), so their bodies only matter once instrumentation is enabled.
_timed_dbg_begin(@nospecialize(label)) = nothing
_timed_dbg_end(@nospecialize(token)) = nothing

"""
    @timed_dbg [level=N] [label] expr

Time `expr` under section `label` when instrumentation is active. Internal;
compiles away to just `expr` unless `DataBrowserProfiling` is loaded.

Forms:

    @timed_dbg expr                  # level 1, label is the expression text
    @timed_dbg "label" expr          # level 1, explicit label
    @timed_dbg level=2 expr          # level 2, label is the expression text
    @timed_dbg level=2 "label" expr  # level 2, explicit label

`level` and `label` are independent. The section records only when
`profile_level() >= level`; otherwise the guard is statically false and the whole
annotation is removed by the compiler, so timing has zero cost when disabled.
"""
macro timed_dbg(args...)
    level = 1
    label = nothing
    rest = collect(args)

    if !isempty(rest) && rest[1] isa Expr && rest[1].head === :(=) && rest[1].args[1] === :level
        level = rest[1].args[2]
        rest = rest[2:end]
    end

    if length(rest) >= 2 && rest[1] isa AbstractString
        label = String(rest[1])
        rest = rest[2:end]
    end

    length(rest) == 1 || error(
        "@timed_dbg expects a single expression, optionally preceded by `level=N` and a string label",
    )
    expr = rest[1]
    lbl = label === nothing ? string(expr) : label

    return quote
        if profile_level() >= $(esc(level))
            local _timed_dbg_token = _timed_dbg_begin($lbl)
            try
                $(esc(expr))
            finally
                _timed_dbg_end(_timed_dbg_token)
            end
        else
            $(esc(expr))
        end
    end
end
