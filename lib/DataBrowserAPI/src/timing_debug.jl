# Internal debug-timing instrumentation.
#
# This file provides ONLY the zero-cost `@timed_dbg` annotation and the dormant
# hooks it calls. It has no dependency on TimerOutputs and no timing machinery.
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
#
# The label-derivation helpers `_timed_label`/`_timed_callee` are dependency-free
# expression parsing, shared by both `@timed_dbg` (here) and the always-on `@timed`
# macro in DataBrowserGUI.

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

# Derive a section label from the timed expression when none is given explicitly.
# For a call `f(...)` this is `"f"`; for `Mod.f(...)` it is `"f"`. Anything else
# falls back to the stringified expression (callers timing a `begin ... end`
# block are expected to pass an explicit label). Shared by `@timed_dbg` and `@timed`.
_timed_label(sym::Symbol) = String(sym)
function _timed_label(ex::Expr)
    if ex.head === :call
        return _timed_callee(ex.args[1])
    elseif ex.head === :. && length(ex.args) == 2 && ex.args[2] isa QuoteNode
        return String(ex.args[2].value)
    end
    return string(ex)
end
_timed_label(x) = string(x)

_timed_callee(sym::Symbol) = String(sym)
function _timed_callee(ex::Expr)
    if ex.head === :. && length(ex.args) == 2 && ex.args[2] isa QuoteNode
        return String(ex.args[2].value)
    end
    return string(ex)
end
_timed_callee(x) = string(x)

"""
    @timed_dbg [level=N] [label] expr

Time `expr` under section `label` when instrumentation is active. Internal;
compiles away to just `expr` unless `DataBrowserProfiling` is loaded.

Forms:

    @timed_dbg expr                  # level 1, label derived from the call
    @timed_dbg "label" expr          # level 1, explicit label
    @timed_dbg level=2 expr          # level 2, label derived from the call
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
    lbl = label === nothing ? _timed_label(expr) : label

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
