# Profiling

DataBrowser is instrumented for timing with `@timed_dbg` markers whose data is aggregated
into native [TimerOutputs](https://github.com/KristofferC/TimerOutputs.jl) timing trees.
Instrumentation is **dev-only and zero-cost when off**: production builds do not depend on
the profiler and the markers compile away.

This document describes how the instrumentation works and how to add it. For the interactive
investigation workflow — drive a workload, read the tree, iterate with Revise — see the
`databrowser-profiling` skill.

> **Two timers, don't confuse them.** This page covers the **dev-only, multi-task** `@timed_dbg`
> engine profiler. There is a separate **always-on** `@timed` timer (`MAIN_TIMER`) that records the
> GUI main-task render loop and is shown live in the Performance window's **Timings** tab (see
> [gui.md](gui.md)). `@timed` is main-task-only and lock-free; `@timed_dbg` spans worker tasks and
> stays off unless `DataBrowserProfiling` is loaded. Same TimerOutputs backend, different scope.

## The `@timed_dbg` marker

`@timed_dbg` lives in `DataBrowserAPI` (`lib/DataBrowserAPI/src/timing_debug.jl`) and is
internal — it is not exported to DataBrowser's public API. Engine code pulls it in with
`using DataBrowserAPI: @timed_dbg`. Forms:

```julia
@timed_dbg some_call(...)              # label taken from the callee
@timed_dbg "label" begin ... end       # explicit label
@timed_dbg level=2 inner_call(...)     # records only when the active level >= 2
```

`level` and `label` are independent; the default level is 1.

### Zero cost when disabled

The macro expands to `if DataBrowserAPI.profile_level() >= level; <timed>; else <expr> end`.
`profile_level()` returns `0` by default, so the guard is statically false and the compiler
removes the timed branch entirely — no timer, no allocation, no call. When the profiler is
not loaded, `@timed_dbg` is exactly the wrapped expression, so production carries no timing
overhead and no dependency on TimerOutputs.

### Levels

Levels are an overhead dial, not an organizing tool. Leave ordinary markers at the default
(level 1). Reserve `level >= 2` for chatty inner-loop markers you do not want recorded in a
default run; they light up only when the active level is raised. Use module granularity and
TimerOutputs' display options (`maxdepth`, subtree indexing) to focus — not extra levels.

## Enabling instrumentation

The dev-only `DataBrowserProfiling` package drives everything. Loading it, from an
environment that depends on it, installs the timing hooks into `DataBrowserAPI`, raises
`profile_level()`, and starts the collector:

```julia
using DataBrowserProfiling   # turns @timed_dbg on, process-wide, for the session
```

That redefinition recompiles the annotated methods with their timed branch live. It runs at
load time (never during precompilation) and survives Revise reloads, so markers added or
edited mid-session stay instrumented.

## Collector model

The collector is task-safe by construction, which matters because engine work runs on a fixed
pool of long-lived worker tasks that migrate across Julia threads:

- Each Julia task records into its **own** `TimerOutput`, keyed on the task (never on
  `threadid()`), so no two tasks touch the same timer and recording needs no lock.
- Nested `@timed_dbg` sections on a task build a native parent/child tree in that task's timer.
- When a task's **outermost** section exits, its completed timer is merged into a shared master
  under a lock, then cleared. The macro's `finally` runs this even when the timed call throws.
- All access to the master — merges, snapshots, resets — is serialized by that one lock, so a
  snapshot never races a merge.

Because concurrent tasks accumulate work that overlaps in wall time, summed row times can exceed
the elapsed wall time. TimerOutputs measures accumulated work, not a thread timeline.

## Reading results

Every entry point returns a native `TimerOutputs.TimerOutput`:

```julia
snapshot_debug_timings()   # independent copy of everything accumulated so far (peek)
take_debug_timings!()      # snapshot AND reset the master, to isolate the next interval
finish_debug_timings!()    # stop recording and return the final tree
reset_debug_timings!()     # discard accumulated timings and (re)enable recording
```

A result contains only work whose outermost section has already been submitted; sections still
running on a worker appear in a later result. Display and analysis are delegated to TimerOutputs
directly — `show`, `print_timer`, `flatten`, subtree indexing (`t["scan_source"]`), and the
Tables.jl interface (`DataFrame(t)`, `CSV.write(path, t)`). This layer adds no custom report
format.

## Sampling profiler and RSS

`DataBrowserProfiling` also ships a bounded Julia CPU sampling profiler
(`start_sampling!` / `stop_sampling!` / `cancel_sampling!`) and a process-RSS helper
(`process_rss_bytes`). These are separate dev-only tools that share no state with the
instrumentation collector above, which is the primary path.
