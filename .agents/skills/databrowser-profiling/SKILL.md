---
name: databrowser-profiling
description: Interactively profile DataBrowser to find where time goes (startup, a scan, interpret/process/cache work, a GUI action, a plot). Use whenever investigating DataBrowser performance in this repo — a slow operation, a hotspot, throughput or latency, or the timing impact of a change. Runs a live Revise session against the RuO2 project and reads TimerOutputs timing trees. Instrumentation, not a sampling profiler.
---

# Profiling DataBrowser

`DataBrowserProfiling` is an instrumentation profiler: `@time_dbg` markers in the
engine are aggregated into native `TimerOutputs.TimerOutput` trees. It is inert
and zero-cost until the package is loaded, task-safe across the async worker pool,
and every result is a plain `TimerOutput` you inspect with TimerOutputs.

You measure by driving a real workload and reading the tree — not by guessing from
the code. This is not a sampling profiler.

## 1. Start a live session

Profile against the copy of the RuO2 project bundled with this skill, under
`project/` next to this file. It defines the real measurement recipes but is a
self-contained copy, so profiling never touches your actual project or data.

Launch Julia with Revise, load that project's setup (its `define_project` +
`register_*` calls), and open the workspace/browser against a data root. Loading
`DataBrowserProfiling` is the whole switch — it turns timing on process-wide:

```julia
using Revise
using DataBrowserProfiling            # loading this enables @time_dbg
# load the bundled project's definitions (project/), then:
ws = open_workspace(project, root)
browser = open_browser(ws)
```

## 2. Isolate one action and measure it

The pattern is **reset -> do one thing -> take**. The returned `TimerOutput` is
exactly that action's accumulated work:

```julia
take_debug_timings!()                 # zero the interval, discard startup noise
# ... perform ONE action: let a scan finish, select items, draw a plot ...
t = take_debug_timings!()             # t = the cost of just that action
show(t)
```

- `snapshot_debug_timings()` — running total so far, **without** resetting (peek).
- `take_debug_timings!()` — snapshot **and** reset, to isolate the next interval.
- `finish_debug_timings!()` — stop recording and return the final tree at the end.

A snapshot contains only work whose outermost `@time_dbg` section has already
finished; sections still running on a worker show up in a later snapshot.

## 3. Read the tree to find the hotspot

`t` is a `TimerOutput`; use TimerOutputs directly (get its functions with
`using DataBrowserProfiling.TimerOutputs`, or `show(t)` works on its own):

```julia
using DataBrowserProfiling.TimerOutputs
print_timer(t; sortby = :time, maxdepth = 3)   # top of the tree by time
flatten(t)                                       # accumulate same-label sections
t["scan_source"]                                 # drill into a subtree, then further
DataFrame(t)                                      # or CSV.write("t.csv", t) — Tables.jl
```

The TimerOutputs README documents the rest (bars, GC time, allocations, `%par`,
complement rows, flame graphs) — this skill does not restate it.

## 4. Iterate with Revise

Timing stays live across Revise reloads, so the loop is measure -> edit engine code
-> let Revise reload -> measure again and compare:

```julia
before = take_debug_timings!()        # baseline
# ... edit engine code; Revise reloads the changed methods ...
after  = take_debug_timings!()        # compare against `before`
```

## 5. Add resolution where you're blind

If a hot section has no inner detail, drop a marker into the engine and let Revise
reload — no restart. `@time_dbg` lives in `DataBrowserAPI`:

```julia
@time_dbg some_call(...)              # label taken from the callee
@time_dbg "phase" begin ... end       # explicit label
@time_dbg level=2 inner_call(...)     # level >= 2: opt-in, for chatty inner loops
```

Sections record when their `level` (default 1) is <= the active level. Markers are
inert in production, so leaving them in the code costs nothing.
