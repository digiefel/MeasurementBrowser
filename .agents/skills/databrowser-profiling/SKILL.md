---
name: databrowser-profiling
description: Interactively profile DataBrowser to find where time goes (startup, a scan, interpret/process/cache work, a GUI action, a plot). Use whenever investigating DataBrowser performance in this repo — a slow operation, a hotspot, throughput or latency, or the timing impact of a change. Runs a live Revise session against the RuO2 project and reads TimerOutputs timing trees. Instrumentation, not a sampling profiler.
---

# Profiling DataBrowser

`DataBrowserProfiling` is an instrumentation profiler: `@timed_dbg` markers in the
engine are aggregated into native `TimerOutputs.TimerOutput` trees. It is inert
and zero-cost until the package is loaded, task-safe across the async worker pool,
and every result is a plain `TimerOutput` you inspect with TimerOutputs.

You measure by driving a real workload and reading the tree — not by guessing from
the code. This is not a sampling profiler.

## 1. Start a live session

Profile against the copy of the RuO2 project bundled with this skill, under
`project/` next to this file. It defines the real measurement recipes but is a
self-contained copy, so profiling never touches your actual project or data.

**Use the `bench` environment.** It is the only environment that has
`DataBrowserProfiling` together with every package the bundled project needs
(CSV, DataFrames, GLMakie, SmoothData, Revise), all dev'd against this repo. A
user's own analysis environment will not have `DataBrowserProfiling`, and
`using DataBrowserProfiling` there fails outright — do not paper over that by
running without instrumentation, and do not add the package to their project.

`jmux -p bench '<code>'` keeps one persistent Julia REPL per environment, which
removes the ~40 s startup and package-load cost from every iteration. Prefer it
over launching a fresh `julia` for each measurement.

**Load `Revise` first, before anything pulls in `DataBrowser`.** Revise cannot
track packages that were already loaded when it starts, and it fails silently:
`Revise.revise()` returns cleanly and you measure the *old* code. If a session is
already loaded without it, `using Revise; using DataBrowserCache, DataBrowserCore;
Revise.track(DataBrowserCache); Revise.track(DataBrowserCore)` recovers it —
track the module by its own name, not as `DataBrowser.DataBrowserCache`.

```julia
using Revise                          # FIRST — see above
using DataBrowserProfiling            # loading this enables @timed_dbg
include(".../databrowser-profiling/project/definitions.jl")   # defines PROJECT
ws = open_workspace(PROJECT, root; metadata_file="device_info.txt")
browser = open_browser(ws)
```

The cache is DuckDB, single-writer: a workspace left open in the persistent
session locks the project cache, and any other process opening the same project
dies with an IO error. Call `close_workspace!(ws)` before running a second
process against the same project.

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

A snapshot contains only work whose outermost `@timed_dbg` section has already
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

`show(t)` pads every row to the width of the longest label, so one marker whose
auto-generated label is a multi-line closure makes the whole tree unreadable in a
terminal. Give such markers an explicit short label (section 5) rather than
fighting the output. Reading a tree while the workload is still running also
yields impossible numbers (children exceeding parents) because the timer mutates
under you — stop the work, or accept that only finished sections are meaningful.

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
reload — no restart. `@timed_dbg` lives in `DataBrowserAPI`:

```julia
@timed_dbg some_call(...)              # label taken from the callee
@timed_dbg "phase" begin ... end       # explicit label
@timed_dbg level=2 inner_call(...)     # level >= 2: opt-in, for chatty inner loops
```

Sections record when their `level` (default 1) is <= the active level. Markers are
inert in production, so leaving them in the code costs nothing.
