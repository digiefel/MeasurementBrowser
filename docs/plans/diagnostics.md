# Diagnostics

## Purpose

Diagnostics are the tools that explain why a workspace or project callback is slow, failing, or using
too much memory. They are not a general telemetry system for the browser itself.

The Performance window is the control surface. `DataBrowserProfiling` should become the small package
that captures and summarizes diagnostic runs for project code: source detection, source reading,
entry creation, item processing, item analysis, plot setup, and plot drawing.

## Scope

Keep three levels separate:

1. **Workspace health**: current scan/cache/work status, visible failures, queue counts, cache buffer
   counts, and process memory. These values are cheap, live, and belong to the workspace or GUI state.
2. **Callback timings**: elapsed time and allocations for project callbacks and selected plot work.
   These are bounded summaries retained for the Performance window.
3. **Explicit diagnostic runs**: user-triggered captures for one selected path, using Julia's standard
   tools where possible.

Everything else is out of scope unless a real debugging need forces it back in.

## What To Remove

Remove the always-available internal tracing product:

- Perfetto JSON export.
- The broad internal span taxonomy for engine implementation details.
- Counter tracks sampled for trace export.
- Crash-forensics JSONL tracing as a normal workspace option.
- A separate engine-profiler UI.
- Makie-based performance figures in the GUI.

The browser can still log errors and expose current workspace status. It should not retain a second
structured event stream for normal use.

## What To Keep

Keep the parts that help diagnose user/project code:

- Per-source-item scan rows: `detect`, `read`, `entries`, item count, total time, and thread ids.
- Per-registration and concrete-type scan summaries.
- Plot redraw phase rows: materialize/load, setup, draw, and total, with allocations.
- A bounded history of recent timing/allocation samples.
- Julia CPU sampling for an explicit user-triggered diagnostic run.
- A reduced hotspot table: self time, total time, function, file, and line.
- Text or CSV export of diagnostic summaries when useful.

Keep process memory and workspace counters, but treat them as app health, not profiling traces.

## Performance Window

The Performance window should have four stable sections:

| Section | Shows | Source |
|---|---|---|
| Workspace | status, failures, queue/cache counts, memory | workspace snapshot and cheap process queries |
| Callbacks | scan rows, plot phases, allocations | callback timing summaries |
| Diagnostic Run | controls for profiling the selected callback path | explicit user action |
| Results | hotspot table, timing table, simple visual summary | last diagnostic result |

The window should not require GLMakie. Use CImGui tables and lightweight visuals:

- progress bars for fractions and relative costs;
- cell bars inside timing tables;
- sparklines for bounded histories;
- flamegraph-style rectangles for sampled hotspots if the table is not enough;
- event lists only for user-triggered diagnostic results.

If these visuals need axes, zooming, or dense interaction later, add ImPlot then. Do not add ImPlot
just to replace a few sparklines and bars.

## Diagnostic Runs

A diagnostic run is explicit and scoped. It answers one question:

- Why is this selected plot slow?
- Why is processing this item slow?
- Why is reading this source item slow?
- Where are allocations coming from in this callback?

The run should reuse the same callback path the browser uses. It should not introduce a parallel debug
execution path with different behavior.

Initial run modes:

- `timed`: elapsed time and allocations around the selected callback path.
- `sampled`: Julia `Profile` sampling around the selected callback path.

Possible later modes, only if needed:

- allocation profiling when Julia tooling makes it practical for the supported Julia version;
- `JET` checks for selected project code;
- handoff helpers for `Cthulhu` or `Debugger.jl`.

Do not implement a debugger inside DataBrowser. If debugging needs grow, integrate with the Julia
debugging tools rather than reproducing them.

## Package Boundary

`DataBrowserProfiling` should own diagnostic data types and capture helpers:

- timing sample rows;
- diagnostic run requests;
- diagnostic run results;
- reduced CPU sample rows;
- formatting-neutral summaries.

It should not own browser windows, CImGui rendering, Makie figures, workspace scheduling, cache writes,
or trace-file formats for engine internals.

`DataBrowserGUI` owns the Performance window and visual presentation. It may depend on
`DataBrowserProfiling` for result types and capture helpers. It should not depend on GLMakie for
performance diagnostics.

`DataBrowserCore` owns workspace health and callback execution. It should call small profiling helpers
at callback boundaries, but the work graph should not be shaped around profiling.

## Migration

1. Keep CImGui tables, bars, and sparklines for normal workspace diagnostics.
2. Keep `DebugTimings` as the explicit benchmark diagnostic surface.
3. Keep `DataBrowserProfiling` focused on timing summaries, CPU sampling, and memory helpers.
4. Extend explicit timing call sites only where a concrete diagnostic question needs them.

The first useful version is small: one button profiles the current plot redraw, one table shows the
hotspots, and the GUI no longer imports GLMakie for the Performance window.
