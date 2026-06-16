# Measurement Workspace Vision

## Goal

MeasurementBrowser should become a live, scriptable measurement workspace. The GUI and the Julia API
should expose the same operations: open a project, select measurements, view tables, plot columns,
compose figures, annotate results, save the workflow, reload it later, and keep working while source
files or project code change.

The target experience is close to working in Makie from a Julia REPL, but with measurement-aware
state:

```text
open a browser window
select measurements or define a matching rule
open a table, plot, or composed figure
add transformations, fits, and annotations
save the workflow
continue using the REPL while the browser stays live
```

This doc is the north-star product model. Near-term migration work lives in
[refactor-gap-report.md](refactor-gap-report.md); focused plans linked from [README.md](README.md)
cover only areas that need extra design detail.

## Principles

- **GUI/API parity.** A useful GUI action should map to a package operation that can also be called
  from Julia code. The GUI is an interface to the workspace model, not a separate implementation.
- **Non-blocking REPL.** Starting the browser should not consume the Julia session. A user should be
  able to keep using the same REPL to inspect data, run helper functions, change project code, and
  drive the app.
- **Live workspace.** Source files, project code, selections, and views should be able to update
  without restarting the app. Long work should have visible progress and should not block unrelated
  interaction.
- **Small project boundary.** Project code should stay close to natural scripts: identify
  measurements, load data, compute project-specific results, and define domain-specific views when
  needed. It should not manage cache files, UI state, background jobs, or package lifecycle.
- **Generic tools first.** Built-in visualizers should cover common work: raw table views, X-vs-Y
  plots, line/scatter overlays, heatmaps, histograms, summaries, and fit inspection. Project code
  should extend those tools, not replace the whole browsing and composition system.
- **Persistent work, not exported scripts.** Figure scripts are a transitional export path. The
  long-term artifact is a saved workflow: a file-saveable set of actions that can be replayed,
  edited, inspected, and driven from either GUI or API.
- **Composable figures.** Figures should be editable objects in the app. Some figures are backed by
  live measurement selections; others are self-contained and do not require a source data root.

## Workspace Concepts

The future app should have a small set of named concepts that apply equally to GUI and API use.

| Concept | Meaning |
|---|---|
| Opened project | A source root plus project implementation plus package-managed state: cache identity, indexes, loaded data, jobs, and watchers. |
| Measurement | A logical measurement entry discovered from source files. It has stable identity, parameters, stats, and access to data. |
| Selection | A concrete list of measurements or a rule that resolves to measurements as the project changes. |
| Visualizer | A reusable view that can render data: raw table, scatter, line, heatmap, summary, domain-specific plot, or debug tool. |
| Figure | A composed visual artifact containing visualizer outputs, layout, styling, annotations, and optional embedded data. |
| Annotation | User-authored tags, notes, coordinates, and spatial positions attached to devices or measurements. |
| Figure annotation | A visual mark attached to a figure, such as text, an arrow, a fit label, or a region. |
| Workflow | A persisted sequence or graph of actions: open/select/transform/visualize/fit/annotate/export. The exact file format and DSL are TBD. |

## Scriptable Shape

The exact API is not decided, but it should feel like normal Julia composition rather than an app
control protocol. This sketch is illustrative:

```julia
workspace = open_workspace(RuO2Project(), root)
browser = open_browser(workspace; wait=false)

selection = select_items!(workspace;
    kind=:pund,
    collection=r"RuO2test/A9",
)

fig = Figure(workspace)
plot!(fig, selection, X(:voltage_V), Y(:polarization_uCcm2); visualizer=LinePlot)
fit = linear_fit!(fig, selection, X(:voltage_V), Y(:current_A))
annotate!(fig, Arrow(fit, label="linear region"))

save_workflow("pund_review.mbflow", fig)
```

The GUI should issue the same kinds of actions when the user clicks, drags, selects plot options, or
adds annotations. The workflow file should store those actions in package-owned terms, not in
private widget state.

## Visualizers

Generic visualizers should be modular and composable. The first useful set is:

- raw table view for any loaded tabular measurement data
- X-vs-Y scatter and line plots
- overlays collected by measurement, device, tag, or parameter
- heatmaps for gridded or pivoted tabular data
- simple summaries and histograms
- fit views for common models, starting with linear fits
- project-specific debug visualizers for analysis tuning

Visualizers should consume measurement data through package data-access functions. They should not
read source files, cache files, or GUI state directly.

Project-specific visualizers remain valuable, but they should look like extensions of the same
system: they declare what data they need, draw into a figure or panel, and let the package own
selection, composition, persistence, and background work.

## Figures And Annotations

Figures should support both live and self-contained modes.

```text
live figure
  references an opened project and one or more selections
  updates when matching measurements, data, or project code change

self-contained figure
  stores enough layout, plotted values, style, and annotations to reopen without source files
  remains editable in the app
```

Plot annotations include user-authored visual marks such as:

- linear fits and fit labels
- arrows
- circles or regions of interest
- text labels
- callouts attached to data points or axes

These are separate from device and measurement annotations. A tag or note describes a device or measurement. A plot
arrow describes a figure.

Open design questions:

- Whether a self-contained figure stores a snapshot of plotted data, only rendered marks, or both.
- How live figures degrade when the source project is unavailable.
- Whether plot annotations are stored inside workflow files, figure files, or a shared annotation
  store.

## Workflows

A workflow is the persistent form of a working session. It should capture intent rather than only
pixels:

```text
open project
select measurements by rule
load or process data
create visualizer
map columns to axes
set collection axis and style
add fit
add annotations
export or save figure
```

The workflow may eventually have a DSL, but the DSL is not the first requirement. The first
requirement is a stable action model that both GUI and API can use. A DSL can be added later as a
text representation of that model.

Workflows should support both concrete selections and live rules. Concrete selections reproduce the
exact same measurements. Live rules keep views attached to matching measurements as files are added,
removed, or changed.

## State And Persistence

The app needs clearer state ownership before workflows and a non-blocking REPL can work well. A
useful split is:

| State | Owner | Persistence |
|---|---|---|
| Project state | Opened-project object | Recomputed from source, cache, and source-root metadata. |
| Data cache | Package cache layer | Generated HDF5 cache outside the source root. |
| GUI state | Browser window | Mostly transient; small preferences can live in app prefs. |
| Annotations | Source root | Hand-editable tags, notes, coordinates, and spatial positions. |
| Workflow state | Workflow model | Saved workflow file, replayable from GUI or API. |
| Figure state | Figure/workflow model | Saved with workflow or as a self-contained figure artifact. |

`ui_state` is a current implementation detail for GUI rendering. It should not become the durable
owner of project data, loaded data, workflow actions, or figure contents. The likely direction is a
package-owned opened-workspace object that the GUI observes and mutates through the same operations
available to Julia callers.

Open design questions:

- Whether opened-project state should be a wrapper around project implementations or fields on
  concrete project values.
- How much GUI layout state belongs in app preferences versus workflow files.
- Whether workflow persistence should store an append-only action log, a normalized object graph, or
  both.
- How self-contained figures share the workflow model without requiring a source project.

## Non-Blocking Browser

The browser should be a window attached to a live Julia workspace, not the only thing the Julia
process can do. The intended REPL flow is:

```text
start Julia in the user's preferred IDE
load MeasurementBrowser and project code
open the browser window without blocking the REPL
use the GUI for inspection and the REPL for scripting
change project code and see views refresh when possible
access the same opened-project data from code and GUI
```

This affects architecture:

- opened-project state needs a real owner that is usable outside `ui_state`
- GUI state should reference workspace objects, not contain the whole data model
- long-running work needs job ownership, progress, cancellation, and failure reporting
- Makie/CImGui work must stay on the appropriate UI/render path
- data access must be thread-aware and clear about what can run in background work

## Relationship To Existing Plans

- [refactor-gap-report.md](refactor-gap-report.md) tracks the current implementation gaps that block
  this model: opened-project ownership, data reuse, progressive analysis, plotting composition, and
  stale code paths.
- [spatial-browser.md](spatial-browser.md) covers annotations and spatial navigation;
  plot annotations are a separate future layer.
- [measurement-parameters-and-stats.md](measurement-parameters-and-stats.md) defines the metadata
  bucket semantics that selections, visualizers, workflows, and saved figures should rely on.
- [figure_scripts.md](../figure_scripts.md) documents the deprecated script export path that
  workflows should replace.
