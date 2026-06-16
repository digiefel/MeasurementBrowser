# Plotting & Visualizer Design

Plotting should be simple for project authors and powerful in the GUI: a project declares *what* a
plot is, and the package owns *how* plots are composed (overlay, faceting, panels). This doc owns the
project-facing plot API and the normal-vs-debug split; the generic built-in visualizers, figures, and
workflows are [workspace-vision.md](workspace-vision.md). It builds on the item/data model in
[data-model-generalization.md](data-model-generalization.md).

## How a project declares a plot

A project registers plots per measurement `kind` with plain callbacks — no types to subtype, no
reflection. A kind may have **several** plots, each distinguished by its `label`:

```julia
register_plot!(project, :pund; label = "PUND loop", setup = …, draw = …,
               debug = (workspace, items) -> …)     # optional: live raw-data tuning variant
register_plot!(project, :pund; label = "Q–V",       setup = …, draw = …)    # same kind, another plot
```

- `setup(workspace, items)` builds and returns the `Figure`; `draw(workspace, items, figure)` fills
  it in. `items::Vector{<:AbstractDataItem}` are the loaded, data-bearing items for the current selection;
  each carries its own `item.data` (see [data-model-generalization.md](data-model-generalization.md)).
  Callbacks never load data or zip a parallel array — data access and caching are package-owned.
- Plots live in `project.plots :: Dict{Symbol, Dict{String, PlotRecipe}}` (kind → label → recipe).
  Re-registering the same `(kind, label)` replaces that one plot in place, so REPL iteration stays
  stable. `registered_plot_kinds(project, kind)` lists every plot for a kind, sorted by label.
- The engine represents each registered plot internally as `RegistryPlot{kind, label}` and dispatches
  `setup_plot`/`plot_data!` to the recipe callbacks; choices persist as `"kind::label"`. **`PlotKind`
  is internal**; projects do not subtype it, and the old "discover plot-kind types by reflection" path
  is gone.

> Current state: the implemented `setup`/`draw` callbacks still receive `(workspace, measurements,
> processed)` with a *parallel* processed-data array. Step 2 of the data-model plan migrates them to
> the `(workspace, items, figure)` + `item.data` shape above; this doc describes that target.

## User flow (GUI)

```text
select items
  -> package shows the plots registered for the selection's common kind
     (via registered_plot_kinds; mixed-kind selection -> none)
  -> user chooses one (persisted choices are validated, unknown ones ignored)
  -> setup_plot once
  -> plot_data once, or many times for overlay/compose
```

Overlay and composition should feel effortless, and stay package-driven:

```text
select compatible items -> Overlay        -> one shared figure, each item drawn into it
select items -> choose plot -> add more    -> optionally facet or split into panels
```

The project defines plot kinds. **The package decides how to compose them** — no project-specific
branching in the composition layer.

## Responsibilities

| Package owns | Project owns |
|---|---|
| selection UI, menus, labels | `register_plot!` `setup`/`draw` per kind |
| overlay, faceting, multi-panel layout | an optional `debug` callback for analysis tuning |
| data access + cache; materializing `item.data` | the figure/axes/layout a plot kind needs |
| calling `setup` once and `draw` as many times as needed | drawing `item.data` into the figure |
| running `recipe.debug` when debug mode is toggled | |

## Normal plots vs. debug plots

A plot may carry an optional **`debug`** callback — the interactive, raw-data variant of *that* plot,
for tuning the analysis rather than viewing the finished result. It is a fourth keyword on the same
registration, not a separate hook (defaults to `nothing`):

```text
draw (normal)                       debug
  show the finished result            explain and tune the analysis
  processed item.data                 RAW data
  static figure                       live figure with sliders / Observables
  cache-friendly                      recomputes inside the figure; never cached
```

- The GUI's debug toggle is **enabled iff the selected plot's recipe has a `debug` callback**;
  toggling runs `debug` instead of `setup`/`draw`. There is no standalone `debug_plot` hook — the
  engine looks up `recipe.debug` exactly as it looks up `recipe.setup`/`recipe.draw`.
- **`debug` receives raw data, deliberately.** `setup`/`draw` get processed `item.data`; `debug`
  needs raw, because moving a slider re-runs the analysis itself. So its `items` expose raw data and
  the callback owns its own Makie `Observable`s and recompute. That asymmetry is the point of the split.
- `debug` is per `(kind, label)`, so each plot can have its own. A debug tool usually tunes the
  *kind's* analysis, so attaching one to the natural plot for that kind is enough; a debug-only view
  with no normal plot is just a plot whose `draw` is the debug view — not solved until it's needed.

The package knows nothing about a project's smoothing windows, pulse thresholds, or similar — it only
loads the raw data and embeds the returned Makie figure:

```text
slider changes -> Observable updates -> project recomputes analysis from raw data
              -> lines, markers, boundaries, labels update in place
```

The normal-plot cache must never store a debug figure as a plot payload — debug figures are
interactive workspaces built from raw data.

*Illustrative (project-specific, not package): a first RuO2 PUND debug figure would show raw
voltage/current and smoothed traces, the derivative + threshold + ramp regions + pulse boundaries,
the I–V loop, and the Q–V / P–V loop.*

> Current state: debug is still a standalone `debug_plot` engine hook toggled globally. Folding it
> into the recipe — `recipe.debug` lookup in the plot bridge + the GUI toggle — is a small pending
> code change.

## Thread discipline

Makie mutation stays on the UI/render thread; data resolution can run in the background:

```text
background / cache work:   data access (read_data, process)  → materialize item.data
UI / render work:          setup, draw, debug                 → Makie figure construction & mutation
```

## Open questions

- **Cross-kind combined plots.** A plot is offered only when the whole selection shares one `kind`
  (mixed-kind selections show no plots), and `draw` already receives the full `items` vector — so a
  combined plot over a *homogeneous* selection (TLM length sweep, fatigue evolution) is just a
  registered plot for that kind. A plot spanning *heterogeneous* kinds has no home yet; decide whether
  that needs a cross-kind registration or stays out of scope.
- **Relationship to generic visualizers.** Registered plots are project-defined visualizers; the
  generic built-ins (raw table, X–Y, heatmap, histogram, fits) are [workspace-vision.md](workspace-vision.md)'s
  domain. The two should share the same `setup`/`draw` + composition contract.

## Relationship to existing plans

This doc owns the **project-facing plot API** (`register_plot!` with `setup`/`draw`/`debug`), the
**normal-vs-debug split**, and the **project↔package composition contract**. It does not own:

- the generic built-in visualizers, figures, annotations, and workflows — [workspace-vision.md](workspace-vision.md);
- the item/data model (`AbstractDataItem`, `DataItem`, `item.data`, `kind`) the callbacks consume —
  [data-model-generalization.md](data-model-generalization.md).
