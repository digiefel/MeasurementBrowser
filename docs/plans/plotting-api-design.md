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
register_plot!(project, :pund; label = "PUND loop",  setup = …, draw = …)
register_plot!(project, :pund; label = "Q–V",        setup = …, draw = …)   # same kind, another plot
```

- `setup(workspace, items)` builds and returns the `Figure`; `draw(workspace, items, figure)` fills
  it in. `items::Vector{<:DataItem}` are the loaded, data-bearing items for the current selection;
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
| overlay, faceting, multi-panel layout | optional `debug_plot` for analysis tuning |
| data access + cache; materializing `item.data` | the figure/axes/layout a plot kind needs |
| calling `setup` once and `draw` as many times as needed | drawing `item.data` into the figure |
| dispatching `debug_plot` for explicit debug workflows | |

## Normal plots vs. debug plots

Two distinct workflows, kept separate:

```text
normal plot                         debug plot
  show the finished result            explain and tune the analysis
  uses processed item.data            uses raw data
  static figure                       live figure with sliders/Observables
  cache-friendly                      recomputes inside the figure; never cached
```

A **debug plot** is a project-specific diagnostic tool for tuning analysis parameters, not a
replacement for normal plotting:

```julia
debug_plot(workspace, items; kwargs...)::Figure   # stub lives in Visualization.jl today
```

The package only loads the relevant raw data and embeds the returned Makie figure — it knows nothing
about a project's smoothing windows, pulse thresholds, or similar. Inside the figure the project
drives Makie `Observable`s:

```text
slider changes -> Observable updates -> project recomputes analysis from raw data
              -> lines, markers, boundaries, labels update in place
```

The normal-plot cache must never store a debug figure as a plot payload — debug figures are
interactive workspaces built from raw data.

*Illustrative (project-specific, not package): a first RuO2 PUND debug figure would show raw
voltage/current and smoothed traces, the derivative + threshold + ramp regions + pulse boundaries,
the I–V loop, and the Q–V / P–V loop.*

## Thread discipline

Makie mutation stays on the UI/render thread; data resolution can run in the background:

```text
background / cache work:   data access (read_data, process)  → materialize item.data
UI / render work:          setup, draw, debug_plot            → Makie figure construction & mutation
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

This doc owns the **project-facing plot API** (`register_plot!`, `setup`/`draw`, `debug_plot`), the
**normal-vs-debug split**, and the **project↔package composition contract**. It does not own:

- the generic built-in visualizers, figures, annotations, and workflows — [workspace-vision.md](workspace-vision.md);
- the item/data model (`DataItem`, `item.data`, `kind`) the callbacks consume —
  [data-model-generalization.md](data-model-generalization.md).
