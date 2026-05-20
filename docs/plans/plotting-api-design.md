# Plotting API Design

## Goal

Plotting should be simple for project authors and powerful in the GUI.

Project code should only need to define:

```julia
setup_plot(project, plot_kind, measurements)
plot_data(project, plot_kind, measurements, fig)
```

The package owns GUI composition: overlay, faceting, panel layouts, and repeated calls into project
plotting hooks.

## User Flow

The GUI starts from a measurement selection.

```text
select measurements
  -> package discovers plot kinds from the project module
  -> show those plot kinds
  -> user chooses one
  -> setup_plot once
  -> plot_data one or more times
```

Overlay should feel effortless:

```text
select compatible measurements
  -> Overlay
  -> one shared plot
  -> selected measurements added to it
```

Composition should feel just as effortless:

```text
select measurements
  -> choose plot kind
  -> add more compatible data to the same plot
  -> optionally facet or split into panels
```

The project defines plot kinds. The package decides how to compose them.

## Plot Kinds

`plot_kind` should be a type, not an unstructured symbol.

```julia
abstract type PlotKind end

struct PUNDLoop <: PlotKind end
struct FatigueSummary <: PlotKind end
struct CVSweepPlot <: PlotKind end
```

The package discovers plot kind types from the project module. A project plot kind is any concrete
type defined in that module that subtypes `PlotKind`.

The GUI receives those type objects:

```julia
Vector{Type{<:PlotKind}}
```

Then calls:

```julia
kind = PUNDLoop
fig = setup_plot(project, kind, measurements)
plot_data(project, kind, measurements, fig)
```

This keeps plot kinds searchable in code and lets multiple dispatch route implementation.

For now, the GUI does not filter the list by selection. If the user chooses a plot kind that does
not support the selected measurements, `setup_plot` or `plot_data` should fail normally.

## Project Responsibilities

Projects define:

```julia
setup_plot(project, ::Type{SomePlotKind}, measurements)::Figure
plot_data(project, ::Type{SomePlotKind}, measurements, fig::Figure)
```

`setup_plot` creates the Makie figure, axes, labels, and layout for that plot kind.

`plot_data` adds the selected measurement data to that figure. It can handle one measurement or a
group, depending on the plot kind.

## Package Responsibilities

The package owns:

- selection UI
- menus and labels
- overlay behavior
- faceting and multi-panel layout
- caching source-file and measurement data
- calling `setup_plot` once and `plot_data` as many times as needed

Package-level composition should not require project-specific branching.

## Data Flow

Plotting uses measurement-level data:

```julia
data_contents(project, measurement)::DataFrame
```

`data_contents` returns the dataframe for one logical measurement. The package can cache it.

Makie mutation stays on the UI thread:

```text
background/cache work:
  file_contents
  data_contents

UI/render work:
  setup_plot
  plot_data
```
