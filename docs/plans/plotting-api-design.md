# Plotting API Design

## Goal

Plotting should be simple for project authors and powerful in the GUI.

Project code should define plot kinds and two hooks:

```julia
setup_plot(
    project::AbstractProject,
    plot_kind::Type{<:PlotKind},
    measurements::Vector{MeasurementInfo},
)::Figure

plot_data(
    project::AbstractProject,
    plot_kind::Type{<:PlotKind},
    measurements::Vector{MeasurementInfo},
    figure::Figure,
)::Nothing
```

The package owns selection UI, plot composition, cache storage, data reuse, and repeated calls into
project plotting hooks.

Plotting hooks are not cache hooks. Projects draw measurement data. The package decides whether
that data comes from memory, preload, HDF5, or disk.

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
kind::Type{<:PlotKind} = PUNDLoop
fig::Figure = setup_plot(project, kind, measurements)
plot_data(project, kind, measurements, fig)
```

This keeps plot kinds searchable in code and lets multiple dispatch route implementation.

For now, the GUI does not filter the list by selection. If the user chooses a plot kind that does
not support the selected measurements, `setup_plot` or `plot_data` should fail normally.

## Project Responsibilities

Projects define:

```julia
setup_plot(
    project::AbstractProject,
    plot_kind::Type{<:PlotKind},
    measurements::Vector{MeasurementInfo},
)::Figure

plot_data(
    project::AbstractProject,
    plot_kind::Type{<:PlotKind},
    measurements::Vector{MeasurementInfo},
    figure::Figure,
)::Nothing
```

`setup_plot` creates the Makie figure, axes, labels, and layout for that plot kind.

`plot_data` mutates the figure by adding the selected measurements. When it needs loaded
measurement data, it calls:

```julia
data_of_measurements(
    project::AbstractProject,
    measurements::Vector{MeasurementInfo};
    should_cancel::Union{Nothing,Function}=nothing,
)::Vector{DataFrame}
```

Project plotting code should not read HDF5 cache storage or know which cache entry supplied the
data.

Projects may also define an interactive debug plot for cases where the user needs to tune analysis
parameters rather than only view the finished result:

```julia
debug_plot(
    project::AbstractProject,
    measurement::MeasurementInfo,
    raw_data::DataFrame,
)::Figure
```

`debug_plot` is for project-specific diagnostic tools. For example, RuO2 PUND debugging needs raw
waveform data, sliders for pulse-detection parameters, and plots that update when those sliders
change.

## Package Responsibilities

The package owns:

- selection UI
- menus and labels
- overlay behavior
- faceting and multi-panel layout
- caching source-file and measurement data
- cache storage, freshness, and repair policy
- calling `setup_plot` once and `plot_data` as many times as needed
- calling `debug_plot` for explicit debug workflows

Package-level composition should not require project-specific branching.

## Data Flow

Normal plots use logical measurement data:

```text
MeasurementInfo
  -> data_of_measurements(project, measurements)
  -> setup_plot(...)
  -> plot_data(...)
```

`data_of_measurements` returns one `DataFrame` per `MeasurementInfo`, in the same order. The package
owns this accessor and can cache its result.

Interactive debug plots use raw data because changing a slider can change the analysis itself:

```text
MeasurementInfo
  -> raw source/measurement dataframe
  -> debug_plot(project, measurement, raw_data)
  -> Makie Observables update the figure in place
```

This keeps the two workflows distinct:

```text
normal plot:
  show the result
  static figure
  cache-friendly

debug plot:
  explain and tune the analysis
  live figure with sliders
  recomputes inside the figure
```

Makie mutation stays on the UI thread:

```text
background/cache work:
  load_source_data
  data_of_measurements

UI/render work:
  setup_plot
  plot_data
  debug_plot
```

## Interactive Debugging

Debug plots are not a replacement for normal plotting. They are explicit tools for project-specific
analysis debugging.

The package should not know about RuO2 PUND smoothing windows, pulse thresholds, or similar
parameters. It only loads the relevant raw data and embeds the returned Makie figure.

Inside the debug figure, the project can use Makie `Observable`s:

```text
slider changes
  -> Observable updates
  -> project recomputes analysis from raw data
  -> lines, markers, boundaries, and labels update in place
```

For RuO2 PUND, the first useful debug figure should show:

- raw voltage/current samples and smoothed traces
- derivative signal, threshold, ramp regions, and pulse boundaries
- I-V loop
- Q-V or P-V loop

The normal plot cache should not store a debug figure as an analyzed plot payload. Debug figures are
interactive workspaces built from raw data.
