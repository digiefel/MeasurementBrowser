# Plotting API Design

## Goal

Plotting should be simple for project authors and powerful in the GUI.

Project code should define plot kinds and two functions:

```julia
setup_plot(
    project::AbstractProject,
    plot_kind::Type{<:PlotKind},
    measurements::Vector{MeasurementInfo},
)::Figure

plot_data!(
    project::AbstractProject,
    plot_kind::Type{<:PlotKind},
    measurements::Vector{MeasurementInfo},
    figure::Figure,
)::Nothing
```

The package owns selection UI, plot composition, cache behavior, data reuse, and repeated calls into
project plotting functions.

Plotting functions are not cache functions. Projects draw measurement data. The package decides
whether that data comes from already loaded data, valid cache, or source files.

## User Flow

The GUI starts from a measurement selection.

```text
select measurements
  -> package discovers loaded PlotKind types
  -> plot window shows all of them
  -> user chooses one in that plot window
  -> setup_plot once
  -> plot_data! one or more times
```

Overlay should feel effortless:

```text
select measurements
  -> choose plot kind in the main plot window
  -> one shared plot
  -> selected measurements drawn together
```

Composition should feel just as effortless:

```text
select measurements
  -> choose plot kind
  -> add more selected data to the same plot
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

The package discovers concrete plot kind types that are already loaded:

```julia
plot_kinds()::Vector{Type{<:PlotKind}}
```

The GUI shows those type objects by type name:

```julia
Vector{Type{<:PlotKind}}
```

Then calls:

```julia
kind::Type{<:PlotKind} = PUNDLoop
fig::Figure = setup_plot(project, kind, measurements)
plot_data!(project, kind, measurements, fig)
```

This keeps plot kinds searchable in code and lets multiple dispatch route implementation.

The UI does not need project-defined labels, descriptions, default choices, compatibility rules, or
minimum selection rules at this stage. Each plot window owns its selected plot kind. The main plot
window follows the browser selection by default. Detached plot windows keep their own measurements
unless Live is enabled in that window. If the selected measurements do not fit that plot kind, the
project plotting method fails normally and the GUI shows the error.

## Project Responsibilities

Projects define concrete plot kind types and implement plotting:

```julia
setup_plot(
    project::AbstractProject,
    plot_kind::Type{<:PlotKind},
    measurements::Vector{MeasurementInfo},
)::Figure

plot_data!(
    project::AbstractProject,
    plot_kind::Type{<:PlotKind},
    measurements::Vector{MeasurementInfo},
    figure::Figure,
)::Nothing
```

`setup_plot` creates the Makie figure, axes, labels, and layout for that plot kind.

`plot_data!` mutates the figure by adding the selected measurements. When it needs direct
measurement data, it calls:

```julia
read_measurement_data(
    project::AbstractProject,
    measurements::Vector{MeasurementInfo},
)::Vector{DataFrame}
```

When it needs processed measurement data, it calls:

```julia
process_measurement_data(
    project::AbstractProject,
    measurements::Vector{MeasurementInfo},
)::Vector{DataFrame}
```

Project plotting code should not read cache storage or know which cache entry supplied the data.

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

The package supplies `raw_data` by asking for the selected measurement's direct data, not by reading
processed plot results. Debug mode may use cached source or measurement data when valid, but it
should not use already processed normal-plot data.

## Package Responsibilities

The package owns:

- selection UI
- menus and labels
- overlay behavior
- faceting and multi-panel layout
- caching source-file and measurement data
- cache storage, freshness, and repair policy
- calling `setup_plot` once and `plot_data!` as many times as needed
- calling `debug_plot` for explicit debug workflows

Package-level composition should not require project-specific branching. The package creates the
figure, calls the selected plot function, and reports normal project errors when a selection does
not fit that plot kind.

## Data Flow

Normal plots use logical measurement data:

```text
MeasurementInfo
  -> setup_plot(...)
  -> plot_data!(...)
  -> read_measurement_data(project, measurements) or process_measurement_data(project, measurements)
```

`read_measurement_data` returns one `DataFrame` per `MeasurementInfo`, in the same order. The package
owns this function and can cache its result.

`process_measurement_data(project, measurements)` returns one processed `DataFrame` per
`MeasurementInfo`, in the same order. The package owns the cached vector form; projects implement
the single-measurement conversion from direct data to processed data.

When the package needs source data, it asks the project for measurement-scoped source data:

```julia
load_source_data(
    project::AbstractProject,
    source_file::SourceFile;
    measurement::MeasurementInfo,
)::DataFrame
```

The project owns any mapping from a physical source file to a logical measurement dataframe. For
example, RuO2 PUND fatigue files contain many cycles in one CSV, and the measurement-scoped
`load_source_data` call returns only the selected cycle. Plotting code should not know how that
selection works.

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
  read_measurement_data
  process_measurement_data

UI/render work:
  setup_plot
  plot_data!
  debug_plot
```

Cancellation belongs to package-run background jobs. It should not appear in project plotting
signatures.

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
