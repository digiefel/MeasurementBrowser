# Project API Rewrite Sketch

## Goal

Adding a new measurement procedure should feel like writing one project module while the package
provides scanning, caching, GUI state, plots, and figure export around it.

The package owns orchestration, opened-project state, cache storage, cache freshness, cache repair,
and UI data access. A project owns source-file meaning, measurement discovery, analysis, and
plotting.

Project-specific cache modules such as `projects/RuO2/Cache.jl` should disappear. Project hooks
produce project results. The package decides whether those results come from memory, HDF5, or disk.

## Naming Decision

The source-file API has two different outputs.

```julia
index_source_file(
    project::AbstractProject,
    source_file::SourceFile;
    should_cancel::Union{Nothing,Function}=nothing,
)::Vector{MeasurementInfo}
```

`index_source_file` returns the logical browser measurements found in one physical source file. It
may inspect file contents when filenames and cheap metadata are not enough. RuO2 fatigue is the
current example: the project must inspect the file to know which cycles exist.

```julia
load_source_data(
    project::AbstractProject,
    source_file::SourceFile;
    should_cancel::Union{Nothing,Function}=nothing,
)::DataFrame
```

`load_source_data` returns reusable table data for one physical source file. Cache building, plot
loading, stats, and exports can use it. It is not a mandatory predecessor to `index_source_file`;
it is the named way to get source data when data is needed.

The current code already has a package helper named `index_source_file(filepath)`, which creates a
cheap `SourceFile` record from a path. Implementation must rename that helper before using
`index_source_file(project, source_file)` as the project hook.

## Core Flow

```text
PACKAGE
  collect_source_files(root::AbstractString)::Vector{SourceFile}

PACKAGE -> PROJECT
  index_source_file(
      project::AbstractProject,
      source_file::SourceFile;
      should_cancel::Union{Nothing,Function}=nothing,
  )::Vector{MeasurementInfo}

PACKAGE -> PROJECT
  compute_and_add_measurement_stats!(
      project::AbstractProject,
      measurements::Vector{MeasurementInfo},
      source_files::Vector{SourceFile},
  )::Nothing

PACKAGE
  data_of_measurements(
      project::AbstractProject,
      measurements::Vector{MeasurementInfo};
      should_cancel::Union{Nothing,Function}=nothing,
  )::Vector{DataFrame}

PACKAGE -> PROJECT, when source data is needed
  load_source_data(
      project::AbstractProject,
      source_file::SourceFile;
      should_cancel::Union{Nothing,Function}=nothing,
  )::DataFrame

PACKAGE -> PROJECT, UI THREAD
  setup_plot(
      project::AbstractProject,
      plot_kind::Type{<:PlotKind},
      measurements::Vector{MeasurementInfo},
  )::Figure

PACKAGE -> PROJECT, UI THREAD
  plot_data(
      project::AbstractProject,
      plot_kind::Type{<:PlotKind},
      measurements::Vector{MeasurementInfo},
      figure::Figure,
  )::Nothing
```

## Data Access

`data_of_measurements` is package-owned. Plotting, stats, exports, and scripts call it when they
need loaded data for logical measurements.

The intended order is:

```text
memory
  -> valid HDF5 cache
  -> source-data load from disk
```

Disk fallback may call `load_source_data`. Repeated source reads are allowed when the project needs
different representations or the cache does not contain the requested data. The package should make
those reads explicit and avoid them when cached data is valid.

Do not add a separate public project hook for selecting measurement rows from source data until the
implementation need is clear. The current code already has project-specific plot loaders; the first
implementation can route through those while the cache boundary is simplified.

## Plotting

Project hooks:

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

`setup_plot` creates the figure, axes, labels, and layout.

`plot_data` mutates `figure`. When it needs loaded measurement data, it calls:

```julia
data_of_measurements(
    project::AbstractProject,
    measurements::Vector{MeasurementInfo};
    should_cancel::Union{Nothing,Function}=nothing,
)::Vector{DataFrame}
```

Project plotting code should not read HDF5 cache storage directly.

## RuO2 Fatigue Example

```julia
function index_source_file(
    project::RuO2Project,
    source_file::SourceFile;
    should_cancel::Union{Nothing,Function}=nothing,
)::Vector{MeasurementInfo}
    kind::Symbol = detect_kind(project, source_file.filename)
    kind === :unknown && return MeasurementInfo[]

    base::MeasurementInfo = measurement_from_filename_and_header(project, source_file, kind)

    if kind === :pund_fatigue
        source_data::DataFrame = load_source_data(project, source_file; should_cancel)
        cycles::Vector = unique(source_data.cycle)
        return [
            MeasurementInfo(base;
                unique_id="$(base.filepath)#fatigue_count=$(Int(cycle))",
                measurement_kind=:pund,
                parameters=merge(
                    deepcopy(base.parameters),
                    Dict{Symbol,Any}(:fatigue_idx => cycle),
                ),
            )
            for cycle in cycles
        ]
    end

    return [base]
end

function load_source_data(
    project::RuO2Project,
    source_file::SourceFile;
    should_cancel::Union{Nothing,Function}=nothing,
)::DataFrame
    if detect_kind(project, source_file.filename) === :pund_fatigue
        return read_pund_fatigue_file(source_file.filepath; should_cancel)
    end
    return read_standard_ruo2_file(project, source_file; should_cancel)
end

function plot_data(
    project::RuO2Project,
    plot_kind::Type{<:PlotKind},
    measurements::Vector{MeasurementInfo},
    figure::Figure,
)::Nothing
    dfs::Vector{DataFrame} = data_of_measurements(project, measurements)
    for (measurement, df) in zip(measurements, dfs)
        # add df to figure
    end
    return nothing
end
```

This keeps the read reasons explicit. `index_source_file` reads fatigue data because cycle discovery
requires it. Later plot calls use `data_of_measurements`, which should reuse valid cache data before
reading the source file again.

## Thread Boundary

These can run in background/cache work:

```text
index_source_file
load_source_data
data_of_measurements
compute_and_add_measurement_stats!
```

These mutate Makie objects and should run on the UI/render path:

```text
setup_plot
plot_data
```

Cache data and selected measurement data if useful. Do not cache or mutate Makie figures from
worker threads.

## Migration Target

```text
rename the current low-level index_source_file(filepath) helper
rename interpret_file(project, source_file) to index_source_file(project, source_file)
add load_source_data(project, source_file)
add data_of_measurements(project, measurements)
make RuO2 fatigue use load_source_data where it currently reads fatigue CSVs directly
route plot data access through data_of_measurements
replace file-oriented plotting hooks with setup_plot / plot_data
remove project-specific cache modules
```
