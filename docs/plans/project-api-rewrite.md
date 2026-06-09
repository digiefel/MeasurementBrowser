# Project API Rewrite Sketch

## Goal

Adding a new measurement procedure should feel like writing one project module while the package
provides scanning, caching, GUI state, plots, and figure export around it.

The package owns orchestration, opened-project data, cache storage, cache freshness, cache repair,
and UI data access. A project owns file meaning, measurement discovery, analysis, and
plotting.

Project-specific cache modules such as `projects/RuO2/Cache.jl` should disappear. Project code
produces project results. The package decides whether those results come from cache or source files.

## What A Project Is

A project is the code object that teaches `MeasurementBrowser` how to understand one measurement
family. `RuO2Project` and `TASEProject` are examples. When the app opens a folder, the package
keeps that opened folder's data on the `AbstractProject` value: root path, cache identity, indexes,
and loaded data. Project scripts receive the project value, but they should not manage those fields
or know how they are stored.

The project API is the set of functions that package code and project code use to talk to each
other. Package code calls project functions when it needs project-specific meaning. Project code
calls package functions when it needs shared services such as cached data. Project authors should
not need to know how cache storage, cache repair, tree data, or GUI data are implemented.

## Source File Functions

Source files are part of the project API because they are the unit the scanner discovers and the
unit the cache repairs. The package can find files and fingerprint them, but only the project knows
what measurements a file contains and how to load its data.

The source-file functions have two different outputs: browser measurements and source data.

```julia
index_source_file(
    project::AbstractProject,
    source_file::SourceFile,
)::Vector{MeasurementInfo}
```

`index_source_file` returns the logical browser measurements found in one physical source file. It
may inspect file contents when filenames and cheap metadata are not enough. RuO2 fatigue is the
current example: the project must inspect the file to know which cycles exist.

```julia
load_source_data(
    project::AbstractProject,
    source_file::SourceFile;
    measurement::Union{Nothing,MeasurementInfo}=nothing,
)::DataFrame
```

`load_source_data(project, source_file)` returns reusable table data for one physical source file.
Cache building, plot loading, stats, and exports can use it. It is not a mandatory predecessor to
`index_source_file`; it is the named way to get source data when data is needed.

`load_source_data(project, source_file; measurement=m)` returns the table data for one logical
browser measurement from that source file. For one-file-one-measurement sources this may be the same
data as the full source file. For virtual measurements, the project owns the
measurement-scoped selection or reshaping. RuO2 PUND fatigue is the current example: the full file
contains many cycles, and the measurement-scoped call returns the selected cycle.

```julia
process_measurement_data(
    project::AbstractProject,
    measurement::MeasurementInfo,
    data::DataFrame,
)::DataFrame
```

`process_measurement_data(project, measurement, data)` returns the processed dataframe for one
logical measurement. The package calls it lazily from
`process_measurement_data(project, measurements)` when processed data is requested and no valid
processed cache entry exists. The default implementation returns `data` unchanged.

The current package helper named `index_source_file(filepath)` creates a cheap `SourceFile` record
from a path. Implementation must rename that helper before using
`index_source_file(project, source_file)` as the project function.

## Project Load Flow

```text
PACKAGE
  collect_source_files(root::AbstractString)::Vector{SourceFile}

PACKAGE -> PROJECT
  index_source_file(
      project::AbstractProject,
      source_file::SourceFile,
  )::Vector{MeasurementInfo}

PACKAGE -> PROJECT
  compute_and_add_measurement_stats!(
      project::AbstractProject,
      measurements::Vector{MeasurementInfo},
      source_files::Vector{SourceFile},
  )::Nothing
```

This is the scan/index path. It creates the browser tree and the source-file list used for cache
comparison.

## Data Flow

```text
PACKAGE -> PROJECT, UI THREAD
  setup_plot(
      project::AbstractProject,
      plot_kind::Type{<:PlotKind},
      measurements::Vector{MeasurementInfo},
  )::Figure

PACKAGE -> PROJECT, UI THREAD
  plot_data!(
      project::AbstractProject,
      plot_kind::Type{<:PlotKind},
      measurements::Vector{MeasurementInfo},
      figure::Figure,
  )::Nothing

PROJECT -> PACKAGE, inside plot_data! or other data-consuming project code
  read_measurement_data(
      project::AbstractProject,
      measurements::Vector{MeasurementInfo},
  )::Vector{DataFrame}

PROJECT -> PACKAGE, when processed measurement data is needed
  process_measurement_data(
      project::AbstractProject,
      measurements::Vector{MeasurementInfo},
  )::Vector{DataFrame}

PACKAGE -> PROJECT, when source data is needed
  load_source_data(
      project::AbstractProject,
      source_file::SourceFile;
      measurement::Union{Nothing,MeasurementInfo}=nothing,
  )::DataFrame
```

This is the data path. The package owns `read_measurement_data` and the cached vector form of
`process_measurement_data`, but project plotting code calls them when it needs data for
already-indexed measurements.

## Data Access

`read_measurement_data` is package-owned. Plotting, stats, exports, and scripts call it when they
need loaded direct data for logical measurements.

`process_measurement_data(project, measurements)` is also package-owned. It returns cached processed
data when available. On a miss, it reads direct data and calls the project method
`process_measurement_data(project, measurement, data)` for each missing measurement.

The intended order is:

```text
already loaded data
  -> valid cache
  -> source file load
```

On a cache miss, the package calls `load_source_data(project, source_file; measurement=m)` to read
data for a logical measurement. This keeps project-specific source-to-measurement mapping in project
code without adding a separate selector function. Repeated source reads are allowed when the project
needs different representations or the cache does not contain the requested data. The package should
make those reads explicit and avoid them when cached data is valid.

## Plotting

Plotting is part of the project API, but the plotting contract lives in
[plotting-api-design.md](plotting-api-design.md). This document only owns the project boundary:
project code defines plot kinds and draws data; the package owns selection, plot windows, caching,
composition, and repeated calls into project plotting code.

## RuO2 Fatigue Example

```julia
function index_source_file(
    project::RuO2Project,
    source_file::SourceFile,
)::Vector{MeasurementInfo}
    kind::Symbol = detect_kind(project, source_file.filename)
    kind === :unknown && return MeasurementInfo[]

    base::MeasurementInfo = measurement_from_filename_and_header(project, source_file, kind)

    if kind === :pund_fatigue
        source_data::DataFrame = load_source_data(project, source_file)
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
    measurement::Union{Nothing,MeasurementInfo}=nothing,
)::DataFrame
    if detect_kind(project, source_file.filename) === :pund_fatigue
        source_data = read_pund_fatigue_file(source_file.filepath)
        measurement === nothing && return source_data
        cycle = Int(measurement.parameters[:fatigue_idx])
        return select_pund_fatigue_cycle(source_data, cycle)
    end
    measurement !== nothing && return read_standard_ruo2_measurement(
        project,
        source_file,
        measurement,
    )
    return read_standard_ruo2_file(project, source_file)
end

function plot_data!(
    project::RuO2Project,
    plot_kind::Type{<:PlotKind},
    measurements::Vector{MeasurementInfo},
    figure::Figure,
)::Nothing
    dfs::Vector{DataFrame} = read_measurement_data(project, measurements)
    for (measurement, df) in zip(measurements, dfs)
        # add df to figure
    end
    return nothing
end
```

This keeps the read reasons explicit. `index_source_file` reads full fatigue data because cycle
discovery requires it. Later plot calls use `read_measurement_data`, which should reuse valid cache
data before calling `load_source_data(project, source_file; measurement=m)` for the selected logical
measurement.

## Thread Rules

These can run in background/cache work:

```text
index_source_file
load_source_data
read_measurement_data
compute_and_add_measurement_stats!
```

These mutate Makie objects and should run on the UI/render path:

```text
setup_plot
plot_data!
```

Cache data and selected measurement data if useful. Do not cache or mutate Makie figures from
worker threads.

Cancellation is package job control, not project API. Background workers may stop between files,
measurements, or cache entries without adding cancellation arguments to project functions.

## Migration Target

```text
rename the current low-level index_source_file(filepath) helper
rename interpret_file(project, source_file) to index_source_file(project, source_file)
add load_source_data(project, source_file; measurement=nothing)
add read_measurement_data(project, measurements)
make RuO2 fatigue use load_source_data where it currently reads fatigue CSVs directly
route plot data access through read_measurement_data
replace file-oriented plotting functions with setup_plot / plot_data!
remove project-specific cache modules
```
