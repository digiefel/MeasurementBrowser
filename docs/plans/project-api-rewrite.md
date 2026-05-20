# Project API Rewrite Sketch

Goal: adding a new measurement procedure should feel like writing one project module while the app
provides scanning, caching, GUI state, plots, and figure export around it.

The package owns orchestration and cache mechanics. A project owns file naming, file meaning,
measurement expansion, stats, analysis, and plotting.

## Core Flow

```text
PACKAGE
  collect_source_files(root)
    -> SourceFile metadata

PACKAGE -> PROJECT
  read_source_file(project, source_file::SourceFile)
    -> DataFrame for the physical file

PACKAGE
  stores the returned contents
  exposes them through file_contents(source_file)::DataFrame

PACKAGE -> PROJECT
  interpret_file(project, source_file::SourceFile)
    -> Vector{MeasurementInfo}

PACKAGE -> PROJECT
  compute_and_add_measurement_stats!(
      project,
      measurements::Vector{MeasurementInfo},
      source_files::Vector{SourceFile},
  )
    -> history and waveform stats

PACKAGE -> PROJECT, UI THREAD
  setup_plot(project, measurements::Vector{MeasurementInfo})
    -> Makie Figure

PACKAGE -> PROJECT, UI THREAD
  plot_data(project, figure::Figure, measurements::Vector{MeasurementInfo})
    -> add cached data for one or more measurements to the figure
```

## Types

Existing package types:

- `SourceFile`: one physical file, cheap metadata, and package-owned access to cached contents.
- `MeasurementInfo`: one logical browser measurement. One source file may produce many
  measurements. It already carries the source filepath.

Target plotting type:

- Makie `Figure`: returned by `setup_plot` and mutated by `plot_data`.

## Source File Contents

Project-defined function to read the file and return its contents for the first time:

```julia
read_source_file(project, file::SourceFile; should_cancel=nothing)::DataFrame
```

The package then later exposes a cached faster read:

```julia
file_contents(file::SourceFile)::DataFrame
```

The package calls `read_source_file`, stores the dataframe with the `SourceFile`, and later serves
it through `file_contents(file)`. This is the whole physical file dataframe. Project code must not
care whether the dataframe came from disk, memory, async preload, or an HDF5 cache.

DataLoader should provide generic CSV helpers only. Projects use those helpers inside
`read_source_file`, then give semantic names to the columns they care about.

`file_contents` may be populated or refreshed off the UI thread. Makie figures are not part of
this cache.
The package should take care to call `read_source_file` in advance. If `file_contents` is called
before `read_source_file`, it should just error out naturally without much fanfare.

## Measurement Data Contents

Project hook:

```julia
data_contents(project, measurement::MeasurementInfo)::DataFrame
```

`data_contents` returns the dataframe for one logical measurement. For a standalone file this may
be the whole source dataframe. For an expanded file, such as fatigue or wakeup, this is the selected
cycle, amplitude, or block. The package can cache these logical-measurement dataframes separately
from the physical source-file dataframe.

Both `file_contents` and `data_contents` should be cheap views over cached data, as cpu-efficient and
memory-efficient as posisble, thread-safe and fast.

## RuO2 Fatigue Example

```julia
function read_source_file(::RuO2Project, file::SourceFile; should_cancel=nothing)
    detect_kind(RUO2_PROJECT, file.filename) === :pund_fatigue || return nothing
    return read_pund_fatigue_file(file.filepath; should_cancel)
end

function interpret_file(::RuO2Project, file::SourceFile; should_cancel=nothing)
    if detect_kind(RUO2_PROJECT, file.filename) === :pund_fatigue
        fatigue_df = file_contents(file)
        cycles = unique(fatigue_df.cycle)
        # create one logical PUND measurement per cycle
    end
end

function compute_and_add_measurement_stats!(::RuO2Project, measurements, files)
    # reuse file_contents(file) instead of rereading the fatigue CSV
end

function setup_plot(::RuO2Project, measurements::Vector{MeasurementInfo})
    # create the figure and axes for this measurement kind
    return fig
end

function plot_data(::RuO2Project, fig::Figure, measurements::Vector{MeasurementInfo})
    for measurement in measurements
        df = data_contents(RUO2_PROJECT, measurement)
        # add df to fig
    end
end
```

For PUND fatigue, `read_source_file` returns the full physical-file dataframe. `interpret_file`
uses it to create logical measurements, and stats/plotting use the same contents to select the
requested cycle. `setup_plot` creates the shared figure/axes once. `plot_data` can be called for
one or many compatible measurements and adds each measurement's data into that same plot.
The selected data comes from `data_contents`, so repeated overlays do not reread the physical file.

The package decides whether `measurements` contains one measurement or a user-selected overlay
group. The project receives the vector either way.

## Thread Boundary

These can run in background/cache work:

```text
read_source_file
file_contents
data_contents
compute_and_add_measurement_stats!
```

These mutate Makie objects and should run on the UI/render path:

```text
setup_plot
plot_data
```

Cache dataframes and selected measurement data if useful. Do not cache or mutate Makie figures
from worker threads.

## Project Hooks To Standardize

Current hooks to keep and make cleaner:

```julia
read_source_file
file_contents
data_contents
interpret_file
compute_and_add_measurement_stats!
setup_plot
plot_data
```

Overlay plots use one setup call and repeated data calls:

```julia
fig = setup_plot(project, measurements)
plot_data(project, fig, measurements)
```

Debug helpers should mirror the real dispatch path:

```julia
debug_file(file)
debug_measurement(measurement)
```

`debug_file` should show source metadata, `read_source_file` cache status, interpreted
measurements, and any errors with stack traces.

## Migration Target

First migration step:

```text
add read_source_file default hook
add file_contents accessor
add data_contents hook/accessor
store contents on or near SourceFile
make RuO2 PUND fatigue use file_contents
remove duplicate fatigue CSV reads
replace file-oriented plotting hooks with setup_plot / plot_data
```
