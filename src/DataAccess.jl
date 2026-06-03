using DataFrames

"""
    load_source_data(project, source_file; measurement=nothing, should_cancel=nothing)::DataFrame

Project-implemented reader for one physical source file.

When `measurement` is provided, return the table needed for that logical measurement. This lets a
project select a cycle, segment, or sub-measurement from a larger source file while keeping package
code independent of project file formats.
"""
function load_source_data(
    ::AbstractProject,
    source_file::SourceFile;
    measurement::Union{Nothing,MeasurementInfo}=nothing,
    should_cancel::Union{Nothing,Function}=nothing,
)::DataFrame
    error("No source data loader for $(source_file.filepath)")
end

"""
    data_of_measurements(project, measurements; should_cancel=nothing)::Vector{DataFrame}

Package-owned data access for selected measurements.

Returns one `DataFrame` per `MeasurementInfo`, in the same order. Today this opens each measurement's
source file through `load_source_data`; after the cache rework this is the place that will read valid
cached data first and only open source files when needed.
"""
function data_of_measurements(
    project::AbstractProject,
    measurements::Vector{MeasurementInfo};
    should_cancel::Union{Nothing,Function}=nothing,
)::Vector{DataFrame}
    data = DataFrame[]
    sizehint!(data, length(measurements))
    for measurement in measurements
        _check_plot_cancel(should_cancel)
        source_file = index_source_file(measurement.filepath)
        push!(data, load_source_data(project, source_file; measurement, should_cancel))
    end
    return data
end
