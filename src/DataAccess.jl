using DataFrames

"""
    load_source_data(project, source_file; measurement=nothing)::DataFrame

Project-implemented reader for one physical source file.

When `measurement` is provided, return the table needed for that logical measurement. This lets a
project select a cycle, segment, or sub-measurement from a larger source file while keeping package
code independent of project file formats.
"""
function load_source_data(
    ::AbstractProject,
    source_file::SourceFile;
    measurement::Union{Nothing,MeasurementInfo}=nothing,
)::DataFrame
    error("No source data loader for $(source_file.filepath)")
end

"""
    data_of_measurements(project, measurements)::Vector{DataFrame}

Package-owned data access for selected measurements.

Returns one `DataFrame` per `MeasurementInfo`, in the same order. Valid cached data is used first.
When the cache has no fresh data for a measurement, the source file is indexed and loaded through
`load_source_data`.
"""
function data_of_measurements(
    project::AbstractProject,
    measurements::Vector{MeasurementInfo},
)::Vector{DataFrame}
    cached_data = _cached_measurements_data(project, measurements)
    data = DataFrame[]
    sizehint!(data, length(measurements))
    for (index, measurement) in pairs(measurements)
        _check_cancel()
        cached = cached_data[index]
        if cached !== nothing
            push!(data, cached)
            continue
        end
        source_file = index_source_file(measurement.filepath)
        loaded = load_source_data(project, source_file; measurement)
        _write_cached_measurement_data!(project, measurement, loaded)
        push!(data, loaded)
    end
    return data
end
