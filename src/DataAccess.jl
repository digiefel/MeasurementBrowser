using DataFrames

"""
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
Package-owned data access for selected measurements.

Returns one `DataFrame` per `MeasurementInfo`, in the same order. Valid cached data is used first.
When the cache has no fresh data for a measurement, the source file is indexed and loaded through
`load_source_data`.
"""
function read_measurement_data(
    project::AbstractProject,
    measurements::Vector{MeasurementInfo},
)::Vector{DataFrame}
    isempty(measurements) && return DataFrame[]
    cached = cached_measurement_data(project, measurements)
    result = Vector{DataFrame}(undef, length(measurements))
    missing_measurements = MeasurementInfo[]
    missing_data = DataFrame[]
    source_files = Dict{String,SourceFile}()
    for (index, measurement) in pairs(measurements)
        _check_cancel()
        if cached[index] === nothing
            source_file = get!(source_files, measurement.filepath) do
                index_source_file(measurement.filepath)
            end
            data = load_source_data(
                project,
                source_file;
                measurement,
            )
            push!(missing_measurements, measurement)
            push!(missing_data, data)
            result[index] = data
        else
            result[index] = cached[index]
        end
    end
    write_measurement_data_cache!(project, missing_measurements, missing_data)
    return result
end

"""
Package-owned processed data access for selected measurements.

Returns one processed `DataFrame` per `MeasurementInfo`, in the same order. Valid cached processed
data is used first. Missing processed data is produced from `read_measurement_data` and the
single-measurement `process_measurement_data(project, measurement, data)` method.
"""
function process_measurement_data(
    project::AbstractProject,
    measurements::Vector{MeasurementInfo},
)::Vector{DataFrame}
    isempty(measurements) && return DataFrame[]
    cached = cached_measurement_data(project, measurements; processed=true)
    result = Vector{DataFrame}(undef, length(measurements))
    missing = findall(isnothing, cached)
    missing_measurements = measurements[missing]
    direct_data = read_measurement_data(project, missing_measurements)
    processed_data = DataFrame[]
    sizehint!(processed_data, length(missing))
    for (position, index) in pairs(missing)
        _check_cancel()
        data = process_measurement_data(
            project,
            measurements[index],
            direct_data[position],
        )
        push!(processed_data, data)
        result[index] = data
    end
    for index in eachindex(measurements)
        cached[index] === nothing || (result[index] = cached[index])
    end
    write_measurement_data_cache!(
        project,
        missing_measurements,
        processed_data;
        processed=true,
    )
    return result
end

"""
Project-implemented conversion from direct measurement data to processed measurement data.

The default returns `data` unchanged.
"""
function process_measurement_data(
    ::AbstractProject,
    ::MeasurementInfo,
    data::DataFrame,
)::DataFrame
    return data
end
