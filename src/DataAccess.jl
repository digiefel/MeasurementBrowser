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
    read_measurement_data(project, measurements)::Vector{DataFrame}

Package-owned data access for selected measurements.

Returns one `DataFrame` per `MeasurementInfo`, in the same order. Valid cached data is used first.
When the cache has no fresh data for a measurement, the source file is indexed and loaded through
`load_source_data`.
"""
function read_measurement_data(
    project::AbstractProject,
    measurements::Vector{MeasurementInfo},
)::Vector{DataFrame}
    cached_data = _cached_measurements_data(project, measurements)
    return _fill_measurement_data_cache!(project, measurements, cached_data)
end

"""
    process_measurement_data(project, measurements)::Vector{DataFrame}

Package-owned processed data access for selected measurements.

Returns one processed `DataFrame` per `MeasurementInfo`, in the same order. Valid cached processed
data is used first. Missing processed data is produced from `read_measurement_data` and the
single-measurement `process_measurement_data(project, measurement, data)` method.
"""
function process_measurement_data(
    project::AbstractProject,
    measurements::Vector{MeasurementInfo},
)::Vector{DataFrame}
    cached_data = _cached_measurements_data(project, measurements; processed=true)
    direct_data = any(isnothing, cached_data) ? read_measurement_data(project, measurements) : DataFrame[]
    return _fill_measurement_data_cache!(project, measurements, cached_data; processed=true) do index
        return process_measurement_data(project, measurements[index], direct_data[index])
    end
end

"""
    process_measurement_data(project, measurement, data)::DataFrame

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

function _fill_measurement_data_cache!(
    project::AbstractProject,
    measurements::Vector{MeasurementInfo},
    cached_data::Vector{Union{Nothing,DataFrame}},
    ;
    processed::Bool=false,
)::Vector{DataFrame}
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
        _write_cached_measurement_data!(project, measurement, loaded; processed)
        push!(data, loaded)
    end
    return data
end

function _fill_measurement_data_cache!(
    load_data::Function,
    project::AbstractProject,
    measurements::Vector{MeasurementInfo},
    cached_data::Vector{Union{Nothing,DataFrame}},
    ;
    processed::Bool=false,
)::Vector{DataFrame}
    data = DataFrame[]
    sizehint!(data, length(measurements))
    for (index, measurement) in pairs(measurements)
        _check_cancel()
        cached = cached_data[index]
        if cached !== nothing
            push!(data, cached)
            continue
        end
        loaded = load_data(index)
        _write_cached_measurement_data!(project, measurement, loaded; processed)
        push!(data, loaded)
    end
    return data
end
