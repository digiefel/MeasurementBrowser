"""
Package-owned data access for selected measurements.

Returns one `DataFrame` per `MeasurementInfo`, in the same order. Valid cached data is used first.
When the cache has no fresh data for a measurement, the source file is indexed and loaded through
`load_source_data`.
"""
function read_measurement_data(
    workspace::Workspace,
    measurements::Vector{MeasurementInfo},
)::Vector{DataFrame}
    isempty(measurements) && return DataFrame[]
    result = Vector{DataFrame}(undef, length(measurements))
    unresolved_positions = Int[]
    fingerprints = Dict{String,FileFingerprint}()
    for (position, measurement) in pairs(measurements)
        check_cancel()
        path = normpath(abspath(expanduser(measurement.filepath)))
        fingerprint = get!(fingerprints, path) do
            file_fingerprint(path)
        end
        loaded = get(workspace.direct_data, measurement.unique_id, nothing)
        if loaded === nothing || first(loaded) != fingerprint
            push!(unresolved_positions, position)
        else
            result[position] = last(loaded)
        end
    end

    unresolved = measurements[unresolved_positions]
    cached = cached_measurement_data(workspace.cache.index, unresolved)
    missing_measurements = MeasurementInfo[]
    missing_data = DataFrame[]
    source_files = Dict{String,SourceFile}()
    for (offset, position) in pairs(unresolved_positions)
        check_cancel()
        measurement = measurements[position]
        if cached[offset] === nothing
            source_file = get!(source_files, measurement.filepath) do
                index_source_file(measurement.filepath)
            end
            data = load_source_data(
                workspace.project,
                source_file;
                measurement,
            )
            push!(missing_measurements, measurement)
            push!(missing_data, data)
            result[position] = data
            workspace.direct_data[measurement.unique_id] = (source_file.fingerprint, data)
        else
            data = cached[offset]
            result[position] = data
            workspace.direct_data[measurement.unique_id] = (
                fingerprints[normpath(abspath(expanduser(measurement.filepath)))],
                data,
            )
        end
    end
    write_measurement_data_cache!(
        workspace.cache.index,
        missing_measurements,
        missing_data,
    )
    return result
end

"""
Package-owned processed data access for selected measurements.

Returns one processed `DataFrame` per `MeasurementInfo`, in the same order. Valid cached processed
data is used first. Missing entries are produced by the project's single-measurement
`process_measurement_data(workspace, measurement)` method.
"""
function process_measurement_data(
    workspace::Workspace,
    measurements::Vector{MeasurementInfo},
)::Vector{DataFrame}
    isempty(measurements) && return DataFrame[]
    result = Vector{DataFrame}(undef, length(measurements))
    unresolved_positions = Int[]
    fingerprints = Dict{String,FileFingerprint}()
    for (position, measurement) in pairs(measurements)
        check_cancel()
        path = normpath(abspath(expanduser(measurement.filepath)))
        fingerprint = get!(fingerprints, path) do
            file_fingerprint(path)
        end
        loaded = get(workspace.processed_data, measurement.unique_id, nothing)
        if loaded === nothing || first(loaded) != fingerprint
            push!(unresolved_positions, position)
        else
            result[position] = last(loaded)
        end
    end

    unresolved = measurements[unresolved_positions]
    cached = cached_measurement_data(workspace.cache.index, unresolved; processed=true)
    missing_offsets = findall(isnothing, cached)
    missing_measurements = unresolved[missing_offsets]
    processed_data = DataFrame[]
    sizehint!(processed_data, length(missing_offsets))
    for (position, offset) in pairs(missing_offsets)
        check_cancel()
        measurement = missing_measurements[position]
        data = process_measurement_data(workspace, measurement)
        push!(processed_data, data)
        result[unresolved_positions[offset]] = data
        workspace.processed_data[measurement.unique_id] = (
            fingerprints[normpath(abspath(expanduser(measurement.filepath)))],
            data,
        )
    end
    for offset in eachindex(unresolved)
        cached[offset] === nothing && continue
        measurement = unresolved[offset]
        data = cached[offset]
        result[unresolved_positions[offset]] = data
        workspace.processed_data[measurement.unique_id] = (
            fingerprints[normpath(abspath(expanduser(measurement.filepath)))],
            data,
        )
    end
    write_measurement_data_cache!(
        workspace.cache.index,
        missing_measurements,
        processed_data;
        processed=true,
    )
    return result
end

"""
Return direct data when a project defines no additional processing.
"""
function process_measurement_data(
    workspace::Workspace,
    measurement::MeasurementInfo,
)::DataFrame
    return only(read_measurement_data(workspace, [measurement]))
end
