using HDF5
using SHA
using Serialization
using DataFrames: DataFrame

const _CACHE_DATASET_DEFLATE = 3
const _CACHE_CHUNK_TARGET = 4096
const _CACHE_SCHEMA_VERSION = 1

struct ProjectCacheMissingError <: Exception
    path::String
end

Base.showerror(io::IO, err::ProjectCacheMissingError) =
    print(io, "Project cache does not exist: $(err.path)")

struct ProjectCacheInvalidError <: Exception
    path::String
    message::String
end

Base.showerror(io::IO, err::ProjectCacheInvalidError) =
    print(io, "Invalid project cache $(err.path): $(err.message)")

struct ProjectCacheIdentity
    cache_id::String
    project_name::String
    root_path::String
    cache_path::String
end

struct ProjectCacheStatus
    total_files::Int
    cached_files::Int
    fresh_files::Int
    stale_files::Int
    new_files::Int
    deleted_files::Int
    error_files::Int
end

struct ProjectCacheFileError
    path::String
    message::String
end

struct ProjectCacheSnapshot
    identity::ProjectCacheIdentity
    hierarchy::MeasurementHierarchy
    status::ProjectCacheStatus
    errors::Vector{ProjectCacheFileError}
end

function project_cache_dir()
    depot = isempty(DEPOT_PATH) ? homedir() : DEPOT_PATH[1]
    return joinpath(depot, "measurementbrowser", "cache")
end

function new_project_cache_id()
    return Dates.format(now(), dateformat"yyyymmdd_HHMMSS")
end

function _cache_slug(value::AbstractString)
    normalized = lowercase(String(value))
    return replace(normalized, r"[^a-z0-9_.-]+" => "_")
end

function _cache_normalize_path(path::AbstractString)
    return normpath(abspath(expanduser(String(path))))
end

function project_cache_id(root_path::AbstractString)
    norm_root = _cache_normalize_path(root_path)
    name = basename(norm_root)
    slug = isempty(name) ? "project" : _cache_slug(name)
    digest = bytes2hex(sha1(norm_root))[1:12]
    return "$slug-$digest"
end

function project_cache_path(cache_id::AbstractString, project::AbstractProject)
    return joinpath(
        project_cache_dir(),
        _cache_slug(project_name(project)),
        "$(String(cache_id)).h5",
    )
end

function project_cache_identity(
    cache_id::AbstractString,
    project::AbstractProject,
    root_path::AbstractString,
)
    norm_root = _cache_normalize_path(root_path)
    return ProjectCacheIdentity(
        String(cache_id),
        project_name(project),
        norm_root,
        project_cache_path(cache_id, project),
    )
end

function _hash_key(value::AbstractString)
    return bytes2hex(sha1(String(value)))
end

_file_group_key(path::AbstractString) = _hash_key(_cache_normalize_path(path))
_measurement_group_key(id::AbstractString) = _hash_key(id)

function _dataset_chunk(data)
    data isa AbstractVector || return nothing
    isempty(data) && return nothing
    return (min(length(data), _CACHE_CHUNK_TARGET),)
end

function _write_dataset!(group, name::AbstractString, data; compress::Bool=true)
    haskey(group, name) && HDF5.delete_object(group, name)
    chunk = _dataset_chunk(data)
    if compress && chunk !== nothing && !(eltype(data) <: AbstractString)
        write(group, name, data; chunk, deflate=_CACHE_DATASET_DEFLATE)
    else
        write(group, name, data)
    end
    return nothing
end

function _ensure_group(parent, name::AbstractString)
    return haskey(parent, name) ? parent[name] : create_group(parent, name)
end

function _replace_group(parent, name::AbstractString)
    haskey(parent, name) && HDF5.delete_object(parent, name)
    return create_group(parent, name)
end

function _read_required(parent, path::AbstractString, cache_path::AbstractString)
    haskey(parent, path) || throw(ProjectCacheInvalidError(cache_path, "missing dataset '$path'"))
    return read(parent[path])
end

function _write_string_vector!(group, name::AbstractString, values::AbstractVector{<:AbstractString})
    _write_dataset!(group, name, String.(values); compress=false)
end

function _serialized_bytes(value)
    io = IOBuffer()
    serialize(io, value)
    return take!(io)
end

function _write_serialized!(group, name::AbstractString, value)::Nothing
    _write_dataset!(group, name, _serialized_bytes(value))
    return nothing
end

function _read_serialized(group, name::AbstractString, cache_path::AbstractString)
    bytes = _read_required(group, name, cache_path)
    return deserialize(IOBuffer(bytes))
end

function _write_dataframe!(group, name::AbstractString, data::DataFrame)::Nothing
    _write_serialized!(group, name, data)
    return nothing
end

function _read_dataframe(group, name::AbstractString, cache_path::AbstractString)::DataFrame
    data = _read_serialized(group, name, cache_path)
    data isa DataFrame ||
        throw(ProjectCacheInvalidError(cache_path, "cached dataframe '$name' has type $(typeof(data))"))
    return data
end

function _write_fingerprint!(group, fingerprint::FileFingerprint)
    _write_dataset!(group, "path", fingerprint.path; compress=false)
    _write_dataset!(group, "size_bytes", Int64[fingerprint.size_bytes])
    _write_dataset!(group, "mtime_ns", Int64[fingerprint.mtime_ns])
    return nothing
end

function _read_fingerprint(group, cache_path::AbstractString)
    path = _read_required(group, "path", cache_path)
    sizes = _read_required(group, "size_bytes", cache_path)
    mtimes = _read_required(group, "mtime_ns", cache_path)
    length(sizes) == 1 && length(mtimes) == 1 ||
        throw(ProjectCacheInvalidError(cache_path, "invalid file fingerprint for '$path'"))
    return FileFingerprint(path, Int64(sizes[1]), Int64(mtimes[1]))
end

function _same_fingerprint(a::FileFingerprint, b::FileFingerprint)
    return a.path == b.path && a.size_bytes == b.size_bytes && a.mtime_ns == b.mtime_ns
end

const _ACTIVE_PROJECT_CACHE = Ref{Union{Nothing,ProjectCacheIdentity}}(nothing)

function _set_active_project_cache!(identity::Union{Nothing,ProjectCacheIdentity})::Nothing
    _ACTIVE_PROJECT_CACHE[] = identity
    return nothing
end

function _write_meta!(h5, identity::ProjectCacheIdentity)
    meta = _replace_group(h5, "meta")
    _write_dataset!(meta, "schema_version", Int64[_CACHE_SCHEMA_VERSION])
    _write_dataset!(meta, "project_name", identity.project_name; compress=false)
    _write_dataset!(meta, "cache_id", identity.cache_id; compress=false)
    _write_dataset!(meta, "root_path", identity.root_path; compress=false)
    _write_dataset!(meta, "has_device_metadata", Bool[isfile(joinpath(identity.root_path, "device_info.txt"))])
    _write_dataset!(meta, "updated_at", Dates.format(now(), dateformat"yyyy-mm-ddTHH:MM:SS.s"); compress=false)
    return nothing
end

function _validate_meta!(h5, identity::ProjectCacheIdentity)
    haskey(h5, "meta") || throw(ProjectCacheInvalidError(identity.cache_path, "missing /meta group"))
    meta = h5["meta"]
    haskey(meta, "schema_version") ||
        throw(ProjectCacheInvalidError(identity.cache_path, "cache schema is out of date"))
    schema_values = read(meta["schema_version"])
    length(schema_values) == 1 && Int(schema_values[1]) == _CACHE_SCHEMA_VERSION ||
        throw(ProjectCacheInvalidError(identity.cache_path, "cache schema is out of date"))
    project_name_value = _read_required(meta, "project_name", identity.cache_path)
    project_name_value == identity.project_name ||
        throw(ProjectCacheInvalidError(identity.cache_path, "project mismatch: expected $(identity.project_name), found $project_name_value"))
    cache_id_value = _read_required(meta, "cache_id", identity.cache_path)
    cache_id_value == identity.cache_id ||
        @warn "Cache ID mismatch; expected $(identity.cache_id), found $cache_id_value"
    return nothing
end

function _read_cache_has_device_metadata(h5, cache_path::AbstractString)
    haskey(h5, "meta") || throw(ProjectCacheInvalidError(cache_path, "missing /meta group"))
    meta = h5["meta"]
    values = _read_required(meta, "has_device_metadata", cache_path)
    length(values) == 1 ||
        throw(ProjectCacheInvalidError(cache_path, "invalid has_device_metadata value"))
    return Bool(values[1])
end

function _write_file_group!(
    files_group,
    source::SourceFile,
    project::AbstractProject,
)
    measurements = MeasurementInfo[source.measurements...]
    file_group = _replace_group(files_group, _file_group_key(source.fingerprint.path))
    _write_fingerprint!(file_group, source.fingerprint)
    _write_dataset!(file_group, "source_file_unique_id", source.unique_id; compress=false)

    measurement_keys = String[]
    measurements_group = _ensure_group(file_group, "measurements")
    for measurement in measurements
        measurement_key = _measurement_group_key(measurement.unique_id)
        push!(measurement_keys, measurement_key)
        _replace_group(measurements_group, measurement_key)
    end
    _write_string_vector!(file_group, "measurement_keys", measurement_keys)
    _write_dataset!(file_group, "status", isempty(measurements) ? "skipped" : "ok"; compress=false)
    return length(measurements)
end

function _write_cached_measurement_data!(
    project::AbstractProject,
    measurement::MeasurementInfo,
    data::DataFrame,
    ;
    processed::Bool=false,
)::Nothing
    identity = _ACTIVE_PROJECT_CACHE[]
    identity === nothing && return nothing
    identity.project_name == project_name(project) || return nothing
    isfile(identity.cache_path) || return nothing

    path = _cache_normalize_path(measurement.filepath)
    isfile(path) || return nothing
    current = file_fingerprint(path)

    h5open(identity.cache_path, "r+") do h5
        _validate_meta!(h5, identity)
        haskey(h5, "files") || return nothing
        files_group = h5["files"]
        file_key = _file_group_key(path)
        haskey(files_group, file_key) || return nothing
        file_group = files_group[file_key]
        _file_status(file_group) in ("ok", "analysis_error") || return nothing
        cached = _read_fingerprint(file_group, identity.cache_path)
        _same_fingerprint(current, cached) || return nothing
        haskey(file_group, "measurements") || return nothing
        measurements_group = file_group["measurements"]
        measurement_key = _measurement_group_key(measurement.unique_id)
        haskey(measurements_group, measurement_key) || return nothing
        measurement_group = measurements_group[measurement_key]
        _write_dataframe!(measurement_group, _measurement_data_name(processed), data)
    end
    return nothing
end

function _cached_measurements_data(
    project::AbstractProject,
    measurements::Vector{MeasurementInfo},
    ;
    processed::Bool=false,
)::Vector{Union{Nothing,DataFrame}}
    data = Union{Nothing,DataFrame}[nothing for _ in measurements]
    identity = _ACTIVE_PROJECT_CACHE[]
    identity === nothing && return data
    identity.project_name == project_name(project) || return data
    isfile(identity.cache_path) || return data

    current_fingerprints = Dict{String,FileFingerprint}()
    h5open(identity.cache_path, "r") do h5
        _validate_meta!(h5, identity)
        haskey(h5, "files") || return data
        files_group = h5["files"]
        for (index, measurement) in pairs(measurements)
            _check_cancel()
            path = _cache_normalize_path(measurement.filepath)
            isfile(path) || continue
            current = get!(current_fingerprints, path) do
                file_fingerprint(path)
            end
            file_key = _file_group_key(path)
            haskey(files_group, file_key) || continue
            file_group = files_group[file_key]
            _file_status(file_group) in ("ok", "analysis_error") || continue
            cached = _read_fingerprint(file_group, identity.cache_path)
            _same_fingerprint(current, cached) || continue
            haskey(file_group, "measurements") || continue
            measurements_group = file_group["measurements"]
            measurement_key = _measurement_group_key(measurement.unique_id)
            haskey(measurements_group, measurement_key) || continue
            measurement_group = measurements_group[measurement_key]
            data_name = _measurement_data_name(processed)
            haskey(measurement_group, data_name) || continue
            data[index] = _read_dataframe(measurement_group, data_name, identity.cache_path)
        end
    end
    return data
end

_measurement_data_name(processed::Bool)::String = processed ? "processed_data" : "data"

function _file_status(file_group)
    return haskey(file_group, "status") ? read(file_group["status"]) : "ok"
end

function _cached_file_fingerprints(h5, cache_path::AbstractString)
    haskey(h5, "files") || return Dict{String,FileFingerprint}()
    fingerprints = Dict{String,FileFingerprint}()
    for file_key in keys(h5["files"])
        group = h5["files"][file_key]
        fp = _read_fingerprint(group, cache_path)
        fingerprints[fp.path] = fp
    end
    return fingerprints
end

function _cached_file_statuses(h5)
    haskey(h5, "files") || return Dict{String,String}()
    statuses = Dict{String,String}()
    for file_key in keys(h5["files"])
        group = h5["files"][file_key]
        path = haskey(group, "path") ? read(group["path"]) : file_key
        statuses[path] = _file_status(group)
    end
    return statuses
end

function _read_cache_file_errors(h5, cache_path::AbstractString)
    haskey(h5, "files") || return ProjectCacheFileError[]
    errors = ProjectCacheFileError[]
    for file_key in sort!(collect(keys(h5["files"])))
        group = h5["files"][file_key]
        _file_status(group) in ("error", "analysis_error") || continue
        fp = _read_fingerprint(group, cache_path)
        message = haskey(group, "error_message") ? read(group["error_message"]) :
            "Cache transform failed without an error message"
        push!(errors, ProjectCacheFileError(fp.path, message))
    end
    return errors
end

function _cache_file_errors(identity::ProjectCacheIdentity)
    isfile(identity.cache_path) || throw(ProjectCacheMissingError(identity.cache_path))
    return h5open(identity.cache_path, "r") do h5
        _validate_meta!(h5, identity)
        _read_cache_file_errors(h5, identity.cache_path)
    end
end

function _cache_status_from_fingerprints(
    identity::ProjectCacheIdentity,
    raw::Dict{String,FileFingerprint},
)
    isfile(identity.cache_path) || throw(ProjectCacheMissingError(identity.cache_path))
    return h5open(identity.cache_path, "r") do h5
        _validate_meta!(h5, identity)
        cached = _cached_file_fingerprints(h5, identity.cache_path)
        statuses = _cached_file_statuses(h5)
        return _cache_status_from_fingerprints(raw, cached, statuses)
    end
end

function _source_fingerprints(source::SourceScan)
    fingerprints = Dict{String,FileFingerprint}()
    for file in source.files
        fingerprints[file.fingerprint.path] = file.fingerprint
    end
    return fingerprints
end

function _source_analysis_statuses(source::SourceScan)::Dict{String,String}
    statuses = Dict{String,String}()
    for failure in source.analysis_failures
        statuses[failure.filepath] = "analysis_error"
    end
    return statuses
end

function _cache_status_from_fingerprints(
    raw::Dict{String,FileFingerprint},
    cached::Dict{String,FileFingerprint},
    statuses::Dict{String,String},
)
    stale = count(path -> haskey(raw, path) && !_same_fingerprint(raw[path], cached[path]), keys(cached))
    deleted = count(path -> !haskey(raw, path), keys(cached))
    new_count = count(path -> !haskey(cached, path), keys(raw))
    fresh = count(
        path -> haskey(cached, path) &&
            _same_fingerprint(raw[path], cached[path]) &&
            get(statuses, path, "ok") in ("ok", "skipped"),
        keys(raw),
    )
    error_count = count(status -> status in ("error", "analysis_error"), values(statuses))
    return ProjectCacheStatus(
        length(raw),
        length(cached),
        fresh,
        stale,
        new_count,
        deleted,
        error_count,
    )
end

function _cached_cache_index(identity::ProjectCacheIdentity)
    isfile(identity.cache_path) || throw(ProjectCacheMissingError(identity.cache_path))
    return h5open(identity.cache_path, "r") do h5
        _validate_meta!(h5, identity)
        (_cached_file_fingerprints(h5, identity.cache_path), _cached_file_statuses(h5))
    end
end

function cache_status(identity::ProjectCacheIdentity, source::SourceScan)
    identity.root_path == source.root_path ||
        error("Cache root $(identity.root_path) does not match source scan root $(source.root_path)")
    project_name(source.project) == identity.project_name ||
        error("Cache project $(identity.project_name) does not match source scan project $(project_name(source.project))")
    cached, statuses = _cached_cache_index(identity)
    merge!(statuses, _source_analysis_statuses(source))
    return _cache_status_from_fingerprints(_source_fingerprints(source), cached, statuses)
end

function _cache_status_from_cached_statuses(statuses::Dict{String,String})
    error_count = count(status -> status in ("error", "analysis_error"), values(statuses))
    cached = length(statuses)
    return ProjectCacheStatus(cached, cached, cached - error_count, 0, 0, 0, error_count)
end

function _hierarchy_from_cache_statuses(source::SourceScan, statuses::Dict{String,String})
    measurements = MeasurementInfo[]
    skipped_count = 0
    for source_file in source.files
        status = get(statuses, source_file.fingerprint.path, "missing")
        if status in ("ok", "error", "analysis_error")
            append!(measurements, source_file.measurements)
        elseif status == "skipped"
            skipped_count += 1
        end
    end
    return MeasurementHierarchy(
        measurements,
        source.root_path,
        source.hierarchy.has_device_metadata,
        source.project,
        skipped_count,
    )
end

const _CACHE_INDEX_VERSION = 3

function _cache_index_group(h5, cache_path::AbstractString)
    haskey(h5, "indexes") || throw(ProjectCacheInvalidError(cache_path, "missing /indexes group"))
    index_group = h5["indexes"]
    version_values = _read_required(index_group, "version", cache_path)
    length(version_values) == 1 && Int(version_values[1]) == _CACHE_INDEX_VERSION ||
        throw(ProjectCacheInvalidError(cache_path, "cache index is out of date"))
    return index_group
end

function _cache_startup_blob(hierarchy::MeasurementHierarchy, statuses::Dict{String,String}, errors::Vector{ProjectCacheFileError})
    return (
        measurements=hierarchy.all_measurements,
        statuses=statuses,
        errors=errors,
        skipped_count=hierarchy.skipped_count,
    )
end

function _write_cache_startup_blob!(group, blob)::Nothing
    _write_serialized!(group, "startup_blob", blob)
    return nothing
end

function _read_cache_startup_blob(h5, cache_path::AbstractString)
    index_group = _cache_index_group(h5, cache_path)
    return _read_serialized(index_group, "startup_blob", cache_path)
end

function _rebuild_indexes!(h5, hierarchy::MeasurementHierarchy, statuses::Dict{String,String})
    index_group = _replace_group(h5, "indexes")
    _write_dataset!(index_group, "version", Int64[_CACHE_INDEX_VERSION])

    errors = ProjectCacheFileError[]
    files_group = h5["files"]
    for path in sort!([path for (path, status) in statuses if status in ("error", "analysis_error")])
        file_key = _file_group_key(path)
        haskey(files_group, file_key) || continue
        file_group = files_group[file_key]
        message = haskey(file_group, "error_message") ? read(file_group["error_message"]) : "unknown cache error"
        push!(errors, ProjectCacheFileError(path, message))
    end

    _write_cache_startup_blob!(index_group, _cache_startup_blob(hierarchy, statuses, errors))
    return hierarchy
end

function _write_file_error_group!(
    files_group,
    fingerprint::FileFingerprint,
    source::SourceFile,
    err,
)
    file_group = _replace_group(files_group, _file_group_key(fingerprint.path))
    _write_fingerprint!(file_group, fingerprint)
    _write_dataset!(file_group, "source_file_unique_id", source.unique_id; compress=false)
    measurement_keys = String[]
    measurements_group = _ensure_group(file_group, "measurements")
    for measurement in source.measurements
        measurement_key = _measurement_group_key(measurement.unique_id)
        push!(measurement_keys, measurement_key)
        _replace_group(measurements_group, measurement_key)
    end
    _write_string_vector!(file_group, "measurement_keys", measurement_keys)
    _write_dataset!(file_group, "status", "error"; compress=false)
    _write_dataset!(file_group, "error_type", string(typeof(err)); compress=false)
    _write_dataset!(file_group, "error_message", sprint(showerror, err); compress=false)
    return nothing
end

function _write_file_analysis_failures!(
    files_group,
    source::SourceFile,
    failures::Vector{MeasurementAnalysisFailure},
)::Nothing
    isempty(failures) && return nothing
    file_group = files_group[_file_group_key(source.fingerprint.path)]
    _write_dataset!(file_group, "status", "analysis_error"; compress=false)
    _write_dataset!(
        file_group,
        "error_message",
        join([failure.message for failure in failures], "\n");
        compress=false,
    )
    return nothing
end

function write_project_cache!(
    identity::ProjectCacheIdentity,
    source::SourceScan;
    mode::Symbol=:update,
    on_progress::Union{Nothing,Function}=nothing,
)
    mode in (:build, :rebuild, :update) || error("Unsupported cache write mode '$mode'")
    identity.root_path == source.root_path ||
        error("Cache root $(identity.root_path) does not match source scan root $(source.root_path)")
    project_name(source.project) == identity.project_name ||
        error("Cache project $(identity.project_name) does not match source scan project $(project_name(source.project))")
    mkpath(dirname(identity.cache_path))
    @info "Starting project cache write" root = identity.root_path cache = identity.cache_path mode
    cache_exists = isfile(identity.cache_path)
    h5_mode = (!cache_exists || mode in (:build, :rebuild)) ? "w" : "r+"
    processed = 0
    loaded_measurements = 0
    fingerprints = _source_fingerprints(source)
    source_by_path = Dict(file.fingerprint.path => file for file in source.files)
    analysis_failures = Dict{String,Vector{MeasurementAnalysisFailure}}()
    for failure in source.analysis_failures
        push!(get!(analysis_failures, failure.filepath, MeasurementAnalysisFailure[]), failure)
    end

    hierarchy = h5open(identity.cache_path, h5_mode) do h5
        if cache_exists && mode == :update
            _validate_meta!(h5, identity)
        else
            _write_meta!(h5, identity)
        end
        files_group = _ensure_group(h5, "files")
        cached = _cached_file_fingerprints(h5, identity.cache_path)
        statuses = cache_exists ? _cached_file_statuses(h5) : Dict{String,String}()
        raw_paths = Set(keys(fingerprints))
        deleted_paths = sort!([path for path in keys(cached) if !(path in raw_paths)])
        write_paths = if !cache_exists || mode in (:build, :rebuild)
            sort!(collect(keys(fingerprints)))
        else
            sort!([
                path for path in keys(fingerprints)
                if !haskey(cached, path) ||
                    !_same_fingerprint(fingerprints[path], cached[path]) ||
                    get(statuses, path, "missing") in ("error", "analysis_error") ||
                    haskey(analysis_failures, path)
            ])
        end
        total_changes = length(deleted_paths) + length(write_paths)

        for cached_path in deleted_paths
            _check_cancel()
            cached_path in raw_paths && continue
            key = _file_group_key(cached_path)
            haskey(files_group, key) && HDF5.delete_object(files_group, key)
            delete!(statuses, cached_path)
            processed += 1
            on_progress !== nothing && on_progress((
                phase=:cache_update,
                total_csv=total_changes,
                processed_csv=processed,
                loaded_measurements=loaded_measurements,
                skipped_csv=0,
                current_path=cached_path,
            ))
        end

        for path in write_paths
            _check_cancel()
            source_file = source_by_path[path]
            try
                measurement_count = _write_file_group!(
                    files_group,
                    source_file,
                    source.project,
                )
                loaded_measurements += measurement_count
                statuses[path] = measurement_count == 0 ? "skipped" : "ok"
                failures = get(analysis_failures, path, MeasurementAnalysisFailure[])
                if !isempty(failures)
                    _write_file_analysis_failures!(files_group, source_file, failures)
                    statuses[path] = "analysis_error"
                end
            catch err
                err isa JobCancelled && rethrow()
                _write_file_error_group!(files_group, source_file.fingerprint, source_file, err)
                statuses[path] = "error"
            end
            processed += 1
            on_progress !== nothing && on_progress((
                phase=:cache_update,
                total_csv=total_changes,
                processed_csv=processed,
                loaded_measurements=loaded_measurements,
                skipped_csv=0,
                current_path=path,
            ))
        end
        on_progress !== nothing && on_progress((
            phase=:cache_finalize,
            total_csv=1,
            processed_csv=0,
            loaded_measurements=loaded_measurements,
            skipped_csv=0,
            current_path=identity.cache_path,
        ))
        _write_meta!(h5, identity)
        hierarchy = _hierarchy_from_cache_statuses(source, statuses)
        _rebuild_indexes!(h5, hierarchy, statuses)
        on_progress !== nothing && on_progress((
            phase=:cache_finalize,
            total_csv=1,
            processed_csv=1,
            loaded_measurements=loaded_measurements,
            skipped_csv=0,
            current_path=identity.cache_path,
        ))
        hierarchy
    end
    status = _cache_status_from_fingerprints(identity, fingerprints)
    errors = _cache_file_errors(identity)
    @info "Finished project cache write" cache = identity.cache_path cached_files = status.cached_files error_files = status.error_files
    return ProjectCacheSnapshot(
        identity,
        hierarchy,
        status,
        errors,
    )
end

function _load_project_cache_contents(
    identity::ProjectCacheIdentity,
    ;
    on_progress::Union{Nothing,Function}=nothing,
)
    project = _project_by_name(identity.project_name)
    isfile(identity.cache_path) || throw(ProjectCacheMissingError(identity.cache_path))
    _check_cancel()
    snapshot = h5open(identity.cache_path, "r") do h5
        _validate_meta!(h5, identity)
        has_device_metadata = _read_cache_has_device_metadata(h5, identity.cache_path)
        startup = _read_cache_startup_blob(h5, identity.cache_path)
        measurements = startup.measurements
        statuses = startup.statuses
        status = _cache_status_from_cached_statuses(statuses)
        errors = startup.errors
        skipped_count = count(==("skipped"), values(statuses))
        _emit_progress(on_progress;
            phase=:cache_load,
            total_csv=1,
            processed_csv=1,
            loaded_measurements=length(measurements),
            skipped_csv=skipped_count,
            current_path=identity.cache_path,
        )
        hierarchy = MeasurementHierarchy(
            measurements,
            identity.root_path,
            has_device_metadata,
            project,
            skipped_count,
        )
        return ProjectCacheSnapshot(identity, hierarchy, status, errors)
    end
    _check_cancel()
    return snapshot
end

"""
    load_project_cache(identity; on_progress=nothing)

Load a cache snapshot for browsing without touching the source CSV tree.
"""
function load_project_cache(
    identity::ProjectCacheIdentity,
    ;
    on_progress::Union{Nothing,Function}=nothing,
)
    return _load_project_cache_contents(identity; on_progress)
end

function _project_by_name(name::AbstractString)
    for project in KNOWN_PROJECTS
        project_name(project) == name && return project
    end
    throw(ProjectCacheInvalidError("", "unknown project '$name'"))
end
