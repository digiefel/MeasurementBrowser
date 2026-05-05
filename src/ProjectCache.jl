using HDF5
using SHA
using DataFrames

const CACHE_SCHEMA_VERSION = 1
const _CACHE_DATASET_DEFLATE = 3
const _CACHE_CHUNK_TARGET = 4096

struct ProjectCacheUnsupportedError <: Exception
    project_name::String
end

Base.showerror(io::IO, err::ProjectCacheUnsupportedError) =
    print(io, "Project '$(err.project_name)' does not define an HDF5 cache schema")

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

struct ProjectCacheBuildError <: Exception
    path::String
    cause::Exception
    bt
end

function Base.showerror(io::IO, err::ProjectCacheBuildError)
    print(io, "Could not build cache entry for $(err.path): ")
    showerror(io, err.cause)
end

struct FileFingerprint
    path::String
    size_bytes::Int64
    mtime_ns::Int64
end

struct ProjectCacheIdentity
    cache_id::String
    project_name::String
    project_schema_version::Int
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
    semantic_fields::Dict{Symbol,Vector{Symbol}}
    errors::Vector{ProjectCacheFileError}
end

project_cache_schema_version(project::AbstractProject) =
    throw(ProjectCacheUnsupportedError(project_name(project)))

project_cache_semantic_fields(::AbstractProject) = Dict{Symbol,Vector{Symbol}}()

project_cache_file_matches(::AbstractProject, ::IndexedCsvFile) = true

function project_cache_transform(
    project::AbstractProject,
    indexed::IndexedCsvFile,
    meta::Union{Nothing,Dict{Tuple{Vararg{String}},Dict{Symbol,Any}}};
    should_cancel::Union{Nothing,Function}=nothing,
)
    throw(ProjectCacheUnsupportedError(project_name(project)))
end

function project_cache_write_measurement_payload!(
    group,
    project::AbstractProject,
    measurement::MeasurementInfo;
    should_cancel::Union{Nothing,Function}=nothing,
)
    throw(ProjectCacheUnsupportedError(project_name(project)))
end

function project_cache_write_file_payload!(
    file_group,
    project::AbstractProject,
    indexed::IndexedCsvFile,
    measurements::Vector{MeasurementInfo};
    should_cancel::Union{Nothing,Function}=nothing,
)
    measurements_group = file_group["measurements"]
    for measurement in measurements
        measurement_group = measurements_group[_measurement_group_key(measurement.id)]
        project_cache_write_measurement_payload!(measurement_group, project, measurement; should_cancel)
    end
    return nothing
end

function project_cache_read_plot_payload(
    project::AbstractProject,
    measurement::MeasurementInfo,
    file_group,
    measurement_group,
)
    throw(ProjectCacheUnsupportedError(project_name(project)))
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
    norm_path = abspath(expanduser(String(path)))
    return ispath(norm_path) ? realpath(norm_path) : norm_path
end

function project_cache_path(cache_id::AbstractString, project::AbstractProject)
    schema = project_cache_schema_version(project)
    return joinpath(
        project_cache_dir(),
        _cache_slug(project_name(project)),
        "schema$(schema)",
        "$(String(cache_id)).h5",
    )
end

function project_cache_identity(
    cache_id::AbstractString,
    project::AbstractProject,
    root_path::AbstractString,
)
    norm_root = _cache_normalize_path(root_path)
    schema = project_cache_schema_version(project)
    return ProjectCacheIdentity(
        String(cache_id),
        project_name(project),
        schema,
        norm_root,
        project_cache_path(cache_id, project),
    )
end

function _hash_key(value::AbstractString)
    return bytes2hex(sha1(String(value)))
end

_file_group_key(path::AbstractString) = _hash_key(_cache_normalize_path(path))
_measurement_group_key(id::AbstractString) = _hash_key(id)

function file_fingerprint(path::AbstractString)
    normalized = _cache_normalize_path(path)
    stat_info = stat(normalized)
    return FileFingerprint(normalized, Int64(stat_info.size), Int64(stat_info.mtime * 1_000_000_000))
end

function collect_csv_fingerprints(
    root_path::AbstractString;
    should_cancel::Union{Nothing,Function}=nothing,
    on_progress::Union{Nothing,Function}=nothing,
    phase::Symbol=:cache_discovery,
)
    root = _cache_normalize_path(root_path)
    paths = String[]
    for (dir, _, names) in walkdir(root)
        for name in names
            _check_cancel(should_cancel)
            endswith(lowercase(name), ".csv") || continue
            path = joinpath(dir, name)
            push!(paths, path)
            _emit_progress(on_progress;
                phase,
                total_csv=0,
                processed_csv=length(paths),
                loaded_measurements=0,
                skipped_csv=0,
                current_path=path,
            )
        end
    end

    fingerprints = Dict{String,FileFingerprint}()
    fingerprint_lock = ReentrantLock()
    progress_lock = ReentrantLock()
    processed = Base.Threads.Atomic{Int}(0)
    worker_limit = Base.Semaphore(max(1, Base.Threads.nthreads()))

    @sync for path in paths
        Base.acquire(worker_limit)
        Base.Threads.@spawn begin
            try
                _check_cancel(should_cancel)
                fp = file_fingerprint(path)
                lock(fingerprint_lock) do
                    fingerprints[fp.path] = fp
                end
                processed_now = Base.Threads.atomic_add!(processed, 1) + 1
                lock(progress_lock) do
                    _emit_progress(on_progress;
                        phase,
                        total_csv=length(paths),
                        processed_csv=processed_now,
                        loaded_measurements=0,
                        skipped_csv=0,
                        current_path=fp.path,
                    )
                end
            finally
                Base.release(worker_limit)
            end
            _check_cancel(should_cancel)
        end
    end
    return fingerprints
end

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

function _write_symbol_vector!(group, name::AbstractString, values::AbstractVector{Symbol})
    _write_string_vector!(group, name, String.(values))
end

function _read_symbol_vector(group, name::AbstractString, cache_path::AbstractString)
    return Symbol.(_read_required(group, name, cache_path))
end

function _encode_cache_value(value)
    if value === nothing
        return "nothing", ""
    elseif value isa Bool
        return "bool", value ? "true" : "false"
    elseif value isa Integer
        return "int", string(Int64(value))
    elseif value isa AbstractFloat
        return "float", string(Float64(value))
    elseif value isa Symbol
        return "symbol", String(value)
    elseif value isa DateTime
        return "datetime", Dates.format(value, dateformat"yyyy-mm-ddTHH:MM:SS.s")
    elseif value isa Date
        return "date", Dates.format(value, dateformat"yyyy-mm-dd")
    elseif value isa AbstractString
        return "string", String(value)
    end
    error("Unsupported cache metadata value of type $(typeof(value))")
end

function _decode_cache_value(type_name::AbstractString, encoded::AbstractString)
    if type_name == "nothing"
        return nothing
    elseif type_name == "bool"
        encoded in ("true", "false") || error("Invalid cached Bool value '$encoded'")
        return encoded == "true"
    elseif type_name == "int"
        return parse(Int, encoded)
    elseif type_name == "float"
        return parse(Float64, encoded)
    elseif type_name == "symbol"
        return Symbol(encoded)
    elseif type_name == "datetime"
        return DateTime(encoded, dateformat"yyyy-mm-ddTHH:MM:SS.s")
    elseif type_name == "date"
        return Date(encoded, dateformat"yyyy-mm-dd")
    elseif type_name == "string"
        return String(encoded)
    end
    error("Unsupported cached metadata value type '$type_name'")
end

function _write_parameters!(group, name::AbstractString, parameters::Dict{Symbol,Any})
    param_group = _replace_group(group, name)
    keys_sorted = sort!(collect(keys(parameters)); by=String)
    names = String[]
    types = String[]
    values = String[]
    for key in keys_sorted
        type_name, encoded = _encode_cache_value(parameters[key])
        push!(names, String(key))
        push!(types, type_name)
        push!(values, encoded)
    end
    _write_string_vector!(param_group, "names", names)
    _write_string_vector!(param_group, "types", types)
    _write_string_vector!(param_group, "values", values)
    return nothing
end

function _read_parameters(group, name::AbstractString, cache_path::AbstractString)
    haskey(group, name) || return Dict{Symbol,Any}()
    param_group = group[name]
    names = _read_required(param_group, "names", cache_path)
    types = _read_required(param_group, "types", cache_path)
    values = _read_required(param_group, "values", cache_path)
    length(names) == length(types) == length(values) ||
        throw(ProjectCacheInvalidError(cache_path, "parameter arrays for '$name' have mismatched lengths"))
    params = Dict{Symbol,Any}()
    for index in eachindex(names, types, values)
        params[Symbol(names[index])] = _decode_cache_value(types[index], values[index])
    end
    return params
end

function _write_measurement_metadata!(group, measurement::MeasurementInfo)
    _write_dataset!(group, "id", measurement.id; compress=false)
    _write_dataset!(group, "filename", measurement.filename; compress=false)
    _write_dataset!(group, "filepath", measurement.filepath; compress=false)
    _write_dataset!(group, "clean_title", measurement.clean_title; compress=false)
    _write_dataset!(group, "measurement_kind", String(measurement.measurement_kind); compress=false)
    timestamp = measurement.timestamp === nothing ? "" :
        Dates.format(measurement.timestamp, dateformat"yyyy-mm-ddTHH:MM:SS.s")
    _write_dataset!(group, "timestamp", timestamp; compress=false)
    _write_string_vector!(group, "device_path", measurement.device_info.location)
    _write_parameters!(group, "device_parameters", measurement.device_info.parameters)
    _write_parameters!(group, "measurement_parameters", measurement.parameters)
    wakeup = measurement.wakeup_pulse_count === nothing ? Int64[-1] : Int64[Int64(measurement.wakeup_pulse_count)]
    _write_dataset!(group, "wakeup_pulse_count", wakeup)
    return nothing
end

function _read_measurement_metadata(group, cache_path::AbstractString)
    id = _read_required(group, "id", cache_path)
    filename = _read_required(group, "filename", cache_path)
    filepath = _read_required(group, "filepath", cache_path)
    clean_title = _read_required(group, "clean_title", cache_path)
    kind = Symbol(_read_required(group, "measurement_kind", cache_path))
    timestamp_text = _read_required(group, "timestamp", cache_path)
    timestamp = isempty(timestamp_text) ? nothing :
        DateTime(timestamp_text, dateformat"yyyy-mm-ddTHH:MM:SS.s")
    device_path = String.(_read_required(group, "device_path", cache_path))
    device_parameters = _read_parameters(group, "device_parameters", cache_path)
    measurement_parameters = _read_parameters(group, "measurement_parameters", cache_path)
    wakeup_values = _read_required(group, "wakeup_pulse_count", cache_path)
    length(wakeup_values) == 1 ||
        throw(ProjectCacheInvalidError(cache_path, "invalid wakeup_pulse_count dataset for '$id'"))
    wakeup = wakeup_values[1] < 0 ? nothing : Int(wakeup_values[1])
    return MeasurementInfo(
        id,
        filename,
        filepath,
        clean_title,
        kind,
        timestamp,
        DeviceInfo(device_path, device_parameters),
        measurement_parameters,
        wakeup,
    )
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

function _write_dataframe!(group, name::AbstractString, df::DataFrame)
    df_group = _replace_group(group, name)
    column_names = names(df)
    _write_string_vector!(df_group, "column_names", column_names)
    columns_group = _ensure_group(df_group, "columns")
    for column_name in column_names
        values = df[!, column_name]
        values isa AbstractVector || error("DataFrame column '$column_name' is not a vector")
        _write_dataset!(columns_group, column_name, collect(values))
    end
    return nothing
end

function _read_dataframe(group, name::AbstractString, cache_path::AbstractString)
    haskey(group, name) || throw(ProjectCacheInvalidError(cache_path, "missing DataFrame group '$name'"))
    df_group = group[name]
    column_names = String.(_read_required(df_group, "column_names", cache_path))
    columns_group = df_group["columns"]
    columns = Pair{Symbol,Any}[]
    for column_name in column_names
        haskey(columns_group, column_name) ||
            throw(ProjectCacheInvalidError(cache_path, "missing cached DataFrame column '$column_name'"))
        push!(columns, Symbol(column_name) => read(columns_group[column_name]))
    end
    return DataFrame(columns...)
end

function _write_meta!(h5, identity::ProjectCacheIdentity)
    meta = _replace_group(h5, "meta")
    _write_dataset!(meta, "cache_schema_version", Int64[CACHE_SCHEMA_VERSION])
    _write_dataset!(meta, "project_name", identity.project_name; compress=false)
    _write_dataset!(meta, "project_schema_version", Int64[identity.project_schema_version])
    _write_dataset!(meta, "cache_id", identity.cache_id; compress=false)
    _write_dataset!(meta, "root_path", identity.root_path; compress=false)
    _write_dataset!(meta, "has_device_metadata", Bool[isfile(joinpath(identity.root_path, "device_info.txt"))])
    _write_dataset!(meta, "updated_at", Dates.format(now(), dateformat"yyyy-mm-ddTHH:MM:SS.s"); compress=false)
    return nothing
end

function _validate_meta!(h5, identity::ProjectCacheIdentity)
    haskey(h5, "meta") || throw(ProjectCacheInvalidError(identity.cache_path, "missing /meta group"))
    meta = h5["meta"]
    schema_values = _read_required(meta, "cache_schema_version", identity.cache_path)
    length(schema_values) == 1 && schema_values[1] == CACHE_SCHEMA_VERSION ||
        throw(ProjectCacheInvalidError(identity.cache_path, "unsupported cache schema version"))
    project_name_value = _read_required(meta, "project_name", identity.cache_path)
    project_name_value == identity.project_name ||
        throw(ProjectCacheInvalidError(identity.cache_path, "project mismatch: expected $(identity.project_name), found $project_name_value"))
    project_schema_values = _read_required(meta, "project_schema_version", identity.cache_path)
    length(project_schema_values) == 1 && project_schema_values[1] == identity.project_schema_version ||
        throw(ProjectCacheInvalidError(identity.cache_path, "project schema version mismatch"))
    cache_id_value = _read_required(meta, "cache_id", identity.cache_path)
    cache_id_value == identity.cache_id ||
        throw(ProjectCacheInvalidError(identity.cache_path, "cache id mismatch"))
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

function _write_semantic_fields!(h5, semantic_fields::Dict{Symbol,Vector{Symbol}})
    group = _replace_group(h5, "semantics")
    names_sorted = sort!(collect(keys(semantic_fields)); by=String)
    _write_symbol_vector!(group, "groups", names_sorted)
    for name in names_sorted
        _write_symbol_vector!(group, String(name), semantic_fields[name])
    end
    return nothing
end

function _read_semantic_fields(h5, cache_path::AbstractString)
    haskey(h5, "semantics") || return Dict{Symbol,Vector{Symbol}}()
    group = h5["semantics"]
    names = _read_symbol_vector(group, "groups", cache_path)
    fields = Dict{Symbol,Vector{Symbol}}()
    for name in names
        fields[name] = _read_symbol_vector(group, String(name), cache_path)
    end
    return fields
end

function _write_file_group!(
    files_group,
    project::AbstractProject,
    fingerprint::FileFingerprint,
    indexed::IndexedCsvFile,
    meta::Union{Nothing,Dict{Tuple{Vararg{String}},Dict{Symbol,Any}}};
    should_cancel::Union{Nothing,Function}=nothing,
)
    file_group = _replace_group(files_group, _file_group_key(fingerprint.path))
    _write_fingerprint!(file_group, fingerprint)
    _write_dataset!(file_group, "source_file_id", indexed.id; compress=false)
    if !project_cache_file_matches(project, indexed)
        _write_string_vector!(file_group, "measurement_keys", String[])
        _write_dataset!(file_group, "status", "skipped"; compress=false)
        return 0
    end
    items = project_cache_transform(project, indexed, meta; should_cancel)
    measurements = [_measurement_info_from_item(item) for item in items]

    measurement_keys = String[]
    measurements_group = _ensure_group(file_group, "measurements")
    for measurement in measurements
        measurement_key = _measurement_group_key(measurement.id)
        push!(measurement_keys, measurement_key)
        measurement_group = _replace_group(measurements_group, measurement_key)
        _write_measurement_metadata!(measurement_group, measurement)
    end
    _write_string_vector!(file_group, "measurement_keys", measurement_keys)
    project_cache_write_file_payload!(file_group, project, indexed, measurements; should_cancel)
    _write_dataset!(file_group, "status", isempty(measurements) ? "skipped" : "ok"; compress=false)
    return length(measurements)
end

function _read_file_measurements(file_group, cache_path::AbstractString)
    haskey(file_group, "measurements") || return MeasurementInfo[]
    measurement_keys = String.(_read_required(file_group, "measurement_keys", cache_path))
    measurements_group = file_group["measurements"]
    measurements = MeasurementInfo[]
    sizehint!(measurements, length(measurement_keys))
    for measurement_key in measurement_keys
        haskey(measurements_group, measurement_key) ||
            throw(ProjectCacheInvalidError(cache_path, "missing measurement group '$measurement_key'"))
        push!(measurements, _read_measurement_metadata(measurements_group[measurement_key], cache_path))
    end
    return measurements
end

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
        _file_status(group) == "error" || continue
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
            get(statuses, path, "ok") != "error",
        keys(raw),
    )
    error_count = count(==("error"), values(statuses))
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

function cache_status(
    identity::ProjectCacheIdentity;
    should_cancel::Union{Nothing,Function}=nothing,
    on_progress::Union{Nothing,Function}=nothing,
)
    raw = collect_csv_fingerprints(identity.root_path; should_cancel, on_progress, phase=:cache_check)
    return _cache_status_from_fingerprints(identity, raw)
end

function _cache_status_from_cached_statuses(statuses::Dict{String,String})
    error_count = count(==("error"), values(statuses))
    cached = length(statuses)
    return ProjectCacheStatus(cached, cached, cached - error_count, 0, 0, 0, error_count)
end

function _rebuild_indexes!(h5, identity::ProjectCacheIdentity, project::AbstractProject)
    files_group = _ensure_group(h5, "files")
    measurements = MeasurementInfo[]
    skipped_count = 0
    for file_key in sort!(collect(keys(files_group)))
        file_group = files_group[file_key]
        status = _file_status(file_group)
        status in ("skipped", "error") && (skipped_count += 1)
        append!(measurements, _read_file_measurements(file_group, identity.cache_path))
    end
    hierarchy = MeasurementHierarchy(
        measurements,
        identity.root_path,
        _read_cache_has_device_metadata(h5, identity.cache_path),
        project,
        skipped_count,
    )
    index_group = _replace_group(h5, "indexes")
    _write_string_vector!(index_group, "measurement_ids", [m.id for m in hierarchy.all_measurements])
    _write_string_vector!(index_group, "measurement_keys", [_measurement_group_key(m.id) for m in hierarchy.all_measurements])
    _write_string_vector!(index_group, "device_keys", [device_path_key(m.device_info) for m in hierarchy.all_measurements])
    _write_string_vector!(index_group, "measurement_kinds", [String(m.measurement_kind) for m in hierarchy.all_measurements])
    _write_string_vector!(index_group, "timestamps", [
        m.timestamp === nothing ? "" : Dates.format(m.timestamp, dateformat"yyyy-mm-ddTHH:MM:SS.s")
        for m in hierarchy.all_measurements
    ])
    _write_dataset!(index_group, "skipped_count", Int64[skipped_count])
    return hierarchy
end

function _write_file_error_group!(
    files_group,
    fingerprint::FileFingerprint,
    indexed::IndexedCsvFile,
    err,
)
    file_group = _replace_group(files_group, _file_group_key(fingerprint.path))
    _write_fingerprint!(file_group, fingerprint)
    _write_dataset!(file_group, "source_file_id", indexed.id; compress=false)
    _write_string_vector!(file_group, "measurement_keys", String[])
    _write_dataset!(file_group, "status", "error"; compress=false)
    _write_dataset!(file_group, "error_type", string(typeof(err)); compress=false)
    _write_dataset!(file_group, "error_message", sprint(showerror, err); compress=false)
    return nothing
end

function build_project_cache!(
    root_path::AbstractString,
    project::AbstractProject,
    cache_id::AbstractString;
    full_rebuild::Bool=false,
    should_cancel::Union{Nothing,Function}=nothing,
    on_progress::Union{Nothing,Function}=nothing,
)
    identity = project_cache_identity(cache_id, project, root_path)
    mkpath(dirname(identity.cache_path))
    @info "Starting project cache build" root = identity.root_path cache = identity.cache_path full_rebuild
    fingerprints = collect_csv_fingerprints(identity.root_path; should_cancel, on_progress)
    @info "Discovered CSV files for project cache" count = length(fingerprints) root = identity.root_path
    cache_exists = isfile(identity.cache_path)
    mode = (!cache_exists || full_rebuild) ? "w" : "r+"
    meta = _load_scan_metadata(identity.root_path)
    processed = 0
    loaded_measurements = 0

    hierarchy = h5open(identity.cache_path, mode) do h5
        if cache_exists && !full_rebuild
            _validate_meta!(h5, identity)
        else
            _write_meta!(h5, identity)
            _write_semantic_fields!(h5, project_cache_semantic_fields(project))
        end
        files_group = _ensure_group(h5, "files")
        cached = _cached_file_fingerprints(h5, identity.cache_path)
        statuses = _cached_file_statuses(h5)
        raw_paths = Set(keys(fingerprints))
        for cached_path in keys(cached)
            cached_path in raw_paths && continue
            key = _file_group_key(cached_path)
            haskey(files_group, key) && HDF5.delete_object(files_group, key)
        end

        for path in sort!(collect(keys(fingerprints)))
            _check_cancel(should_cancel)
            fp = fingerprints[path]
            if !full_rebuild &&
               haskey(cached, path) &&
               _same_fingerprint(fp, cached[path]) &&
               get(statuses, path, "ok") != "error"
                processed += 1
                on_progress !== nothing && on_progress((
                    phase=:cache_update,
                    total_csv=length(fingerprints),
                    processed_csv=processed,
                    loaded_measurements=loaded_measurements,
                    skipped_csv=0,
                    current_path=path,
                ))
                continue
            end
            indexed = index_csv_file(path)
            try
                loaded_measurements += _write_file_group!(
                    files_group,
                    project,
                    fp,
                    indexed,
                    meta;
                    should_cancel,
                )
            catch err
                (err isa ScanCancelled || err isa PlotCancelled) && rethrow()
                _write_file_error_group!(files_group, fp, indexed, err)
            end
            processed += 1
            on_progress !== nothing && on_progress((
                phase=:cache_update,
                total_csv=length(fingerprints),
                processed_csv=processed,
                loaded_measurements=loaded_measurements,
                skipped_csv=0,
                current_path=path,
            ))
        end
        _write_meta!(h5, identity)
        _write_semantic_fields!(h5, project_cache_semantic_fields(project))
        _rebuild_indexes!(h5, identity, project)
    end
    status = _cache_status_from_fingerprints(identity, fingerprints)
    errors = _cache_file_errors(identity)
    @info "Finished project cache build" cache = identity.cache_path cached_files = status.cached_files error_files = status.error_files
    return ProjectCacheSnapshot(
        identity,
        hierarchy,
        status,
        project_cache_semantic_fields(project),
        errors,
    )
end

function _load_project_cache_contents(
    root_path::AbstractString,
    project::AbstractProject,
    cache_id::AbstractString,
    ;
    should_cancel::Union{Nothing,Function}=nothing,
    on_progress::Union{Nothing,Function}=nothing,
    on_file_loaded::Union{Nothing,Function}=nothing,
)
    measurements = MeasurementInfo[]
    metadata = _stream_project_cache_contents(
        root_path,
        project,
        cache_id;
        should_cancel,
        on_progress,
        on_file_loaded=(file_measurements) -> begin
            on_file_loaded !== nothing && on_file_loaded(file_measurements)
            append!(measurements, file_measurements)
        end,
    )
    hierarchy = MeasurementHierarchy(
        measurements,
        metadata.identity.root_path,
        metadata.has_device_metadata,
        project,
        metadata.skipped_count,
    )
    return ProjectCacheSnapshot(
        metadata.identity,
        hierarchy,
        metadata.status,
        metadata.semantic_fields,
        metadata.errors,
    )
end

function _stream_project_cache_contents(
    root_path::AbstractString,
    project::AbstractProject,
    cache_id::AbstractString,
    ;
    should_cancel::Union{Nothing,Function}=nothing,
    on_progress::Union{Nothing,Function}=nothing,
    on_file_loaded::Union{Nothing,Function}=nothing,
)
    identity = project_cache_identity(cache_id, project, root_path)
    isfile(identity.cache_path) || throw(ProjectCacheMissingError(identity.cache_path))
    _check_cancel(should_cancel)
    metadata = h5open(identity.cache_path, "r") do h5
        _validate_meta!(h5, identity)
        files_group = haskey(h5, "files") ? h5["files"] :
            throw(ProjectCacheInvalidError(identity.cache_path, "missing /files group"))
        file_keys = sort!(collect(keys(files_group)))
        skipped_count = 0
        loaded_measurements = 0
        processed = 0
        semantic_fields = _read_semantic_fields(h5, identity.cache_path)
        has_device_metadata = _read_cache_has_device_metadata(h5, identity.cache_path)
        cached_status = _cache_status_from_cached_statuses(_cached_file_statuses(h5))
        for file_key in file_keys
            _check_cancel(should_cancel)
            file_group = files_group[file_key]
            status = _file_status(file_group)
            status in ("skipped", "error") && (skipped_count += 1)
            file_measurements = _read_file_measurements(file_group, identity.cache_path)
            if on_file_loaded !== nothing && !isempty(file_measurements)
                on_file_loaded(file_measurements)
            end
            loaded_measurements += length(file_measurements)
            processed += 1
            current_path = haskey(file_group, "path") ? read(file_group["path"]) : file_key
            _emit_progress(on_progress;
                phase=:cache_load,
                total_csv=length(file_keys),
                processed_csv=processed,
                loaded_measurements,
                skipped_csv=skipped_count,
                current_path,
            )
        end
        (
            identity=identity,
            status=cached_status,
            semantic_fields=semantic_fields,
            skipped_count=skipped_count,
            has_device_metadata=has_device_metadata,
        )
    end
    _check_cancel(should_cancel)
    errors = _cache_file_errors(identity)
    return (
        identity=metadata.identity,
        status=metadata.status,
        semantic_fields=metadata.semantic_fields,
        errors=errors,
        skipped_count=metadata.skipped_count,
        has_device_metadata=metadata.has_device_metadata,
    )
end

function load_project_cache(
    root_path::AbstractString,
    project::AbstractProject,
    cache_id::AbstractString,
    ;
    should_cancel::Union{Nothing,Function}=nothing,
    on_progress::Union{Nothing,Function}=nothing,
)
    snapshot = _load_project_cache_contents(root_path, project, cache_id; should_cancel, on_progress)
    status = cache_status(snapshot.identity; should_cancel, on_progress)
    return ProjectCacheSnapshot(
        snapshot.identity,
        snapshot.hierarchy,
        status,
        snapshot.semantic_fields,
        snapshot.errors,
    )
end

function _measurement_group_for_cached_plot(
    identity::ProjectCacheIdentity,
    measurement::MeasurementInfo,
)
    isfile(identity.cache_path) || throw(ProjectCacheMissingError(identity.cache_path))
    return h5open(identity.cache_path, "r") do h5
        _validate_meta!(h5, identity)
        file_key = _file_group_key(measurement.filepath)
        measurement_key = _measurement_group_key(measurement.id)
        files_group = h5["files"]
        haskey(files_group, file_key) ||
            throw(ProjectCacheInvalidError(identity.cache_path, "missing cached file group for $(measurement.filepath)"))
        measurements_group = files_group[file_key]["measurements"]
        haskey(measurements_group, measurement_key) ||
            throw(ProjectCacheInvalidError(identity.cache_path, "missing cached measurement group for $(measurement.id)"))
        return project_cache_read_plot_payload(
            _project_by_name(identity.project_name),
            measurement,
            files_group[file_key],
            measurements_group[measurement_key],
        )
    end
end

function _project_by_name(name::AbstractString)
    for project in KNOWN_PROJECTS
        project_name(project) == name && return project
    end
    throw(ProjectCacheInvalidError("", "unknown project '$name'"))
end
