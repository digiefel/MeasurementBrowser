using HDF5
using SHA
using Serialization
using DataFrames: DataFrame

const PROJECT_CACHE_COMPRESSION = 3
const PROJECT_CACHE_SCHEMA_VERSION = 4
const PROJECT_CACHE_INDEX_DATASET = "index"
const PROJECT_CACHE_LOCK = ReentrantLock()

"""Error raised when a project cache is missing or does not match the current project."""
struct ProjectCacheError <: Exception
    path::String
    message::String
end

"""Print the cache path and validation failure."""
Base.showerror(io::IO, err::ProjectCacheError)::Nothing =
    print(io, "Invalid project cache $(err.path): $(err.message)")

"""The project, source root, and HDF5 file belonging to one cache."""
struct ProjectCacheIdentity
    cache_id::String
    project_name::String
    root_path::String
    cache_path::String
end

"""Counts describing the differences between a source scan and its cache."""
struct ProjectCacheStatus
    total_files::Int
    cached_files::Int
    fresh_files::Int
    stale_files::Int
    new_files::Int
    deleted_files::Int
    error_files::Int
end

"""All non-dataframe content needed to restore and compare a project cache."""
struct ProjectCacheIndex
    identity::ProjectCacheIdentity
    source::SourceScan
    files::Dict{String,SourceFile}
    analysis_errors::Dict{String,String}
end

const ACTIVE_PROJECT_CACHE = Ref{Union{Nothing,ProjectCacheIndex}}(nothing)

"""
Return the deterministic cache id for a source root.

The folder name keeps paths recognizable; the digest distinguishes equal folder names at different
locations.
"""
function project_cache_id(root_path::AbstractString)::String
    root = normpath(abspath(expanduser(String(root_path))))
    name = replace(lowercase(basename(root)), r"[^a-z0-9_.-]+" => "_")
    isempty(name) && (name = "project")
    return "$name-$(bytes2hex(sha1(root))[1:12])"
end

"""
Bind one cache id to a project and normalized source root.

The returned identity also contains the package-owned HDF5 path.
"""
function project_cache_identity(
    cache_id::AbstractString,
    project::AbstractProject,
    root_path::AbstractString,
)::ProjectCacheIdentity
    root = normpath(abspath(expanduser(String(root_path))))
    depot = isempty(DEPOT_PATH) ? homedir() : DEPOT_PATH[1]
    project_dir = replace(lowercase(project_name(project)), r"[^a-z0-9_.-]+" => "_")
    return ProjectCacheIdentity(
        String(cache_id),
        project_name(project),
        root,
        joinpath(
            depot,
            "measurementbrowser",
            "cache",
            project_dir,
            "$(String(cache_id)).h5",
        ),
    )
end

"""
Return a fixed-length HDF5 name for a source path or measurement id.
"""
function cache_object_key(value::AbstractString)::String
    return bytes2hex(sha1(String(value)))
end

"""
Replace one HDF5 dataset with a Julia serialization of `value`.

Cache indexes are contiguous for fast startup reads. Measurement payloads use larger compressed
chunks to reduce disk use without fragmenting them into many tiny reads.
"""
function write_serialized_dataset!(
    parent::Union{HDF5.File,HDF5.Group},
    name::AbstractString,
    value::Any,
    ;
    compress::Bool=true,
)::Nothing
    io = IOBuffer()
    serialize(io, value)
    bytes = take!(io)
    haskey(parent, name) && HDF5.delete_object(parent, name)
    if compress
        write(
            parent,
            name,
            bytes;
            chunk=(min(length(bytes), 262_144),),
            deflate=PROJECT_CACHE_COMPRESSION,
        )
    else
        write(parent, name, bytes)
    end
    return nothing
end

"""
Read and deserialize one required HDF5 dataset.

Missing datasets identify an invalid cache rather than a cache miss.
"""
function read_serialized_dataset(
    parent::Union{HDF5.File,HDF5.Group},
    name::AbstractString,
    cache_path::AbstractString,
)::Any
    haskey(parent, name) ||
        throw(ProjectCacheError(String(cache_path), "missing dataset '$name'"))
    return deserialize(IOBuffer(read(parent[name])))
end

"""
Open one HDF5 cache while holding the process-wide cache lock.

All cache reads and writes use this path so background updates cannot overlap plotting reads.
"""
function with_project_cache_file(
    read_or_write::Function,
    identity::ProjectCacheIdentity,
    mode::AbstractString,
)::Any
    return lock(PROJECT_CACHE_LOCK) do
        h5open(read_or_write, identity.cache_path, mode)
    end
end

"""
Open and validate one existing project cache index.

Validation happens before deserialization so an old schema is reported clearly.
"""
function project_cache_index(identity::ProjectCacheIdentity)::ProjectCacheIndex
    isfile(identity.cache_path) ||
        throw(ProjectCacheError(identity.cache_path, "file does not exist"))
    return with_project_cache_file(identity, "r") do h5
        file_attributes = HDF5.attributes(h5)
        haskey(file_attributes, "schema_version") ||
            throw(ProjectCacheError(identity.cache_path, "cache schema is out of date"))
        read(file_attributes["schema_version"]) == PROJECT_CACHE_SCHEMA_VERSION ||
            throw(ProjectCacheError(identity.cache_path, "cache schema is out of date"))
        index = read_serialized_dataset(h5, PROJECT_CACHE_INDEX_DATASET, identity.cache_path)
        index isa ProjectCacheIndex ||
            throw(ProjectCacheError(identity.cache_path, "cache index has the wrong type"))
        index.identity.cache_id == identity.cache_id &&
            index.identity.project_name == identity.project_name &&
            index.identity.root_path == identity.root_path ||
            throw(ProjectCacheError(identity.cache_path, "cache identity does not match"))
        return index
    end
end

"""
Build the complete cache index represented by one finished source scan.

Analysis failures remain attached to their physical source files while successful measurements stay
available.
"""
function ProjectCacheIndex(
    identity::ProjectCacheIdentity,
    source::SourceScan,
)::ProjectCacheIndex
    errors = Dict{String,Vector{String}}()
    for failure in source.analysis_failures
        push!(get!(errors, failure.filepath, String[]), failure.message)
    end
    return ProjectCacheIndex(
        identity,
        source,
        Dict(file.fingerprint.path => file for file in source.files),
        Dict(path => join(messages, "\n") for (path, messages) in errors),
    )
end

"""
Compare a cached index with the index produced by the current source scan.

The result counts matching, changed, new, deleted, and failed source files.
"""
function cache_status(
    cached::ProjectCacheIndex,
    source::SourceScan,
)::ProjectCacheStatus
    current = ProjectCacheIndex(cached.identity, source)
    stale = 0
    fresh = 0
    for (path, source_file) in current.files
        previous = get(cached.files, path, nothing)
        previous === nothing && continue
        changed = previous.fingerprint != source_file.fingerprint ||
            get(cached.analysis_errors, path, "") != get(current.analysis_errors, path, "")
        stale += changed
        fresh += !changed && !haskey(current.analysis_errors, path)
    end
    return ProjectCacheStatus(
        length(current.files),
        length(cached.files),
        fresh,
        stale,
        count(path -> !haskey(cached.files, path), keys(current.files)),
        count(path -> !haskey(current.files, path), keys(cached.files)),
        length(current.analysis_errors),
    )
end

"""
Make one loaded cache index available to cache-aware measurement-data access.
"""
function set_active_project_cache!(
    index::Union{Nothing,ProjectCacheIndex},
)::Nothing
    lock(PROJECT_CACHE_LOCK) do
        ACTIVE_PROJECT_CACHE[] = index
    end
    return nothing
end

"""
Create or update a cache from a finished source scan. `replace=true` discards all existing payloads.
Normal updates preserve payloads for unchanged files and remove only data made invalid by changed or
deleted files.
"""
function write_project_cache!(
    identity::ProjectCacheIdentity,
    source::SourceScan;
    replace::Bool=false,
    on_progress::Union{Nothing,Function}=nothing,
)::ProjectCacheIndex
    identity.root_path == source.root_path ||
        error("Cache root does not match source scan root")
    identity.project_name == project_name(source.project) ||
        error("Cache project does not match source scan project")

    new_index = ProjectCacheIndex(identity, source)
    cache_exists = isfile(identity.cache_path)
    old_index = cache_exists && !replace ? project_cache_index(identity) : nothing

    mkpath(dirname(identity.cache_path))
    with_project_cache_file(identity, !cache_exists || replace ? "w" : "r+") do h5
        if old_index !== nothing
            invalid_files = filter(collect(keys(old_index.files))) do path
                current = get(new_index.files, path, nothing)
                return current === nothing ||
                    old_index.files[path].fingerprint != current.fingerprint ||
                    get(old_index.analysis_errors, path, "") !=
                        get(new_index.analysis_errors, path, "")
            end
            for group_name in ("direct", "processed")
                haskey(h5, group_name) || continue
                group = h5[group_name]
                for path in invalid_files
                    key = cache_object_key(path)
                    haskey(group, key) && HDF5.delete_object(group, key)
                end
            end
        end

        check_cancel()
        loaded_measurements = sum(
            length(file.measurements) for file in values(new_index.files)
        )
        file_attributes = HDF5.attributes(h5)
        haskey(file_attributes, "schema_version") ||
            (file_attributes["schema_version"] = PROJECT_CACHE_SCHEMA_VERSION)
        write_serialized_dataset!(
            h5,
            PROJECT_CACHE_INDEX_DATASET,
            new_index;
            compress=false,
        )
        emit_progress(
            on_progress;
            phase=:cache_finalize,
            total_csv=1,
            processed_csv=1,
            loaded_measurements,
            skipped_csv=0,
            current_path=identity.cache_path,
        )
    end

    return new_index
end

"""
Load the complete browser index without scanning or opening any source file.
"""
function load_project_cache(
    identity::ProjectCacheIdentity;
    on_progress::Union{Nothing,Function}=nothing,
)::ProjectCacheIndex
    index = project_cache_index(identity)
    emit_progress(
        on_progress;
        phase=:cache_load,
        total_csv=1,
        processed_csv=1,
        loaded_measurements=length(index.source.hierarchy.all_measurements),
        skipped_csv=index.source.hierarchy.skipped_count,
        current_path=identity.cache_path,
    )
    return index
end

"""
Return the cached source entry when the file still exists and is unchanged. Reuse `fingerprints`
when several measurements belong to the same physical file.
"""
function valid_cached_source_file(
    index::ProjectCacheIndex,
    measurement::MeasurementInfo,
    fingerprints::Dict{String,FileFingerprint},
)::Union{Nothing,SourceFile}
    path = normpath(abspath(expanduser(measurement.filepath)))
    isfile(path) || return nothing
    source_file = get(index.files, path, nothing)
    source_file === nothing && return nothing
    current = get!(fingerprints, path) do
        file_fingerprint(path)
    end
    current == source_file.fingerprint || return nothing
    return source_file
end

"""
Read every valid requested payload in one HDF5 operation. Missing or stale entries return `nothing`.
"""
function cached_measurement_data(
    project::AbstractProject,
    measurements::Vector{MeasurementInfo};
    processed::Bool=false,
)::Vector{Union{Nothing,DataFrame}}
    data = Union{Nothing,DataFrame}[nothing for _ in measurements]
    index = ACTIVE_PROJECT_CACHE[]
    index === nothing && return data
    identity = index.identity
    identity.project_name == project_name(project) || return data
    isfile(identity.cache_path) || return data

    fingerprints = Dict{String,FileFingerprint}()
    group_name = processed ? "processed" : "direct"
    with_project_cache_file(identity, "r") do h5
        haskey(h5, group_name) || return data
        group = h5[group_name]
        for (position, measurement) in pairs(measurements)
            check_cancel()
            valid_cached_source_file(index, measurement, fingerprints) === nothing && continue
            file_key = cache_object_key(measurement.filepath)
            haskey(group, file_key) || continue
            file_group = group[file_key]
            data_key = cache_object_key(measurement.unique_id)
            haskey(file_group, data_key) || continue
            value = read_serialized_dataset(
                file_group,
                data_key,
                identity.cache_path,
            )
            value isa DataFrame ||
                throw(ProjectCacheError(
                    identity.cache_path,
                    "measurement data has type $(typeof(value))",
                ))
            data[position] = value
        end
    end
    return data
end

"""
Write valid measurement payloads to the active cache in one HDF5 operation.
"""
function write_measurement_data_cache!(
    project::AbstractProject,
    measurements::Vector{MeasurementInfo},
    data::Vector{DataFrame};
    processed::Bool=false,
)::Nothing
    length(measurements) == length(data) ||
        throw(DimensionMismatch("measurements and data must have equal lengths"))
    isempty(measurements) && return nothing
    index = ACTIVE_PROJECT_CACHE[]
    index === nothing && return nothing
    identity = index.identity
    identity.project_name == project_name(project) || return nothing
    isfile(identity.cache_path) || return nothing

    fingerprints = Dict{String,FileFingerprint}()
    group_name = processed ? "processed" : "direct"
    with_project_cache_file(identity, "r+") do h5
        group = haskey(h5, group_name) ? h5[group_name] : create_group(h5, group_name)
        for (measurement, value) in zip(measurements, data)
            check_cancel()
            valid_cached_source_file(index, measurement, fingerprints) === nothing && continue
            file_key = cache_object_key(measurement.filepath)
            file_group = haskey(group, file_key) ? group[file_key] : create_group(group, file_key)
            write_serialized_dataset!(
                file_group,
                cache_object_key(measurement.unique_id),
                value,
            )
        end
    end
    return nothing
end
