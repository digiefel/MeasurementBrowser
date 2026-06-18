using HDF5
using SHA
using Serialization
const PROJECT_CACHE_COMPRESSION = 3
const PROJECT_CACHE_SCHEMA_VERSION = 11
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

"""The source and HDF5 file belonging to one cache."""
struct ProjectCacheIdentity
    cache_id::String
    source_id::String
    source_label::String
    cache_path::String
end

"""Counts describing the differences between a source scan and its cache."""
struct ProjectCacheStatus
    total_source_items::Int
    cached_source_items::Int
    fresh_source_items::Int
    stale_source_items::Int
    new_source_items::Int
    deleted_source_items::Int
    error_source_items::Int
end

"""All non-dataframe content needed to restore and compare a project cache."""
struct ProjectCacheIndex
    identity::ProjectCacheIdentity
    source::SourceScan
    analysis_errors::Dict{String,String}
end

"""
Return the deterministic cache id for a source.

The folder name keeps ids recognizable; the digest distinguishes equal labels at different origins.
"""
function project_cache_id(source)::String
    id = source_id(source)
    name = replace(lowercase(source_label(source)), r"[^a-z0-9_.-]+" => "_")
    isempty(name) && (name = "source")
    return "$name-$(bytes2hex(sha1(id))[1:12])"
end

"""
Bind one cache id to a source.

The returned identity also contains the package-owned HDF5 path.
"""
function project_cache_identity(
    cache_id::AbstractString,
    source,
)::ProjectCacheIdentity
    depot = isempty(DEPOT_PATH) ? homedir() : DEPOT_PATH[1]
    source_dir = replace(lowercase(source_label(source)), r"[^a-z0-9_.-]+" => "_")
    return ProjectCacheIdentity(
        String(cache_id),
        source_id(source),
        source_label(source),
        joinpath(
            depot,
            "measurementbrowser",
            "cache",
            source_dir,
            "$(String(cache_id)).h5",
        ),
    )
end

"""
Return a fixed-length HDF5 name for a source-item or item key.
"""
function cache_object_key(value::AbstractString)::String
    return bytes2hex(sha1(String(value)))
end

"""
Replace one HDF5 dataset with a Julia serialization of `value`.

Cache indexes are contiguous for fast startup reads. Item payloads use larger compressed
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
    bytes = read(parent[name])
    # A cache written by an incompatible struct layout deserializes into garbage or throws. Treat any
    # such failure as an invalid (rebuildable) cache rather than a fatal error, so the browser
    # self-heals instead of refusing to load forever.
    try
        return deserialize(IOBuffer(bytes))
    catch err
        throw(ProjectCacheError(
            String(cache_path),
            "could not deserialize dataset '$name' (incompatible cache): " *
                sprint(showerror, err),
        ))
    end
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
            index.identity.source_id == identity.source_id &&
            index.identity.source_label == identity.source_label ||
            throw(ProjectCacheError(identity.cache_path, "cache identity does not match"))
        return index
    end
end

"""
Build the complete cache index represented by one finished source scan.

Analysis failures remain attached to their source items while successful items stay available.
"""
function ProjectCacheIndex(
    identity::ProjectCacheIdentity,
    source::SourceScan,
)::ProjectCacheIndex
    errors = Dict{String,Vector{String}}()
    for failure in source.analysis_failures
        push!(get!(errors, failure.source_item_id, String[]), failure.message)
    end
    return ProjectCacheIndex(
        identity,
        source,
        Dict(id => join(messages, "\n") for (id, messages) in errors),
    )
end

"""
Compare a cached index with the index produced by the current source scan.

The result counts matching, changed, new, deleted, and failed source items.
"""
function cache_status(
    cached::ProjectCacheIndex,
    source::SourceScan,
)::ProjectCacheStatus
    current = ProjectCacheIndex(cached.identity, source)
    stale = 0
    fresh = 0
    for (id, fingerprint) in source.source_item_fingerprints
        previous = get(cached.source.source_item_fingerprints, id, nothing)
        previous === nothing && continue
        changed = previous != fingerprint ||
            get(cached.analysis_errors, id, "") != get(current.analysis_errors, id, "")
        stale += changed
        fresh += !changed && !haskey(current.analysis_errors, id)
    end
    current_ids = keys(source.source_item_fingerprints)
    cached_ids = keys(cached.source.source_item_fingerprints)
    return ProjectCacheStatus(
        length(current_ids),
        length(cached_ids),
        fresh,
        stale,
        count(id -> !haskey(cached.source.source_item_fingerprints, id), current_ids),
        count(id -> !haskey(source.source_item_fingerprints, id), cached_ids),
        length(current.analysis_errors),
    )
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
    identity.source_id == source.source_id ||
        error("Cache source does not match source scan")

    new_index = ProjectCacheIndex(identity, source)
    cache_exists = isfile(identity.cache_path)
    old_index = cache_exists && !replace ? project_cache_index(identity) : nothing

    mkpath(dirname(identity.cache_path))
    with_project_cache_file(identity, !cache_exists || replace ? "w" : "r+") do h5
        if old_index !== nothing
            invalid_items = filter(collect(keys(old_index.source.source_item_fingerprints))) do id
                current = get(new_index.source.source_item_fingerprints, id, nothing)
                return current === nothing ||
                    old_index.source.source_item_fingerprints[id] != current ||
                    get(old_index.analysis_errors, id, "") !=
                        get(new_index.analysis_errors, id, "")
            end
            for group_name in ("items",)
                haskey(h5, group_name) || continue
                group = h5[group_name]
                for id in invalid_items
                    key = cache_object_key(id)
                    haskey(group, key) && HDF5.delete_object(group, key)
                end
            end
        end

        check_cancel()
        loaded_items = length(new_index.source.hierarchy.all_items)
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
            total_source_items=1,
            processed_source_items=1,
            loaded_items,
            skipped_source_items=0,
            current_source_item=identity.cache_path,
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
        total_source_items=1,
        processed_source_items=1,
        loaded_items=length(index.source.hierarchy.all_items),
        skipped_source_items=index.source.hierarchy.skipped_count,
        current_source_item=identity.cache_path,
    )
    return index
end

"""
Read every valid requested payload in one HDF5 operation. Missing or stale entries return `nothing`.
"""
function cached_item_data(
    index::Union{Nothing,ProjectCacheIndex},
    items::Vector{ItemRecord},
)::Vector{Any}
    data = Any[nothing for _ in items]
    index === nothing && return data
    identity = index.identity
    isfile(identity.cache_path) || return data

    group_name = "items"
    with_project_cache_file(identity, "r") do h5
        haskey(h5, group_name) || return data
        group = h5[group_name]
        for (position, item) in pairs(items)
            check_cancel()
            item.source_item_fingerprint === nothing && continue
            item.item_fingerprint === nothing && continue
            get(index.source.source_item_fingerprints, item.source_item_id, nothing) ==
                item.source_item_fingerprint || continue
            group_key = cache_object_key(item.source_item_id)
            haskey(group, group_key) || continue
            item_group = group[group_key]
            data_key = cache_object_key(item.id)
            haskey(item_group, data_key) || continue
            value = read_serialized_dataset(
                item_group,
                data_key,
                identity.cache_path,
            )
            data[position] = value
        end
    end
    return data
end

"""
Write valid item payloads to the specified loaded cache in one HDF5 operation.
"""
function write_item_data_cache!(
    index::Union{Nothing,ProjectCacheIndex},
    items::Vector{ItemRecord},
    data::Vector,
)::Nothing
    length(items) == length(data) ||
        throw(DimensionMismatch("items and data must have equal lengths"))
    isempty(items) && return nothing
    index === nothing && return nothing
    identity = index.identity
    isfile(identity.cache_path) || return nothing

    group_name = "items"
    with_project_cache_file(identity, "r+") do h5
        group = haskey(h5, group_name) ? h5[group_name] : create_group(h5, group_name)
        for (item, value) in zip(items, data)
            check_cancel()
            item.source_item_fingerprint === nothing && continue
            item.item_fingerprint === nothing && continue
            get(index.source.source_item_fingerprints, item.source_item_id, nothing) ==
                item.source_item_fingerprint || continue
            group_key = cache_object_key(item.source_item_id)
            item_group = haskey(group, group_key) ? group[group_key] : create_group(group, group_key)
            write_serialized_dataset!(
                item_group,
                cache_object_key(item.id),
                value,
            )
        end
    end
    return nothing
end
