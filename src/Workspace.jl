module Workspace

using DataFrames: DataFrame

import ..Cache:
    CacheDB,
    ProjectCacheError,
    ProjectCacheIdentity,
    ProjectCacheIndex,
    ProjectCacheStatus,
    clear_cache_index!,
    close_cache_db!,
    finalize_scan!,
    load_cache_index,
    open_cache_db,
    persist_stats!,
    project_cache_id,
    project_cache_identity,
    read_item_payloads,
    reconcile_source_item!,
    write_item_payloads!,
    write_scan_identity!
import ..ItemIndex:
    DataItem,
    Hierarchy,
    ItemFailure,
    ItemRecord,
    SourceScan,
    apply_collection_parameters!,
    check_cancel,
    insert_item!,
    is_job_cancelled,
    metadata_dict,
    scan_source,
    with_cancel
import ..Projects:
    AbstractDataSource,
    AbstractDataItem,
    Project,
    close_source!,
    collection_stats,
    compute_item_stats,
    cacheable,
    id,
    item_data,
    load_data_item,
    project_name,
    source_id,
    source_label

"""
Progress reported by one workspace job.
"""
Base.@kwdef mutable struct WorkspaceProgress
    phase::Symbol = :idle
    total_source_items::Int = 0
    processed_source_items::Int = 0
    loaded_items::Int = 0
    skipped_source_items::Int = 0
    current_source_item::String = ""
end

"""
Convert one scan or cache progress event into workspace-owned progress.
"""
function WorkspaceProgress(progress::NamedTuple)::WorkspaceProgress
    return WorkspaceProgress(
        phase=progress.phase,
        total_source_items=progress.total_source_items,
        processed_source_items=progress.processed_source_items,
        loaded_items=progress.loaded_items,
        skipped_source_items=progress.skipped_source_items,
        current_source_item=progress.current_source_item,
    )
end

"""
One cancellable workspace operation and its latest visible state.
"""
mutable struct WorkspaceJob
    id::Int
    state::Symbol
    progress::WorkspaceProgress
    error::String
    events::Union{Nothing,Channel{NamedTuple}}
    cancel_token::Union{Nothing,Base.Threads.Atomic{Bool}}
end

WorkspaceJob()::WorkspaceJob =
    WorkspaceJob(0, :idle, WorkspaceProgress(), "", nothing, nothing)

"""
The progressively populated item index for one open source.
"""
mutable struct WorkspaceIndex
    hierarchy::Hierarchy
    items::Dict{String,ItemRecord}
    collection_parameter_keys::Vector{Symbol}
    source::Union{Nothing,SourceScan}
    analysis_errors::Dict{String,String}
end

"""
Stable selection identities owned by a workspace.
"""
Base.@kwdef mutable struct WorkspaceSelection
    collection_paths::Vector{String} = String[]
    item_ids::Vector{String} = String[]
end

"""
Loaded cache state for one workspace.
"""
mutable struct WorkspaceCache
    identity::ProjectCacheIdentity
    db::CacheDB
    index::Union{Nothing,ProjectCacheIndex}
    status::Union{Nothing,ProjectCacheStatus}
    operation::Symbol
end

"""How many materialized items one workspace keeps resident before evicting the oldest."""
const DEFAULT_ITEM_CACHE_CAPACITY = 256

"""
Bounded least-recently-used cache of materialized items, keyed by item id.

Each entry stores the `(source_item_fingerprint, item_fingerprint)` it was loaded at, so a stale entry
is rejected on read. Eviction is cheap: an evicted item is re-materialized from the DuckDB payload
cache without re-reading the origin.
"""
mutable struct ItemCache
    capacity::Int
    entries::Dict{String,Tuple{Any,Any}}
    order::Vector{String}
end

ItemCache(capacity::Integer=DEFAULT_ITEM_CACHE_CAPACITY)::ItemCache =
    ItemCache(Int(capacity), Dict{String,Tuple{Any,Any}}(), String[])

"""Mark `id` most-recently-used."""
function _touch_item!(cache::ItemCache, id::String)::Nothing
    index = findfirst(==(id), cache.order)
    index === nothing || deleteat!(cache.order, index)
    push!(cache.order, id)
    return nothing
end

"""Return the cached `(fingerprint, item)` for `id`, marking it used, or `default` when absent."""
function lookup_item(cache::ItemCache, id::String, default)
    entry = get(cache.entries, id, nothing)
    entry === nothing && return default
    _touch_item!(cache, id)
    return entry
end

"""Insert or refresh one cached item, evicting the least-recently-used entries past capacity."""
function cache_item!(cache::ItemCache, id::String, fingerprint, item)
    if !haskey(cache.entries, id)
        while length(cache.order) >= cache.capacity && !isempty(cache.order)
            delete!(cache.entries, popfirst!(cache.order))
        end
    end
    cache.entries[id] = (fingerprint, item)
    _touch_item!(cache, id)
    return item
end

"""
One open project/source pair and all package-managed state belonging to it.
"""
mutable struct Workspace{S<:AbstractDataSource}
    project::Project
    source::S
    index::WorkspaceIndex
    selection::WorkspaceSelection
    cache::WorkspaceCache
    scan::WorkspaceJob
    analysis::WorkspaceJob
    cache_job::WorkspaceJob
    loaded_items::ItemCache
    background_tasks::Vector{Task}
    closed::Bool
end

"""
Create the empty state for one project-owned source.
"""
function Workspace(
    project::Project,
    source::S,
)::Workspace{S} where {S<:AbstractDataSource}
    hierarchy = Hierarchy(source_id(source), false)
    identity = project_cache_identity(project_cache_id(source), source)
    cache_db = open_cache_db(identity)
    return Workspace(
        project,
        source,
        WorkspaceIndex(
            hierarchy,
            Dict{String,ItemRecord}(),
            Symbol[],
            nothing,
            Dict{String,String}(),
        ),
        WorkspaceSelection(),
        WorkspaceCache(identity, cache_db, nothing, nothing, :load),
        WorkspaceJob(),
        WorkspaceJob(),
        WorkspaceJob(),
        ItemCache(),
        Task[],
        false,
    )
end

# ---------------------------------------------------------------------------
# Pretty printing
# ---------------------------------------------------------------------------

"""One-line description of where the source scan currently stands."""
function _scan_summary(workspace::Workspace)::String
    state = workspace.scan.state
    progress = workspace.scan.progress
    if state in (:counting, :discovering)
        return "$state ($(progress.processed_source_items) source items found)"
    elseif state == :scanning
        return "$state ($(progress.processed_source_items)/$(progress.total_source_items) source items)"
    else
        return string(state)
    end
end

Base.show(io::IO, workspace::Workspace) = print(
    io,
    "Workspace(", project_name(workspace.project), ", ", source_label(workspace.source),
    ", ", length(workspace.index.items), " items)",
)

function Base.show(io::IO, ::MIME"text/plain", workspace::Workspace)
    println(io, "Workspace · ", project_name(workspace.project))
    println(io, "  source label: ", source_label(workspace.source))
    println(io, "  source:       ", source_id(workspace.source))
    println(io, "  scan:         ", _scan_summary(workspace))
    println(io, "  cache:        ", workspace.cache_job.state)
    println(io, "  items:        ", length(workspace.index.items))
    print(io,   "  failures:     ", length(workspace.index.analysis_errors))
    workspace.closed && print(io, "\n  (closed)")
end

include("Workspace/Operations.jl")
include("Workspace/DataAccess.jl")

end
