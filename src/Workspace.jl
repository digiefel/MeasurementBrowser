module Workspace

using DataFrames: DataFrame

import ..Cache:
    ProjectCacheError,
    ProjectCacheIdentity,
    ProjectCacheIndex,
    ProjectCacheStatus,
    cache_status,
    load_project_cache,
    project_cache_id,
    project_cache_identity,
    write_project_cache!
import ..ItemIndex:
    Hierarchy,
    ItemRecord,
    SourceScan,
    check_cancel,
    insert_item!,
    item_record_key,
    is_job_cancelled,
    scan_source,
    with_cancel
import ..Projects:
    AbstractDataSource,
    AbstractDataItem,
    Project,
    RegisteredProjectSource,
    close_source!,
    item_data,
    item_id,
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
    collection_metadata_keys::Vector{Symbol}
    source::Union{Nothing,SourceScan}
    analysis_errors::Dict{String,String}
end

"""
Stable selection identities owned by a workspace.
"""
Base.@kwdef mutable struct WorkspaceSelection
    collection_paths::Vector{String} = String[]
    item_keys::Vector{String} = String[]
end

"""
Loaded cache state for one workspace.
"""
mutable struct WorkspaceCache
    identity::ProjectCacheIdentity
    index::Union{Nothing,ProjectCacheIndex}
    status::Union{Nothing,ProjectCacheStatus}
    operation::Symbol
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
    cache_job::WorkspaceJob
    loaded_items::Dict{String,Tuple{Any,Any}}
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
        WorkspaceCache(identity, nothing, nothing, :load),
        WorkspaceJob(),
        WorkspaceJob(),
        Dict{String,Tuple{Any,Any}}(),
        Task[],
        false,
    )
end

Workspace(project::Project, root_path::AbstractString) =
    Workspace(project, RegisteredProjectSource(project, root_path))

# ---------------------------------------------------------------------------
# Pretty printing
# ---------------------------------------------------------------------------

"""One-line description of where the source scan currently stands."""
function _scan_summary(workspace::Workspace)::String
    state = workspace.scan.state
    progress = workspace.scan.progress
    if state in (:counting, :discovering)
        return "$state ($(progress.processed_source_items) source items found)"
    elseif state in (:scanning, :analyzing)
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
