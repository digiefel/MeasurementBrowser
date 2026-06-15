module Workspace

using DataFrames: DataFrame

import ..Cache:
    ProjectCacheError,
    ProjectCacheIdentity,
    ProjectCacheIndex,
    ProjectCacheStatus,
    cache_status,
    cached_measurement_data,
    load_project_cache,
    project_cache_id,
    project_cache_identity,
    write_measurement_data_cache!,
    write_project_cache!
import ..MeasurementIndex:
    FileFingerprint,
    MeasurementHierarchy,
    MeasurementInfo,
    SourceFile,
    SourceScan,
    check_cancel,
    file_fingerprint,
    has_device_metadata,
    index_source_file,
    insert_measurement!,
    is_job_cancelled,
    scan_source,
    with_cancel
import ..Projects:
    AbstractProject,
    load_source_data,
    process_measurement_data,
    project_name

"""
Progress reported by one workspace job.
"""
Base.@kwdef mutable struct WorkspaceProgress
    phase::Symbol = :idle
    total_files::Int = 0
    processed_files::Int = 0
    loaded_measurements::Int = 0
    skipped_files::Int = 0
    current_path::String = ""
end

"""
Convert one scan or cache progress event into workspace-owned progress.
"""
function WorkspaceProgress(progress::NamedTuple)::WorkspaceProgress
    return WorkspaceProgress(
        phase=progress.phase,
        total_files=progress.total_csv,
        processed_files=progress.processed_csv,
        loaded_measurements=progress.loaded_measurements,
        skipped_files=progress.skipped_csv,
        current_path=progress.current_path,
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
The progressively populated measurement index for one open source root.
"""
mutable struct WorkspaceIndex
    hierarchy::MeasurementHierarchy
    measurements::Dict{String,MeasurementInfo}
    device_metadata_keys::Vector{Symbol}
    source::Union{Nothing,SourceScan}
    analysis_errors::Dict{String,String}
end

"""
Stable selection identities owned by a workspace.
"""
Base.@kwdef mutable struct WorkspaceSelection
    device_paths::Vector{String} = String[]
    measurement_ids::Vector{String} = String[]
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
One open measurement source root and all package-managed state belonging to it.
"""
mutable struct Workspace{P<:AbstractProject}
    project::P
    root_path::String
    index::WorkspaceIndex
    selection::WorkspaceSelection
    cache::WorkspaceCache
    scan::WorkspaceJob
    cache_job::WorkspaceJob
    direct_data::Dict{String,Tuple{FileFingerprint,DataFrame}}
    processed_data::Dict{String,Tuple{FileFingerprint,DataFrame}}
    background_tasks::Vector{Task}
    closed::Bool
end

"""
Create the empty state for one project and normalized source root.
"""
function Workspace(
    project::P,
    root_path::AbstractString,
)::Workspace{P} where {P<:AbstractProject}
    root = normpath(abspath(expanduser(String(root_path))))
    hierarchy = MeasurementHierarchy(root, has_device_metadata(root), project)
    identity = project_cache_identity(project_cache_id(root), project, root)
    return Workspace(
        project,
        root,
        WorkspaceIndex(
            hierarchy,
            Dict{String,MeasurementInfo}(),
            Symbol[],
            nothing,
            Dict{String,String}(),
        ),
        WorkspaceSelection(),
        WorkspaceCache(identity, nothing, nothing, :load),
        WorkspaceJob(),
        WorkspaceJob(),
        Dict{String,Tuple{FileFingerprint,DataFrame}}(),
        Dict{String,Tuple{FileFingerprint,DataFrame}}(),
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
        return "$state ($(progress.processed_files) files found)"
    elseif state in (:scanning, :analyzing)
        return "$state ($(progress.processed_files)/$(progress.total_files) files)"
    else
        return string(state)
    end
end

Base.show(io::IO, workspace::Workspace) = print(
    io,
    "Workspace(", project_name(workspace.project),
    ", ", length(workspace.index.measurements), " measurements)",
)

function Base.show(io::IO, ::MIME"text/plain", workspace::Workspace)
    println(io, "Workspace · ", project_name(workspace.project))
    println(io, "  root:         ", Base.contractuser(workspace.root_path))
    println(io, "  scan:         ", _scan_summary(workspace))
    println(io, "  cache:        ", workspace.cache_job.state)
    println(io, "  measurements: ", length(workspace.index.measurements))
    print(io,   "  failures:     ", length(workspace.index.analysis_errors))
    workspace.closed && print(io, "\n  (closed)")
end

include("Workspace/Operations.jl")
include("Workspace/DataAccess.jl")

end
