module Workspace

using DataFrames: DataFrame
using DBInterface
using Printf
import ..Profiling

import ..Cache:
    CacheDB,
    ProjectCacheIdentity,
    ProjectCacheIndex,
    ProjectCacheStatus,
    cached_item_data_ids,
    clear_cache_index!,
    close_cache_db!,
    finalize_scan!,
    load_cache_index,
    open_cache_db,
    persist_item_failure!,
    persist_stats!,
    project_cache_identity,
    read_cached_item_data,
    reconcile_source_items!,
    set_cache_memory_limit!,
    with_reader,
    with_writer_transaction,
    write_cached_item_data!,
    write_scan_identity!
import ..ItemIndex:
    DataItem,
    Hierarchy,
    ItemFailure,
    ItemRecord,
    SourceItemInterpretation,
    SourceScan,
    apply_collection_parameters!,
    check_cancel,
    insert_item!,
    is_job_cancelled,
    metadata_dict,
    interpret_source_item,
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
    process,
    project_name,
    source_id,
    source_item_id,
    source_items,
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
    item_stats::Dict{String,Dict{Symbol,Any}}
    collection_parameter_keys::Vector{Symbol}
    source::Union{Nothing,SourceScan}
    analysis_errors::Dict{String,String}
end

"""One deduplicated processing request and the selected callers waiting for its result."""
mutable struct ProcessingJob
    record::ItemRecord
    interpreted::Union{Nothing,AbstractDataItem}
    state::Symbol
    priority::Int
    waiters::Vector{Channel{Any}}
end

"""One completed processed value waiting for its small shared cache batch."""
struct ProcessedWriteRequest
    record::ItemRecord
    item::Union{Nothing,AbstractDataItem}
    completion::Channel{Any}
end

"""Workspace-owned processing work; completed values are not retained after delivery."""
mutable struct ProcessingQueue
    lock::ReentrantLock
    condition::Base.Threads.Condition
    jobs::Dict{String,ProcessingJob}
    selected::Vector{String}
    background::Vector{String}
    background_index::Int
    source_locks::Dict{String,ReentrantLock}
    events::Channel{NamedTuple}
    workers::Vector{Task}
    pending_writes::Vector{ProcessedWriteRequest}
    limit::Int
    total::Int
    completed::Int
    closed::Bool
end

function ProcessingQueue()::ProcessingQueue
    queue_lock = ReentrantLock()
    return ProcessingQueue(
        queue_lock,
        Base.Threads.Condition(queue_lock),
        Dict{String,ProcessingJob}(),
        String[],
        String[],
        1,
        Dict{String,ReentrantLock}(),
        Channel{NamedTuple}(Inf),
        Task[],
        ProcessedWriteRequest[],
        max(4, 2 * Base.Threads.nthreads()),
        0,
        0,
        false,
    )
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

"""
A single snapshot of everything a watcher needs to show about a workspace's background work.

This is the stable contract between the engine and any watcher (the GUI today, scripts and workflows
later): watchers read `WorkspaceStatus` and nothing else about jobs, progress, or the cache. It is
recomputed only when [`poll_workspace!`](@ref) observes a change or work is active, so an idle render
loop reads a cached value instead of rebuilding strings every frame.

- `level` drives the watcher's color/emphasis: `:none`, `:busy`, `:fresh`, `:stale`, `:missing`,
  `:error`.
- `label` is a short word for a button or chip ("Building", "Fresh", "Errors").
- `detail` is the one merged human line — the live activity while `busy`, otherwise a cache summary.
- `busy` is true while any scan, analysis, or cache work runs.
- `progress` is a determinate fraction when counts are known, or `nothing` for an indeterminate or
  absent bar.
- `errors` lists source-item failures as `id => first message line`, streamed as they occur.
"""
struct WorkspaceStatus
    level::Symbol
    label::String
    detail::String
    busy::Bool
    progress::Union{Nothing,Float32}
    errors::Vector{Pair{String,String}}
end

WorkspaceStatus() =
    WorkspaceStatus(:none, "Opening", "Opening the source…", true, nothing, Pair{String,String}[])

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
    # The cache is the persistence side of the scan, not a separately cancellable job, so it needs
    # only a visible state and last error rather than a full WorkspaceJob.
    cache_state::Symbol
    cache_error::String
    sampling_profile::Union{Nothing,Profiling.SamplingProfile}
    sampling_active::Bool
    processing::ProcessingQueue
    background_tasks::Vector{Task}
    status::WorkspaceStatus
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
    identity = project_cache_identity(project_name(project), source)
    cache_db = open_cache_db(identity)
    workspace = Workspace(
        project,
        source,
        WorkspaceIndex(
            hierarchy,
            Dict{String,ItemRecord}(),
            Dict{String,Dict{Symbol,Any}}(),
            Symbol[],
            nothing,
            Dict{String,String}(),
        ),
        WorkspaceSelection(),
        WorkspaceCache(identity, cache_db, nothing, nothing, :load),
        WorkspaceJob(),
        WorkspaceJob(),
        :idle,
        "",
        nothing,
        false,
        ProcessingQueue(),
        Task[],
        WorkspaceStatus(),
        false,
    )
    start_processing_workers!(workspace)
    return workspace
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
    println(io, "  cache:        ", workspace.cache_state)
    println(io, "  items:        ", length(workspace.index.items))
    print(io,   "  failures:     ", length(workspace.index.analysis_errors))
    workspace.closed && print(io, "\n  (closed)")
end

include("Workspace/Operations.jl")
include("Workspace/Status.jl")
include("Workspace/DataAccess.jl")
include("Workspace/Processing.jl")

end
