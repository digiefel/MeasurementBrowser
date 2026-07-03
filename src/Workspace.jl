module Workspace

using Printf
import ..Profiling

import ..Cache:
    AbstractCacheDB,
    BuildMetrics,
    CacheResultKey,
    CacheResultKind,
    CacheResultStatus,
    CacheDB,
    COLLECTION_STATS_RESULT,
    ITEM_STATS_RESULT,
    PROCESSING_RESULT,
    ProjectCacheSchemaError,
    ProjectCacheIdentity,
    ProjectCacheIndex,
    ProjectCacheStatus,
    RESULT_FAILED,
    RESULT_READY,
    cache_built,
    cache_has_pending_writes,
    cache_pending_counts,
    clear_cache_index!,
    close_cache_db!,
    delete_collection_stats!,
    delete_source_item!,
    invalidate_item_results!,
    load_cache_index,
    open_memory_cache_db,
    open_cache_db,
    project_cache_identity,
    read_item_data,
    record_cache_phase!,
    result_input_fingerprint,
    reset_build_metrics!,
    set_cache_memory_limit!,
    start_cache!,
    stop_cache!,
    store_interpreted_data!,
    store_interpreted_records!,
    store_interpreted!,
    store_item_stats!,
    store_collection_stats!,
    store_processed!,
    store_result_failure!,
    wait_condition_deadline,
    write_meta_header!
import ..ItemIndex:
    DataItem,
    Hierarchy,
    ItemFailure,
    ItemRecord,
    JobCancelled,
    MetadataDict,
    SourceItemInterpretation,
    SourceScan,
    check_cancel,
    collection_path_key,
    collection_path_tuple,
    effective_parameters,
    effective_record,
    is_job_cancelled,
    metadata_dict,
    interpret_source_item,
    with_cancel
import ..Projects:
    AbstractDataSource,
    AbstractDataSourceItem,
    AbstractDataItem,
    Project,
    SourceChanges,
    SourceError,
    close_source!,
    collection_stats,
    compute_item_stats,
    cacheable,
    fingerprint,
    has_collection_parameters,
    id,
    item_data,
    process,
    open_source,
    project_name,
    record_scan_phase!,
    reset_scan_profile!,
    scan_profile_summary,
    source_id,
    source_item_id,
    source_items,
    source_label,
    watch_source

"""
One cancellable workspace operation and its latest visible state.
"""
mutable struct WorkspaceJob
    id::Int
    state::Symbol
    error::String
    cancel_token::Union{Nothing,Base.Threads.Atomic{Bool}}
end

WorkspaceJob()::WorkspaceJob = WorkspaceJob(0, :idle, "", nothing)

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

@enum WorkKind begin
    INTERPRET_SOURCE
    PROCESS_ITEM
    ITEM_STATS
    COLLECTION_STATS
end

"""Stable identity of one result-producing operation."""
struct WorkKey
    kind::WorkKind
    entity::String
end

"""One revision of a result-producing operation."""
mutable struct WorkNode
    key::WorkKey
    revision::UInt64
    state::Symbol
    priority::Int
    dependencies::Vector{WorkKey}
    waiters::Vector{Channel{Any}}
    queued_ns::UInt64
end

"""Workspace-owned dependency state and scheduling queue."""
mutable struct WorkDependencyGraph
    lock::ReentrantLock
    condition::Base.Threads.Condition
    nodes::Dict{WorkKey,WorkNode}
    dependents::Dict{WorkKey,Set{WorkKey}}
    queue::Vector{Tuple{WorkKey,UInt64}}
    source_items::Dict{String,AbstractDataSourceItem}
    source_locks::Dict{String,ReentrantLock}
    workers::Vector{Task}
    total::Int
    completed::Int
    source_batch::Int
    source_batch_open::Bool
    closed::Bool
end

function WorkDependencyGraph()::WorkDependencyGraph
    work_lock = ReentrantLock()
    return WorkDependencyGraph(
        work_lock,
        Base.Threads.Condition(work_lock),
        Dict{WorkKey,WorkNode}(),
        Dict{WorkKey,Set{WorkKey}}(),
        Tuple{WorkKey,UInt64}[],
        Dict{String,AbstractDataSourceItem}(),
        Dict{String,ReentrantLock}(),
        Task[],
        0,
        0,
        0,
        false,
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
    db::AbstractCacheDB
    disk_error::Union{Nothing,Exception}
    status::Union{Nothing,ProjectCacheStatus}
    operation::Symbol
end

"""
A single snapshot of everything a watcher needs to show about a workspace's background work.

This is the stable contract between the engine and any watcher (the GUI today, scripts and workflows
later): watchers read `WorkspaceStatus` and nothing else about jobs, progress, or the cache. It is
recomputed only when [`refresh_status!`](@ref) observes active or just-stopped work, so an idle
render loop reads a cached value instead of rebuilding strings every frame.

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
    # The cache is the persistence side of the scan, not a separately cancellable job, so it needs
    # only a visible state and last error rather than a full WorkspaceJob.
    cache_state::Symbol
    cache_error::String
    # Last recoverable failure reported by the source watcher; cleared by the next good update.
    source_error::String
    work::WorkDependencyGraph
    background_processing::Bool
    background_tasks::Vector{Task}
    metrics::BuildMetrics
    profiler::Profiling.ProfileSession
    profile_restart_pending::Bool
    publish_lock::ReentrantLock
    idle_condition::Base.Threads.Condition
    status::WorkspaceStatus
    closed::Bool
end

"""
Create the empty state for one project-owned source.
"""
function Workspace(
    project::Project,
    source::S;
    profile_internal::Bool=Profiling.environment_flag("MB_PROFILE_INTERNAL"),
    profile_cpu::Bool=Profiling.environment_flag("MB_PROFILE_CPU"),
    profile_output::Union{Nothing,AbstractString}=
        Profiling.environment_path("MB_PROFILE_OUTPUT"),
    crash_trace::Union{Nothing,AbstractString}=
        Profiling.environment_path("MB_CRASH_TRACE"),
    rebuild::Bool=false,
    cache::Bool=true,
    background_processing::Bool=false,
)::Workspace{S} where {S<:AbstractDataSource}
    hierarchy = Hierarchy(source_id(source), false)
    identity = project_cache_identity(project_name(project), source)
    profiler = Profiling.ProfileSession(
        profile_internal, profile_cpu, profile_output, crash_trace)
    metrics = BuildMetrics()
    disk_error::Union{Nothing,Exception} = nothing
    cache_db::AbstractCacheDB = try
        cache ? open_cache_db(identity, profiler, metrics; rebuild) :
            open_memory_cache_db(identity, profiler, metrics)
    catch error
        if cache && !rebuild && error isa ProjectCacheSchemaError
            disk_error = error
            @warn(
                "Generated project cache is unavailable; continuing without disk cache",
                cache=identity.cache_path,
                error,
            )
            open_memory_cache_db(identity, profiler, metrics)
        else
            Profiling.close!(profiler)
            rethrow()
        end
    end
    publish_lock = ReentrantLock()
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
        WorkspaceCache(identity, cache_db, disk_error, nothing, :load),
        WorkspaceJob(),
        :idle,
        "",
        "",
        WorkDependencyGraph(),
        cache_db isa CacheDB && background_processing,
        Task[],
        metrics,
        profiler,
        false,
        publish_lock,
        Base.Threads.Condition(publish_lock),
        WorkspaceStatus(),
        false,
    )
    start_cache!(cache_db)
    start_work_workers!(workspace)
    return workspace
end

# ---------------------------------------------------------------------------
# Pretty printing
# ---------------------------------------------------------------------------

Base.show(io::IO, workspace::Workspace) = print(
    io,
    "Workspace(", project_name(workspace.project), ", ", source_label(workspace.source),
    ", ", length(workspace.index.items), " items)",
)

function Base.show(io::IO, ::MIME"text/plain", workspace::Workspace)
    println(io, "Workspace · ", project_name(workspace.project))
    println(io, "  source label: ", source_label(workspace.source))
    println(io, "  source:       ", source_id(workspace.source))
    println(io, "  scan:         ", workspace.scan.state)
    println(io, "  cache:        ", workspace.cache_state)
    println(io, "  items:        ", length(workspace.index.items))
    print(io,   "  failures:     ", length(workspace.index.analysis_errors))
    workspace.closed && print(io, "\n  (closed)")
end

include("Workspace/Operations.jl")
include("Workspace/Status.jl")
include("Workspace/DataAccess.jl")
include("Workspace/Processing.jl")
include("Workspace/MemoryDiagnostics.jl")

end
