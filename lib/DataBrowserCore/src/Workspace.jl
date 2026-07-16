module Workspace

using Printf
import DataBrowserProfiling as Profiling
using DataBrowserSources
using CancellationTokens:
    CancellationToken,
    CancellationTokenSource,
    OperationCanceledException,
    cancel,
    get_token,
    is_cancellation_requested

using ..WorkGraph:
    WorkDependencyGraph,
    WorkKey,
    WorkKind,
    WorkNode,
    bump_revision!,
    current_revision,
    dependencies_ready,
    pop_queued_node!,
    push_queue_entry!,
    queue_ready_node!,
    seed_node_dependencies!,
    take_work!,
    wake_ready_dependents!
import ..WorkGraph:
    COLLECTION_ANALYZE,
    COLLECTION_PROCESS,
    ITEM_ANALYZE,
    ITEM_PROCESS,
    SOURCE_INTERPRET

using DataBrowserCache:
    AbstractCacheDB,
    BuildMetrics,
    CacheResultKey,
    CacheResultKind,
    CacheResultStatus,
    CacheDB,
    CacheStageSummary,
    COLLECTION_ANALYSIS_RESULT,
    COLLECTION_PROCESS_RESULT,
    ITEM_ANALYSIS_RESULT,
    PROCESSING_RESULT,
    ProjectCacheSchemaError,
    ProjectCacheDataError,
    ProjectCacheIdentity,
    ProjectCacheIndex,
    ProjectCacheStatus,
    RESULT_FAILED,
    RESULT_READY,
    cache_stage_summary,
    cache_built,
    cache_has_pending_writes,
    cache_pending_counts,
    clear_cache_index!,
    clear_cached_result_state!,
    clear_cached_source_state!,
    close_cache_db!,
    delete_collection_records!,
    delete_collection_metadata!,
    delete_source_item!,
    load_cache_index,
    _load_source_item_fingerprints,
    open_memory_cache_db,
    open_cache_db,
    project_cache_identity,
    query_items,
    read_item_data,
    record_cache_phase!,
    reset_build_metrics!,
    set_cache_memory_limit!,
    start_cache!,
    stop_cache!,
    store_collection_metadata!,
    store_collection_index!,
    store_collection_process_result!,
    store_interpreted!,
    store_item_metadata!,
    store_processed!,
    store_result_failure!,
    store_source_item_failure!,
    edit_source_item_metadata!,
    wait_condition_deadline,
    write_meta_header!
import DataBrowserCache: query_items, read_item_data, set_cache_memory_limit!
using DataBrowserAPI.ItemIndex:
    CollectionIndex,
    CollectionInput,
    ItemFailure,
    ItemRecord,
    MetadataDict,
    RegisteredDataItem,
    SourceScan,
    append_item!,
    clear_collection_analysis!,
    collection_inputs,
    collection_item_ids,
    collection_path_keys,
    effective_metadata,
    effective_record,
    metadata_dict,
    registered_collection_path,
    registration_names,
    remove_item!,
    resolve_collection_path!,
    resolve_collection_paths!,
    set_collection_analysis!
using ..DataBrowserCore:
    SourceItemInterpretation,
    interpret_source_item
import DataBrowserAPI
import DataBrowserAPI:
    AbstractCollection,
    AbstractDataSource,
    AbstractDataSourceItem,
    AbstractDataItem,
    Project,
    SourceChanges,
    SourceError,
    _analyze_collection,
    _analyze_item,
    close_source!,
    _has_collection_analysis,
    _has_collection_process,
    _process_collection,
    cacheable,
    fingerprint,
    id,
    item_data,
    label,
    metadata,
    process,
    open_source,
    project_name,
    record_scan_phase!,
    reset_scan_profile!,
    scan_profile_summary,
    source_id,
    source_items,
    source_item_noun,
    source_label,
    source_open_options,
    watch_source


"""
One cancellable workspace operation and its latest visible state.
"""
mutable struct WorkspaceJob
    id::Int
    state::Symbol
    error::String
    cancel_token::Union{Nothing,CancellationTokenSource}
    discovered::Base.Threads.Atomic{Int}
end

WorkspaceJob()::WorkspaceJob = WorkspaceJob(0, :idle, "", nothing, Base.Threads.Atomic{Int}(0))

"""
The progressively populated item index for one open source.
"""
mutable struct WorkspaceIndex
    collections::CollectionIndex
    items::Dict{String,ItemRecord}
    # The computed metadata layers per item (analyze output merged with any collection-process
    # overwrite); the entries layer stays on the record.
    item_metadata::Dict{String,Dict{Symbol,Any}}
    collection_metadata_keys::Vector{Symbol}
    source::Union{Nothing,SourceScan}
    analysis_errors::Dict{Union{String,Int64},String}
    # Published item ids per source item, so per-publish lookups avoid scanning every item.
    items_by_source::Dict{String,Vector{String}}
end

"""
Stable selection identities owned by a workspace.
"""
Base.@kwdef mutable struct WorkspaceSelection
    collection_ids::Vector{String} = String[]
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

"""Counts row shown by status watchers."""
struct WorkspaceStageCounts
    source_noun::String
    sources_found::Int
    sources_pending::Int
    cache::CacheStageSummary
end

WorkspaceStageCounts()::WorkspaceStageCounts =
    WorkspaceStageCounts("source items", 0, 0, CacheStageSummary())

"""
A single snapshot of everything a watcher needs to show about a workspace's background work.

This is the stable contract between the engine and any watcher (the GUI today, scripts and workflows
later): watchers read `WorkspaceStatus` and nothing else about jobs, progress, or the cache. It is
recomputed only after engine publications and discovery progress, so an idle render loop reads a
cached value instead of rebuilding strings every frame.

- `level` drives the watcher's color/emphasis: `:none`, `:busy`, `:fresh`, `:stale`, `:missing`,
  `:error`.
- `label` is a short word for a button or chip ("Building", "Fresh", "Errors").
- `detail` is the one merged human line: the live activity while `busy`, otherwise a short state.
- `busy` is true while any scan, analysis, or cache work runs.
- `progress` is a determinate fraction when counts are known, or `nothing` for an indeterminate or
  absent bar.
- `counts` carries the cache/source numbers shown by status watchers.
- `errors` lists source-item failures as `id => first message line`, streamed as they occur.
"""
struct WorkspaceStatus
    level::Symbol
    label::String
    detail::String
    busy::Bool
    progress::Union{Nothing,Float32}
    counts::WorkspaceStageCounts
    errors::Vector{Pair{String,String}}
end

WorkspaceStatus() =
    WorkspaceStatus(
        :none, "Opening", "Opening the source…", true, nothing,
        WorkspaceStageCounts(), Pair{String,String}[])


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
    # The effective construction options `open_workspace` was called with, replayed verbatim
    # (splatted) when the browser reopens an equivalent workspace on the same source root.
    open_options::NamedTuple
    background_tasks::Vector{Task}
    metrics::BuildMetrics
    profiler::Profiling.ProfileSession
    profile_restart_pending::Bool
    publish_lock::ReentrantLock
    idle_condition::Base.Threads.Condition
    status::WorkspaceStatus
    status_dirty::Base.Threads.Atomic{Bool}
    cancel_source::CancellationTokenSource
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
    collections = CollectionIndex(source_id(source))
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
    open_options = (;
        profile_internal,
        profile_cpu,
        profile_output,
        crash_trace,
        cache,
        background_processing,
        source_open_options(source)...,
    )
    workspace = Workspace(
        project,
        source,
        WorkspaceIndex(
            collections,
            Dict{String,ItemRecord}(),
            Dict{String,Dict{Symbol,Any}}(),
            Symbol[],
            nothing,
            Dict{Union{String,Int64},String}(),
            Dict{String,Vector{String}}(),
        ),
        WorkspaceSelection(),
        WorkspaceCache(identity, cache_db, disk_error, nothing, :load),
        WorkspaceJob(),
        :idle,
        "",
        "",
        WorkDependencyGraph(),
        background_processing,
        open_options,
        Task[],
        metrics,
        profiler,
        false,
        publish_lock,
        Base.Threads.Condition(publish_lock),
        WorkspaceStatus(),
        Base.Threads.Atomic{Bool}(true),
        CancellationTokenSource(),
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

function open_workspace(
    project::Project,
    root_path::AbstractString;
    recursive::Bool=true,
    metadata_file::Union{Nothing,AbstractString}=DataBrowserSources.DEFAULT_DIRECTORY_METADATA_FILE,
    profile_internal::Bool=Profiling.environment_flag("MB_PROFILE_INTERNAL"),
    profile_cpu::Bool=Profiling.environment_flag("MB_PROFILE_CPU"),
    profile_output::Union{Nothing,AbstractString}=
        Profiling.environment_path("MB_PROFILE_OUTPUT"),
    crash_trace::Union{Nothing,AbstractString}=
        Profiling.environment_path("MB_CRASH_TRACE"),
    rebuild::Bool=false,
    cache::Bool=true,
    background_processing::Bool=false,
)::Workspace
    return open_workspace(
        project,
        DataBrowserSources.DirectorySource(root_path; recursive, metadata_file);
        profile_internal,
        profile_cpu,
        profile_output,
        crash_trace,
        rebuild,
        cache,
        background_processing,
    )
end

end
