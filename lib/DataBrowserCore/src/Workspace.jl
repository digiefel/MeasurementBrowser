module Workspace

using Printf
using DataBrowserAPI: @timed_dbg
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
using DataBrowserAPI:
    PipelineStage,
    COLLECTION_ANALYZE,
    COLLECTION_PROCESS,
    ITEM_ANALYZE,
    ITEM_PROCESS,
    SOURCE_READ,
    SOURCE_INTERPRET

using DataBrowserCache:
    CacheDB,
    BuildMetrics,
    CacheResultKey,
    CacheResultStatus,
    CacheStageSummary,
    ProjectCacheSchemaError,
    ProjectCacheDataError,
    ProjectCacheIdentity,
    ProjectCacheIndex,
    ProjectCacheStatus,
    RESULT_FAILED,
    RESULT_READY,
    cache_stage_summary,
    cache_has_pending_writes,
    cache_pending_counts,
    cached_result_state,
    clear_cache_index!,
    clear_cached_result_state!,
    close_cache_db!,
    delete_collection_records!,
    delete_collection_metadata!,
    delete_source_item!,
    delete_source_output!,
    has_payload,
    load_cache_index,
    cached_source_fingerprints,
    open_cache_db,
    project_cache_identity,
    read_payload,
    reset_build_metrics!,
    store_collection_metadata!,
    store_collection_index!,
    store_collection_process_result!,
    store_interpreted!,
    store_interpreted_data!,
    store_item_metadata!,
    store_item_metadata_layer!,
    store_processed!,
    store_result_failure!,
    store_source_identity!,
    store_source_read!,
    source_complete,
    cache_knows_source,
    source_item_key!
import DataBrowserCache
import DataBrowserCache: query_items, set_cache_memory_limit!
using DataBrowserAPI.ItemIndex:
    CollectionIndex,
    CollectionInput,
    ItemRecord,
    MetadataDict,
    append_item!,
    clear_collection_analysis!,
    collection_inputs,
    collection_item_ids,
    collection_path_keys,
    effective_metadata,
    effective_record,
    metadata_dict,
    collection_value_path,
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
    AbstractProject,
    SourceChanges,
    SourceError,
    analyze,
    annotate_collection_path,
    attach_record,
    close_source!,
    collection,
    fingerprint,
    id,
    item_data,
    label,
    metadata,
    process,
    read,
    reconstruct,
    open_source,
    project_name,
    source_id,
    source_items,
    source_item_noun,
    source_label,
    watch_source


"""
State of the workspace's source scan: running or idle, cancel token, and how many sources were found.

`epoch` increments each time a scan starts. A background scan task captures that value and ignores
its own publishes if a newer scan has begun.
"""
mutable struct WorkspaceJob
    epoch::Int
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
    analysis_errors::Dict{Union{String,Int64},String}
    # Published item ids per source-item key, so per-publish lookups avoid scanning every item.
    items_by_source::Dict{Int64,Vector{String}}
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
    db::CacheDB
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
mutable struct Workspace{P<:AbstractProject}
    project::P
    source::AbstractDataSource
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
    # Requested disk cache. Kept even when a schema error falls back to a memory cache, so a
    # rebuild can open a real disk cache again.
    disk_cache::Bool
    background_tasks::Vector{Task}
    metrics::BuildMetrics
    publish_lock::ReentrantLock
    # Held for the whole of `modify_workspace!` / `close_workspace!` so a frame cannot read a
    # half-rebuilt runtime. The GUI trylocks this at the start of a frame.
    lifecycle_lock::ReentrantLock
    idle_condition::Base.Threads.Condition
    status::WorkspaceStatus
    status_dirty::Base.Threads.Atomic{Bool}
    cancel_source::CancellationTokenSource
    closed::Bool
end

"""
Open the cache for one project/source pair, falling back to memory on a schema error.
"""
function _open_workspace_cache(
    project::AbstractProject,
    source::AbstractDataSource,
    metrics::BuildMetrics;
    rebuild::Bool,
    cache::Bool,
)::Tuple{CacheDB,Union{Nothing,Exception},ProjectCacheIdentity}
    identity = project_cache_identity(project_name(project), source)
    disk_error::Union{Nothing,Exception} = nothing
    cache_db::CacheDB = try
        @timed_dbg open_cache_db(identity, metrics; rebuild, persistent=cache)
    catch error
        if cache && !rebuild && error isa ProjectCacheSchemaError
            disk_error = error
            @warn(
                "Generated project cache is unavailable; continuing without disk cache",
                cache=identity.cache_path,
                error,
            )
            open_cache_db(identity, metrics; persistent=false)
        else
            rethrow()
        end
    end
    return cache_db, disk_error, identity
end

"""
Create the empty state for one project-owned source.
"""
function Workspace(
    project::P,
    source::AbstractDataSource;
    rebuild::Bool=false,
    cache::Bool=true,
    background_processing::Bool=false,
)::Workspace{P} where {P<:AbstractProject}
    collections = CollectionIndex(source_id(source))
    metrics = BuildMetrics()
    cache_db, disk_error, identity =
        _open_workspace_cache(project, source, metrics; rebuild, cache)
    publish_lock = ReentrantLock()
    workspace = Workspace(
        project,
        source,
        WorkspaceIndex(
            collections,
            Dict{String,ItemRecord}(),
            Dict{String,Dict{Symbol,Any}}(),
            Symbol[],
            Dict{Union{String,Int64},String}(),
            Dict{Int64,Vector{String}}(),
        ),
        WorkspaceSelection(),
        WorkspaceCache(identity, cache_db, disk_error, nothing, :load),
        WorkspaceJob(),
        :idle,
        "",
        "",
        WorkDependencyGraph(),
        background_processing,
        cache,
        Task[],
        metrics,
        publish_lock,
        ReentrantLock(),
        Base.Threads.Condition(publish_lock),
        WorkspaceStatus(),
        Base.Threads.Atomic{Bool}(true),
        CancellationTokenSource(),
        false,
    )
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
    project::AbstractProject,
    root_path::AbstractString;
    recursive::Bool=true,
    metadata_file::Union{Nothing,AbstractString}=DataBrowserSources.DEFAULT_DIRECTORY_METADATA_FILE,
    rebuild::Bool=false,
    cache::Bool=true,
    background_processing::Bool=false,
)::Workspace
    return open_workspace(
        project,
        DataBrowserSources.DirectorySource(root_path; recursive, metadata_file);
        rebuild,
        cache,
        background_processing,
    )
end

end
