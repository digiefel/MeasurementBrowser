"""
Start cache loading and source scanning for one new workspace.
"""
function open_workspace(
    project::Project,
    source::AbstractDataSource;
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
    opened_source = open_source(source)
    workspace = try
        Workspace(
            project,
            opened_source;
            profile_internal,
            profile_cpu,
            profile_output,
            crash_trace,
            rebuild,
            cache,
            background_processing,
        )
    catch
        close_source!(opened_source)
        rethrow()
    end
    try
        scan_source!(workspace; rebuild)
        watch_source(
            opened_source,
            event -> publish_source_event!(workspace, event),
        )
        return workspace
    catch
        close_workspace!(workspace)
        rethrow()
    end
end

"""
Cancel all work owned by a workspace and wait for it to stop.
"""
function close_workspace!(workspace::Workspace)::Nothing
    workspace.closed && return nothing
    workspace.closed = true
    cancel_scan!(workspace)
    for task in workspace.background_tasks
        istaskdone(task) || wait(task)
    end
    stop_work_workers!(workspace)
    try
        Profiling.close!(workspace.profiler)
    finally
        try
            close_cache_db!(workspace.cache.db)
        finally
            close_source!(workspace.source)
        end
    end
    return nothing
end

"""Replace the selection with the supplied indexed item records."""
function select_items!(
    workspace::Workspace,
    items::AbstractVector{ItemRecord},
)::Nothing
    workspace.selection.item_ids = [item.id for item in items]
    return nothing
end

"""Replace the selection with exact item ids in this workspace."""
function select_items!(
    workspace::Workspace,
    ids::AbstractVector{<:AbstractString},
)::Nothing
    selected = String[id for id in ids]
    missing = String[id for id in selected if !haskey(workspace.index.items, id)]
    isempty(missing) || error(
        "Cannot select $(length(missing)) item id(s) that are not in this workspace: " *
        join(missing, ", "),
    )
    workspace.selection.item_ids = selected
    return nothing
end

"""Replace the selection with data items whose ids are present in the workspace."""
function select_items!(
    workspace::Workspace,
    items::AbstractVector{<:AbstractDataItem},
)::Nothing
    selected = String[]
    for item in items
        selected_id = id(item)
        haskey(workspace.index.items, selected_id) || error(
            "Cannot select item id '$selected_id': no indexed item with that id exists in this workspace",
        )
        push!(selected, selected_id)
    end

    workspace.selection.item_ids = selected
    return nothing
end

"""Clear the selection, or fail clearly for unsupported item selectors."""
function select_items!(
    workspace::Workspace,
    items::AbstractVector,
)::Nothing
    if isempty(items)
        empty!(workspace.selection.item_ids)
        return nothing
    end
    if all(item -> item isa ItemRecord, items)
        return select_items!(workspace, ItemRecord[item for item in items])
    end
    if all(item -> item isa AbstractString, items)
        return select_items!(workspace, String[item for item in items])
    end
    if all(item -> item isa AbstractDataItem, items)
        return select_items!(workspace, AbstractDataItem[item for item in items])
    end
    error(
        "Cannot select values of type $(eltype(items)); pass ItemRecord values, " *
        "exact item ids, or AbstractDataItem values",
    )
end

source_scan_running(workspace::Workspace)::Bool =
    workspace.scan.state in (:discovering, :canceling)

cache_work_running(workspace::Workspace)::Bool =
    workspace.cache_state in (:loading, :writing, :canceling)

analysis_work_running(workspace::Workspace)::Bool =
    work_kind_running(workspace, COLLECTION_PROCESS) ||
    work_kind_running(workspace, COLLECTION_ANALYZE)

"""Whether any item is waiting for or running processing, or its result is still being written."""
function processing_work_running(workspace::Workspace)::Bool
    return work_kind_running(workspace, ITEM_PROCESS) ||
        work_kind_running(workspace, ITEM_ANALYZE) ||
        cache_has_pending_writes(workspace.cache.db)
end

"""Whether the scan or any queued or running work-graph node is still in flight."""
function engine_work_running(workspace::Workspace)::Bool
    _, _, active = work_counts(workspace)
    return active > 0 || source_scan_running(workspace) || cache_work_running(workspace)
end

"""
Whether a watcher should keep refreshing: engine work is in flight or writes are still flushing.

Pending cache writes are display-only busyness — reads overlay the write buffer and appends
backpressure on row pressure, so nothing in the engine gates on flush completion.
"""
workspace_busy(workspace::Workspace)::Bool =
    engine_work_running(workspace) || cache_has_pending_writes(workspace.cache.db)

"""Keep a background task alive until workspace shutdown."""
function track_task!(workspace::Workspace, task::Task)::Task
    filter!(!istaskdone, workspace.background_tasks)
    push!(workspace.background_tasks, task)
    return task
end

"""
Request cancellation of a running source scan.

Only a scan that is actually running moves into `:canceling`; the check and the transition share
the publish lock with the scan task's own settlement, so a settled job can never re-enter
`:canceling`.
"""
function cancel_scan!(workspace::Workspace)::Nothing
    lock(workspace.publish_lock) do
        job = workspace.scan
        job.state === :discovering || return
        Base.Threads.atomic_xchg!(job.cancel_token, true)
        job.state = :canceling
    end
    return nothing
end
cancel_analysis!(workspace::Workspace)::Nothing = (cancel_waiting_work!(workspace); nothing)

"""
Cancel in-flight cache work.

Cache writes are now produced by the source scan itself, so cancelling cache work cancels the scan.
"""
cancel_cache!(workspace::Workspace)::Nothing = cancel_scan!(workspace)

"""Return records currently published for one source item."""
function source_item_records(index::WorkspaceIndex, source_item_id_value::String)::Vector{ItemRecord}
    ids = get(index.items_by_source, source_item_id_value, nothing)
    ids === nothing && return ItemRecord[]
    return ItemRecord[index.items[id] for id in ids if haskey(index.items, id)]
end

"""Return collection keys touched by records and their ancestors."""
function affected_collection_keys(records::Vector{ItemRecord})::Vector{String}
    keys = String[]
    for record in records
        for depth in eachindex(record.collection)
            push!(keys, collection_path_key(record.collection[1:depth]))
        end
    end
    return unique(keys)
end

"""Refresh the collection-metadata key list from the current hierarchy."""
function refresh_collection_metadata_keys!(index::WorkspaceIndex)::Nothing
    metadata_keys = Set{Symbol}()
    for node in values(index.hierarchy.index)
        union!(metadata_keys, keys(node.metadata))
    end
    index.collection_metadata_keys = sort!(collect(metadata_keys); by=String)
    return nothing
end

"""Rebuild hierarchy indexes from the currently published records."""
function rebuild_workspace_hierarchy!(workspace::Workspace)::Nothing
    records = sort!(collect(values(workspace.index.items)); by=record -> (record.source_item_id, record.id))
    workspace.index.hierarchy = Hierarchy(records, workspace.source)
    refresh_collection_metadata_keys!(workspace.index)
    return nothing
end

"""
Reapply source collection metadata, bumping and re-marking only items whose merged values changed.

The per-record diff is the event filter: an item whose effective metadata is unchanged keeps its
results. Changed items bump their work keys' revisions (discarding stale completions) and re-enqueue;
recompute overwrites the cache rows, so nothing invalidates the cache here.
"""
function apply_collection_metadata_changes!(workspace::Workspace)::Bool
    old_hierarchy = workspace.index.hierarchy
    records = sort!(
        collect(values(workspace.index.items));
        by=record -> (record.source_item_id, record.id),
    )
    old_metadata = Dict(
        record.id => effective_metadata(old_hierarchy, record) for record in records
    )
    new_hierarchy = Hierarchy(records, workspace.source, old_hierarchy.skipped_count)
    changed = ItemRecord[
        record for record in records
        if old_metadata[record.id] != effective_metadata(new_hierarchy, record)
    ]
    old_node_analysis = Dict(path => copy(node.analysis) for (path, node) in old_hierarchy.index)
    for (path, analysis) in old_node_analysis
        node = get(new_hierarchy.index, path, nothing)
        node === nothing || merge!(node.analysis, analysis)
    end
    workspace.index.hierarchy = new_hierarchy
    refresh_collection_metadata_keys!(workspace.index)

    isempty(changed) && return false
    invalidate_records_work!(workspace, changed)
    return true
end

"""
Bump and re-mark the work of a set of records and their collections after their inputs changed.

Shared by the in-session metadata event and the reopen source-metadata diff: each record's process
and analyze results are bumped off their prior revisions (a still-:ready ITEM_PROCESS would answer
from a stale payload), and the ancestor collections' process/analyze results are marked stale.
"""
function invalidate_records_work!(workspace::Workspace, changed::Vector{ItemRecord})::Nothing
    affected = affected_collection_keys(changed)
    lock(workspace.work.lock) do
        for record in changed
            delete!(workspace.index.item_metadata, record.id)
            delete!(workspace.index.analysis_errors, record.id)
            bump_work_missing!(workspace, WorkKey(ITEM_PROCESS, record.id))
            bump_work_missing!(
                workspace, WorkKey(ITEM_ANALYZE, record.id);
                dependencies=WorkKey[WorkKey(ITEM_PROCESS, record.id)])
            workspace.background_processing && enqueue_processing!(workspace, record)
        end
    end
    for key in affected
        node = get(workspace.index.hierarchy.index, collection_path_tuple(key), nothing)
        node === nothing || empty!(node.analysis)
    end
    mark_collection_work_missing!(workspace, affected)
    return nothing
end

"""Refresh the authoritative published source snapshot from current index state."""
function refresh_workspace_source!(
    workspace::Workspace,
    fingerprints::AbstractDict{String},
)::Nothing
    metadata_fingerprints = Dict{String,Any}(collection_metadata_fingerprints(workspace.source))
    workspace.index.source = SourceScan(
        source_id(workspace.source),
        source_label(workspace.source),
        Dict{String,Any}(fingerprints),
        metadata_fingerprints,
        workspace.index.hierarchy,
        ItemFailure[],
    )
    store_collection_metadata_fingerprints!(workspace.cache.db, metadata_fingerprints)
    return nothing
end

"""Remove all published output owned by one source item and invalidate affected collections."""
function remove_source_item_output!(
    workspace::Workspace,
    source_item_id_value::String,
)::Tuple{Vector{ItemRecord},Vector{String}}
    old_records = source_item_records(workspace.index, source_item_id_value)
    invalidated = affected_collection_keys(old_records)
    for record in old_records
        delete!(workspace.index.items, record.id)
        delete!(workspace.index.item_metadata, record.id)
        delete!(workspace.index.analysis_errors, record.id)
    end
    delete!(workspace.index.analysis_errors, source_item_id_value)
    delete!(workspace.index.items_by_source, source_item_id_value)
    filter!(id -> haskey(workspace.index.items, id), workspace.selection.item_ids)
    isempty(old_records) && return old_records, invalidated
    edit = edit_hierarchy(workspace.index.hierarchy)
    remove_records!(edit, old_records)
    for key in invalidated
        clear_node_analysis!(edit, collection_path_tuple(key))
    end
    workspace.index.hierarchy = finish_edit!(edit, workspace.source)
    edit_changed_structure(edit) && refresh_collection_metadata_keys!(workspace.index)
    return old_records, invalidated
end

"""Publish one source item's current records after validating duplicate IDs."""
function publish_source_item_records!(
    workspace::Workspace,
    source_item_id_value::String,
    records::Vector{ItemRecord},
)::Tuple{Vector{ItemRecord},Vector{String}}
    old_records = source_item_records(workspace.index, source_item_id_value)
    old_ids = Set(record.id for record in old_records)
    replacement_ids = Set{String}()
    for record in records
        record.source_item_id == source_item_id_value || error(
            "Cannot publish source item '$source_item_id_value': replacement contains " *
            "record '$(record.id)' from source item '$(record.source_item_id)'",
        )
        record.id in replacement_ids && error(
            "Cannot publish source item '$source_item_id_value': duplicate item id '$(record.id)'",
        )
        push!(replacement_ids, record.id)
        haskey(workspace.index.items, record.id) && !(record.id in old_ids) && error(
            "Cannot publish source item '$source_item_id_value': replacement item id " *
            "'$(record.id)' already belongs to another source item",
        )
    end

    old_keys = affected_collection_keys(old_records)
    new_keys = affected_collection_keys(records)
    for record in old_records
        delete!(workspace.index.items, record.id)
        delete!(workspace.index.item_metadata, record.id)
        delete!(workspace.index.analysis_errors, record.id)
    end
    delete!(workspace.index.analysis_errors, source_item_id_value)
    for record in records
        workspace.index.items[record.id] = record
        delete!(workspace.index.item_metadata, record.id)
    end
    if isempty(records)
        delete!(workspace.index.items_by_source, source_item_id_value)
    else
        workspace.index.items_by_source[source_item_id_value] =
            String[record.id for record in records]
    end
    filter!(id -> haskey(workspace.index.items, id), workspace.selection.item_ids)
    invalidated = unique([old_keys; new_keys])
    edit = edit_hierarchy(workspace.index.hierarchy)
    remove_records!(edit, old_records)
    foreach(record -> insert_record!(edit, record), records)
    for key in invalidated
        clear_node_analysis!(edit, collection_path_tuple(key))
    end
    workspace.index.hierarchy = finish_edit!(edit, workspace.source)
    edit_changed_structure(edit) && refresh_collection_metadata_keys!(workspace.index)
    return old_records, invalidated
end

"""
Submit one source-change batch to the work graph and remove deleted published output. Cache
deletes are buffered mutations, so ordering against re-interpretation is kept by the buffers:
a delete enqueued here always precedes the rows a later interpretation appends.
"""
function apply_source_changes!(
    workspace::Workspace,
    changes::SourceChanges,
    fingerprints::AbstractDict{String},
    status::Union{Nothing,ProjectCacheStatus},
)::Bool
    changed = false
    for removed in changes.removals
        old_records, invalidated = remove_source_item_output!(workspace, removed)
        delete_source_item!(workspace.cache.db, removed, old_records)
        delete_collection_metadata!(workspace.cache.db, invalidated)
        mark_collection_work_missing!(workspace, invalidated)
        lock(workspace.work.lock) do
            delete!(workspace.work.source_items, removed)
        end
        changed = true
    end
    isempty(changes.upserts) || lock(workspace.work.lock) do
        for source_item in changes.upserts
            source_item_id_value = source_item_id(source_item)
            old_records = source_item_records(workspace.index, source_item_id_value)
            delete_source_item!(workspace.cache.db, source_item_id_value, old_records)
            key = WorkKey(SOURCE_INTERPRET, source_item_id_value)
            workspace.work.source_items[source_item_id_value] = source_item
            enqueue_work!(
                workspace,
                key,
                bump_revision!(workspace.work, key);
                priority=3,
            )
            changed = true
        end
    end
    changes.metadata_changed && (changed = apply_collection_metadata_changes!(workspace) || changed)
    refresh_workspace_source!(workspace, fingerprints)
    status === nothing || (workspace.cache.status = status)
    return changed
end

"""
Start one source scan that progressively populates the workspace index and cache.

The scan first surfaces any already-cached index for an instant first view, then streams freshly
interpreted source items into the workspace and a bounded cache-write batch. The batch flushes by
source-item count, cached row count, or age, so progressive saves remain frequent without paying for
one DuckDB transaction per source item. `rebuild=true` first wipes the cache and ignores the prior
scan, forcing a full re-interpretation.
"""
function scan_source!(
    workspace::Workspace;
    rebuild::Bool=false,
)::Nothing
    cachedb = workspace.cache.db
    job = workspace.scan
    # Job setup runs under the publish lock so cancel_scan!'s running check can never interleave
    # with the transition into :discovering.
    scan_id, cancel_token = lock(workspace.publish_lock) do
        source_scan_running(workspace) && error("A source scan is already running")
        workspace.closed && error("The workspace is closed")
        job.id += 1
        job.state = :discovering
        job.error = ""
        workspace.cache_state = :loading
        workspace.cache_error = ""
        workspace.cache.operation = rebuild ? :rebuild : :update
        reset_build_metrics!(workspace.metrics)
        reset_scan_profile!(workspace.project)
        # Detach the previous scan so this one streams into a fresh hierarchy: progressive results
        # update the displayed tree live while any still-finishing analysis from the previous scan
        # keeps reading its own, now-detached hierarchy — no shared mutable state, no race. Errors
        # reset so stale ones from the previous scan never linger.
        workspace.index.source = nothing
        empty!(workspace.index.analysis_errors)
        job.cancel_token = Base.Threads.Atomic{Bool}(false)
        (job.id, job.cancel_token)
    end
    task = Base.Threads.@spawn begin
        Profiling.@profile_span workspace.profiler :source :scan Profiling.ProfileAttributes(
            source_id=source_id(workspace.source),
        ) begin
        try
            rebuild && Profiling.@profile_span workspace.profiler :source :cache_clear Profiling.ProfileAttributes() begin
                clear_cache_index!(cachedb)
            end
            cached = !rebuild && cache_built(cachedb) ?
                Profiling.@profile_span(workspace.profiler, :source, :cache_load,
                    Profiling.ProfileAttributes(), load_cache_index(cachedb)) :
                nothing
            publish_cache_state!(
                workspace, scan_id, cached === nothing ? :missing : :ready, cached)

            write_meta_header!(cachedb)
            discovered = Profiling.@profile_span workspace.profiler :source :discover Profiling.ProfileAttributes() begin
                with_cancel(() -> cancel_token[]) do
                    source_items(workspace.source)
                end
            end
            current = Dict(source_item_id(item) => fingerprint(item) for item in discovered)
            previous = cached === nothing ?
                Dict{String,Any}() :
                copy(cached.source.source_item_fingerprints)
            upserts = AbstractDataSourceItem[]
            for item in discovered
                id_value = source_item_id(item)
                # A `nothing` fingerprint means the source cannot prove the item is unchanged,
                # so it is always re-read.
                if !haskey(previous, id_value) || current[id_value] === nothing ||
                   !isequal(previous[id_value], current[id_value])
                    push!(upserts, item)
                end
            end
            present = Set(keys(current))
            removals = String[id for id in keys(previous) if !(id in present)]
            # Reopen replays missed collection-metadata edits from the source-exposed fingerprints:
            # a changed subtree marks its items' and collections' results stale, exactly as an
            # in-session metadata change would.
            previous_metadata = cached === nothing ?
                Dict{String,Any}() :
                copy(cached.source.collection_metadata_fingerprints)
            current_metadata = Dict{String,Any}(
                collection_metadata_fingerprints(workspace.source))
            changed_collections = String[
                key for key in union(keys(previous_metadata), keys(current_metadata))
                if !isequal(
                    get(previous_metadata, key, nothing), get(current_metadata, key, nothing))
            ]
            isempty(changed_collections) ||
                invalidate_changed_metadata!(workspace, scan_id, changed_collections)
            stale = count(id -> haskey(previous, id) && !isequal(previous[id], current[id]), keys(current))
            new_items = count(id -> !haskey(previous, id), keys(current))
            status = ProjectCacheStatus(
                length(current),
                length(previous),
                length(current) - length(workspace.index.analysis_errors),
                stale,
                new_items,
                length(removals),
                length(workspace.index.analysis_errors),
            )
            Profiling.@profile_span workspace.profiler :source :publish_changes Profiling.ProfileAttributes(
                items=length(upserts),
                batch_size=length(removals),
            ) begin
                publish_source_changes!(
                    workspace,
                    scan_id,
                    SourceChanges(upserts, removals; metadata_changed=false),
                    current,
                    status,
                )
            end
            publish_scan_end!(
                workspace, scan_id;
                cache_hit=isempty(upserts) && isempty(removals),
            )
        catch error
            if is_job_cancelled(error)
                publish_scan_end!(workspace, scan_id; canceled=true)
            else
                publish_scan_end!(
                    workspace, scan_id; error, backtrace=catch_backtrace())
            end
        end
        end
    end
    track_task!(workspace, task)
    return nothing
end

"""Wipe the cache and re-scan the source from scratch."""
rebuild_cache!(workspace::Workspace)::Nothing = scan_source!(workspace; rebuild=true)

"""
Return the ids of items whose delivered metadata satisfies a SQL predicate over the query view.

Columns are exactly the metadata names project callbacks emit. Queries committed cache state, so a
value written in the last ≤2s may not yet appear. Memory-only workspaces return an empty vector.
"""
query_items(workspace::Workspace, predicate::AbstractString)::Vector{String} =
    query_items(workspace.cache.db, predicate)

"""Resize this workspace's DuckDB cache buffer-pool limit (MiB) live; suits a GUI control."""
set_cache_memory_limit!(workspace::Workspace, mib::Integer)::Int =
    set_cache_memory_limit!(workspace.cache.db, mib)

"""Return one lock-consistent counter snapshot for the structured profiler."""
function profile_counter_snapshot(workspace::Workspace)::NamedTuple
    completed, total, jobs = work_counts(workspace)
    depth = jobs + cache_pending_counts(workspace.cache.db).items
    return (
        # Scan-level progress counters are not populated; the schema keeps their columns.
        scan_done=0,
        scan_total=0,
        processing_done=completed,
        processing_total=total,
        queue_depth=depth,
    )
end

"""Begin recording immediately when idle, or prepare one clean rebuild after active work stops."""
function start_internal_profile!(workspace::Workspace)::Nothing
    profiler = workspace.profiler
    Profiling.validate_start(profiler)
    # Deciding busy, arming the restart, and canceling must be one atomic step against
    # finish_publish!: a completion landing in between would either strand the profiler in
    # :preparing forever or let the cancels kill the freshly restarted profiled rebuild.
    lock(workspace.publish_lock) do
        busy = source_scan_running(workspace) || processing_work_running(workspace) ||
            analysis_work_running(workspace) || cache_work_running(workspace)
        if busy
            profiler.state = :preparing
            workspace.profile_restart_pending = true
            cancel_scan!(workspace)
            cancel_analysis!(workspace)
            cancel_waiting_work!(workspace)
        else
            Profiling.start!(profiler, () -> profile_counter_snapshot(workspace))
        end
    end
    return nothing
end

"""Stop an active internal profile without canceling application work."""
function stop_internal_profile!(
    workspace::Workspace,
)::Union{Nothing,Profiling.ProfileReport}
    profiler = workspace.profiler
    # Aborting a pending restart must be atomic against finish_publish!: a completion landing
    # between the state check and the clear would start a profiled rebuild this stop then clobbers.
    aborted_preparation = lock(workspace.publish_lock) do
        profiler.state === :preparing || return false
        workspace.profile_restart_pending = false
        profiler.state = :idle
        return true
    end
    aborted_preparation && return nothing
    return Profiling.stop!(profiler)
end

"""Clear one stopped internal profile report."""
reset_internal_profile!(workspace::Workspace)::Nothing =
    Profiling.reset!(workspace.profiler)

"""Export one stopped internal profile as Chrome/Perfetto JSON."""
function export_internal_profile!(
    workspace::Workspace,
    path::AbstractString,
)::String
    return Profiling.export!(workspace.profiler, path)
end

"""
Status of a cache that now mirrors `source` exactly.

Every source item is cached; `fresh` counts those without an analysis error. Relative to the cache
that was loaded before this scan, `new`/`stale`/`deleted` count the source items that were added,
re-read because their fingerprint changed, and removed — the work an incremental reopen actually did.
"""
function _post_write_status(
    identity::ProjectCacheIdentity,
    source::SourceScan,
    cached::Union{Nothing,ProjectCacheIndex},
)::ProjectCacheStatus
    current = source.source_item_fingerprints
    total = length(current)
    errors = length(ProjectCacheIndex(identity, source).analysis_errors)
    fresh = total - errors
    if cached === nothing
        return ProjectCacheStatus(total, total, fresh, 0, total, 0, errors)
    end
    previous = cached.source.source_item_fingerprints
    new_items = count(id -> !haskey(previous, id), keys(current))
    stale = count(
        id -> haskey(previous, id) && !isequal(previous[id], current[id]), keys(current))
    deleted = count(id -> !haskey(current, id), keys(previous))
    return ProjectCacheStatus(total, total, fresh, stale, new_items, deleted, errors)
end

"""Replace the workspace item index with one complete hierarchy."""
function replace_item_index!(
    workspace::Workspace,
    hierarchy::Hierarchy,
)::Nothing
    items = Dict{String,ItemRecord}()
    sizehint!(items, length(hierarchy.all_items))
    items_by_source = Dict{String,Vector{String}}()
    metadata_keys = Set{Symbol}()
    for item in hierarchy.all_items
        haskey(items, item.id) && error("Duplicate item id generated during scan: $(item.id)")
        items[item.id] = item
        push!(get!(() -> String[], items_by_source, item.source_item_id), item.id)
    end
    for node in values(hierarchy.index)
        union!(metadata_keys, keys(node.metadata))
    end
    item_metadata = Dict{String,Dict{Symbol,Any}}()
    for item in hierarchy.all_items
        previous = get(workspace.index.item_metadata, item.id, nothing)
        previous === nothing || (item_metadata[item.id] = previous)
    end
    workspace.index = WorkspaceIndex(
        hierarchy,
        items,
        item_metadata,
        sort!(collect(metadata_keys); by=String),
        workspace.index.source,
        workspace.index.analysis_errors,
        items_by_source,
    )
    return nothing
end

"""Surface the cached index as an instant, independent first view, before the fresh scan refines it."""
function apply_cache_index!(
    workspace::Workspace,
    index::ProjectCacheIndex,
)::Nothing
    index.identity.source_id == source_id(workspace.source) ||
        error("Loaded cache belongs to $(index.identity.source_id), not $(source_id(workspace.source))")
    workspace.cache.identity = index.identity

    # The cache owns item-local records only. Rebuild the hierarchy with the open source's current
    # metadata so cached and memory-only startup produce the same effective inputs.
    hierarchy = Hierarchy(
        index.source.hierarchy.all_items,
        workspace.source,
        index.source.hierarchy.skipped_count,
    )
    replace_item_index!(workspace, hierarchy)
    # The cache restores the computed metadata layer for every item; the source-fingerprint diff in
    # scan_source! then invalidates the ones downstream of changed sources or changed metadata.
    workspace.index.item_metadata = Dict{String,Dict{Symbol,Any}}(
        id => copy(dict) for (id, dict) in index.item_metadata)
    for (path, cached_node) in index.source.hierarchy.index
        node = get(hierarchy.index, path, nothing)
        node === nothing && continue
        merge!(node.analysis, cached_node.analysis)
    end
    workspace.index.source = SourceScan(
        index.source.source_id,
        index.source.source_label,
        index.source.source_item_fingerprints,
        index.source.collection_metadata_fingerprints,
        hierarchy,
        index.source.analysis_failures,
    )

    cached_items = length(index.source.source_item_fingerprints)
    workspace.cache.status = ProjectCacheStatus(
        cached_items, cached_items, 0, 0, 0, 0, length(index.analysis_errors))
    workspace.index.analysis_errors = copy(index.analysis_errors)
    seed_cached_work_nodes!(workspace, index)
    for state in values(index.result_states)
        CacheResultStatus(state.status) === RESULT_FAILED || continue
        kind = CacheResultKind(state.kind)
        key = WorkKey(
            kind === PROCESSING_RESULT ? ITEM_PROCESS :
            kind === ITEM_ANALYSIS_RESULT ? ITEM_ANALYZE :
            kind === COLLECTION_PROCESS_RESULT ? COLLECTION_PROCESS : COLLECTION_ANALYZE,
            state.entity,
        )
        valid = lock(workspace.work.lock) do
            node = get(workspace.work.nodes, key, nothing)
            node !== nothing && node.state === :failed
        end
        valid || continue
        workspace.index.analysis_errors[state.entity] =
            something(state.message, "cached work failed")
    end
    if workspace.background_processing
        missing_analysis = lock(workspace.work.lock) do
            ItemRecord[
                record for record in workspace.index.hierarchy.all_items
                if !haskey(workspace.work.nodes, WorkKey(ITEM_ANALYZE, record.id))
            ]
        end
        for record in missing_analysis
            enqueue_processing!(workspace, record)
            enqueue_item_analysis!(workspace, record)
        end
    end
    return nothing
end

"""
Seed work graph completion nodes from a loaded cache index, each at counter revision 1.

Reopen validity is derived from position downstream of unchanged sources, not from a stored per-result
claim: everything seeds valid here, and `scan_source!`'s source-fingerprint diff bumps and re-marks
whatever changed while the workspace was closed.
"""
function seed_cached_work_nodes!(
    workspace::Workspace,
    index::ProjectCacheIndex,
)::Nothing
    records = Dict(record.id => record for record in workspace.index.hierarchy.all_items)
    for (result_key, state) in index.result_states
        work_state = CacheResultStatus(state.status) === RESULT_READY ? :ready : :failed
        kind = result_key.kind
        if kind === PROCESSING_RESULT || kind === ITEM_ANALYSIS_RESULT
            haskey(records, result_key.entity) || continue
            seed_work_node!(
                workspace,
                WorkKey(kind === PROCESSING_RESULT ? ITEM_PROCESS : ITEM_ANALYZE,
                    result_key.entity),
                UInt64(1),
                work_state,
            )
        else
            path = collection_path_tuple(result_key.entity)
            haskey(workspace.index.hierarchy.index, path) || continue
            seed_work_node!(
                workspace,
                WorkKey(kind === COLLECTION_PROCESS_RESULT ? COLLECTION_PROCESS :
                    COLLECTION_ANALYZE, result_key.entity),
                UInt64(1),
                work_state,
            )
        end
    end
    return nothing
end

"""
Invalidate the results of items and collections under changed collection-metadata subtrees at reopen.

The changed collection keys name the metadata rows the source reassigned while closed; every item at
or below one of those paths, and every ancestor collection, has its work bumped and re-marked so the
event graph recomputes them, matching an in-session metadata change.
"""
function invalidate_changed_metadata!(
    workspace::Workspace,
    scan_id::Int,
    changed_collections::Vector{String},
)::Nothing
    lock(workspace.publish_lock) do
        scan_id == workspace.scan.id || return
        prefixes = Vector{String}[collect(collection_path_tuple(key))
            for key in changed_collections]
        affected_records = ItemRecord[
            record for record in workspace.index.hierarchy.all_items
            if any(prefix -> length(prefix) <= length(record.collection) &&
                record.collection[1:length(prefix)] == prefix, prefixes)
        ]
        invalidate_records_work!(workspace, affected_records)
    end
    return nothing
end

"""Mark invalidated collection graph nodes stale by bumping their process/analyze revisions."""
function mark_collection_work_missing!(
    workspace::Workspace,
    collection_keys::Vector{String},
)::Nothing
    for collection_key in unique(collection_keys)
        path = collection_path_tuple(collection_key)
        haskey(workspace.index.hierarchy.index, path) || continue
        bump_work_missing!(workspace, WorkKey(COLLECTION_PROCESS, collection_key))
        bump_work_missing!(workspace, WorkKey(COLLECTION_ANALYZE, collection_key))
    end
    return nothing
end

"""Return the current node if the completion revision still owns the result."""
function current_completion_node(
    workspace::Workspace,
    key::WorkKey,
    revision::UInt64,
)::Union{Nothing,WorkNode}
    return lock(workspace.work.lock) do
        node = get(workspace.work.nodes, key, nothing)
        node === nothing && return nothing
        node.revision == revision && node.state === :running || return nothing
        return node
    end
end

"""Finish a current work node and wake selected waiters."""
function finish_work_node!(
    workspace::Workspace,
    node::WorkNode,
    state::Symbol,
)::Vector{Channel{Any}}
    return lock(workspace.work.lock) do
        current = get(workspace.work.nodes, node.key, nothing)
        current === node || return Channel{Any}[]
        node.state = state
        workspace.work.completed += 1
        wake_ready_dependents!(workspace.work, node.key)
        waiters = copy(node.waiters)
        empty!(node.waiters)
        waiters
    end
end

"""Queue item analysis after processing succeeds."""
function enqueue_item_analysis!(workspace::Workspace, record::ItemRecord)::Nothing
    key = WorkKey(ITEM_ANALYZE, record.id)
    enqueue_work!(
        workspace,
        key,
        current_revision(workspace.work, key);
        priority=2,
        dependencies=WorkKey[WorkKey(ITEM_PROCESS, record.id)],
    )
    return nothing
end

"""
Queue a collection's terminal work when all member item-analyze nodes are terminal.

A collection with a registered `process` recipe queues COLLECTION_PROCESS; its COLLECTION_ANALYZE
then depends on that node, else directly on the member ITEM_ANALYZE nodes.
"""
function enqueue_ready_collection_work!(
    workspace::Workspace,
    collection_keys::Vector{String},
)::Nothing
    for collection_key in collection_keys
        path = collection_path_tuple(collection_key)
        node = get(workspace.index.hierarchy.index, path, nothing)
        node === nothing && continue
        isempty(node.items) && continue
        member_kinds = unique(record.kind for record in node.items)
        member_analyze = WorkKey[WorkKey(ITEM_ANALYZE, record.id) for record in node.items]
        has_process = any(k -> has_collection_process(workspace.project, k), member_kinds)
        has_analyze = any(k -> has_collection_analysis(workspace.project, k), member_kinds)
        if has_process
            process_key = WorkKey(COLLECTION_PROCESS, collection_key)
            enqueue_work!(
                workspace,
                process_key,
                current_revision(workspace.work, process_key);
                priority=1,
                dependencies=member_analyze,
            )
        end
        has_analyze || continue
        analyze_key = WorkKey(COLLECTION_ANALYZE, collection_key)
        enqueue_work!(
            workspace,
            analyze_key,
            current_revision(workspace.work, analyze_key);
            priority=1,
            dependencies=has_process ?
                WorkKey[WorkKey(COLLECTION_PROCESS, collection_key)] : member_analyze,
        )
    end
    return nothing
end

"""
Apply one successful work completion to the WorkspaceIndex.

Runs under the publish lock and stays cheap: the worker already persisted its own payload before
completing, and every cache mutation reachable from here is a non-blocking buffered enqueue.
"""
function publish_work_success!(
    workspace::Workspace,
    node::WorkNode,
    result,
)::Nothing
    key = node.key
    if key.kind === SOURCE_INTERPRET
        records = result.records::Vector{ItemRecord}
        old_records, invalidated = Profiling.@profile_span workspace.profiler :source :publish_records Profiling.ProfileAttributes(
            source_id=key.entity,
            batch_size=length(records),
        ) begin
            publish_source_item_records!(workspace, key.entity, records)
        end
        delete_collection_metadata!(workspace.cache.db, invalidated)
        mark_collection_work_missing!(workspace, invalidated)
        isempty(result.conflicts) ||
            (workspace.index.analysis_errors[key.entity] = join(result.conflicts, "; "))
        # One lock acquisition for the whole batch: the per-record calls re-enter it cheaply, and
        # each newly queued node wakes exactly one worker instead of the pool per record.
        lock(workspace.work.lock) do
            for record in records
                # Re-interpretation replaced the source data and deleted the old payloads, so the
                # item's process and analyze results are stale: bump both off their prior revisions
                # (a still-:ready ITEM_PROCESS would otherwise answer from a deleted payload).
                bump_work_missing!(workspace, WorkKey(ITEM_PROCESS, record.id))
                bump_work_missing!(
                    workspace,
                    WorkKey(ITEM_ANALYZE, record.id);
                    dependencies=WorkKey[WorkKey(ITEM_PROCESS, record.id)],
                )
                workspace.background_processing && enqueue_processing!(workspace, record)
            end
        end
    elseif key.kind === ITEM_PROCESS
        # The worker stored the processed payload before completing, so item analysis becomes
        # runnable the moment this node finishes and late joiners reading the cache find the data.
        record = get(workspace.index.items, key.entity, nothing)
        record === nothing || enqueue_item_analysis!(workspace, record)
    elseif key.kind === ITEM_ANALYZE
        record = workspace.index.items[key.entity]
        publish_item_metadata_layer!(workspace, record, metadata_dict(result))
        enqueue_ready_collection_work!(workspace, affected_collection_keys([record]))
    elseif key.kind === COLLECTION_PROCESS
        # Members' rewritten payloads and metadata were persisted worker-side; publish their layers.
        for output in result.outputs
            member_id = id(output)
            record = get(workspace.index.items, member_id, nothing)
            record === nothing && continue
            publish_item_metadata_layer!(workspace, record, result.metadata_by_id[member_id])
        end
    else
        path = collection_path_tuple(key.entity)
        hierarchy_node = get(workspace.index.hierarchy.index, path, nothing)
        hierarchy_node === nothing && return nothing
        analysis = metadata_dict(result)
        empty!(hierarchy_node.analysis)
        merge!(hierarchy_node.analysis, analysis)
        dropped = store_collection_metadata!(workspace.cache.db, key.entity, analysis)
        isempty(dropped) ||
            (workspace.index.analysis_errors[key.entity] = join(dropped, "; "))
    end
    return nothing
end

"""
Publish one item's computed metadata layer to the index and cache, surfacing any wide-cache type
conflict through the analysis-error channel.
"""
function publish_item_metadata_layer!(
    workspace::Workspace,
    record::ItemRecord,
    computed,
)::Nothing
    workspace.index.item_metadata[record.id] = Dict{Symbol,Any}(computed)
    dropped = store_item_metadata!(
        workspace.cache.db, record, delivered_metadata(workspace, record))
    isempty(dropped) && return nothing
    workspace.index.analysis_errors[record.id] = join(dropped, "; ")
    return nothing
end

"""Apply one failed work completion to WorkspaceIndex and ProjectCache."""
function publish_work_failure!(
    workspace::Workspace,
    node::WorkNode,
    failure::CapturedException,
)::Bool
    key = node.key
    if is_job_cancelled(failure.ex)
        # Cancellation is not an analysis error: finish the node and answer waiters with the
        # cancellation, but record nothing in the index or the cache failure results.
        waiters = finish_work_node!(workspace, node, :failed)
        for waiter in waiters
            put!(waiter, ProcessingResult(nothing, failure))
        end
        return false
    end
    message = sprint(showerror, failure.ex)
    if key.kind === SOURCE_INTERPRET
        old_records, invalidated = remove_source_item_output!(workspace, key.entity)
        delete_source_item!(workspace.cache.db, key.entity, old_records)
        delete_collection_metadata!(workspace.cache.db, invalidated)
        mark_collection_work_missing!(workspace, invalidated)
        workspace.index.analysis_errors[key.entity] = "interpret_source_item: " * message
        finish_work_node!(workspace, node, :failed)
        return true
    end
    source_item_id_value = ""
    record = key.kind in (ITEM_PROCESS, ITEM_ANALYZE) ?
        get(workspace.index.items, key.entity, nothing) : nothing
    record === nothing || (source_item_id_value = record.source_item_id)
    store_result_failure!(
        workspace.cache.db,
        key.kind === ITEM_PROCESS ? PROCESSING_RESULT :
        key.kind === ITEM_ANALYZE ? ITEM_ANALYSIS_RESULT :
        key.kind === COLLECTION_PROCESS ? COLLECTION_PROCESS_RESULT : COLLECTION_ANALYSIS_RESULT,
        key.entity,
        source_item_id_value,
        message,
    )
    workspace.index.analysis_errors[key.entity] = message
    waiters = finish_work_node!(workspace, node, :failed)
    for waiter in waiters
        put!(waiter, ProcessingResult(nothing, failure))
    end
    if key.kind === ITEM_ANALYZE && record !== nothing
        enqueue_ready_collection_work!(workspace, affected_collection_keys([record]))
    end
    return true
end

"""Commit one finished work revision the moment its worker completes; stale revisions are ignored."""
function publish_work_completion!(
    workspace::Workspace,
    key::WorkKey,
    revision::UInt64,
    result,
)::Nothing
    lock(workspace.publish_lock) do
        node = current_completion_node(workspace, key, revision)
        node === nothing && return
        if result isa CapturedException
            publish_work_failure!(workspace, node, result)
        else
            publish_work_success!(workspace, node, result)
            waiters = finish_work_node!(workspace, node, :ready)
            item = key.kind === ITEM_PROCESS ? result.item : nothing
            for waiter in waiters
                put!(waiter, ProcessingResult(item, nothing))
            end
        end
        finish_publish!(workspace)
    end
    return nothing
end

"""Advance idle-gated work and wake blocked waiters after one publication (publish lock held)."""
function finish_publish!(workspace::Workspace)::Nothing
    if workspace.profile_restart_pending && !engine_work_running(workspace)
        reset_work_graph!(workspace)
        workspace.profile_restart_pending = false
        workspace.profiler.state = :idle
        Profiling.start!(
            workspace.profiler,
            () -> profile_counter_snapshot(workspace),
        )
        scan_source!(workspace; rebuild=true)
    end
    notify(workspace.idle_condition; all=true)
    return nothing
end

"""
Block until the workspace has no running scan and no queued or running work, or until `timeout`
seconds pass. Purely event-driven: every publication notifies the workspace idle condition. Cache
writes still flushing do not delay idleness — reads overlay the write buffer, and closing the
workspace drains it.
"""
function wait_workspace_idle!(workspace::Workspace; timeout::Real=60)::Workspace
    deadline = time() + Float64(timeout)
    lock(workspace.publish_lock) do
        while engine_work_running(workspace)
            wait_condition_deadline(workspace.idle_condition, deadline) || break
        end
    end
    return workspace
end

"""Publish the loaded (or missing) cache index as the instant first view of one scan."""
function publish_cache_state!(
    workspace::Workspace,
    scan_id::Int,
    state::Symbol,
    index::Union{Nothing,ProjectCacheIndex},
)::Nothing
    lock(workspace.publish_lock) do
        scan_id == workspace.scan.id || return
        if index === nothing
            replace_item_index!(workspace, Hierarchy(
                source_id(workspace.source),
                has_collection_metadata(workspace.source),
            ))
            empty!(workspace.index.analysis_errors)
        else
            apply_cache_index!(workspace, index)
        end
        workspace.cache_state = state
        state == :missing && workspace.cache.operation != :rebuild &&
            (workspace.cache.operation = :build)
        finish_publish!(workspace)
    end
    return nothing
end

"""Publish one scan's discovered source changes and queue their interpretation."""
function publish_source_changes!(
    workspace::Workspace,
    scan_id::Int,
    changes::SourceChanges,
    fingerprints::AbstractDict{String},
    status::ProjectCacheStatus,
)::Nothing
    lock(workspace.publish_lock) do
        scan_id == workspace.scan.id || return
        apply_source_changes!(workspace, changes, fingerprints, status)
        workspace.cache_state = :ready
        finish_publish!(workspace)
    end
    return nothing
end

"""Publish one source-owned watcher update through the same path as a scan; clears any source error."""
function publish_source_event!(
    workspace::Workspace,
    changes::SourceChanges,
)::Nothing
    lock(workspace.publish_lock) do
        workspace.closed && return
        workspace.source_error = ""
        fingerprints = workspace.index.source === nothing ?
            Dict{String,Any}() :
            workspace.index.source.source_item_fingerprints
        apply_source_changes!(workspace, changes, fingerprints, nothing)
        finish_publish!(workspace)
    end
    return nothing
end

"""Publish one recoverable source failure; the next good update clears it."""
function publish_source_event!(
    workspace::Workspace,
    error::SourceError,
)::Nothing
    lock(workspace.publish_lock) do
        workspace.closed && return
        workspace.source_error = error.message
        finish_publish!(workspace)
    end
    return nothing
end

"""
Publish the end of one scan job: finished normally, canceled, or failed.

This is the single settlement point for `scan.state`: the scan task reaches it on every exit path,
and it assigns a terminal state unconditionally, so a cancel request that lands between the last
progress publication and this one can never leave the job stuck in `:canceling`.
"""
function publish_scan_end!(
    workspace::Workspace,
    scan_id::Int;
    cache_hit::Bool=false,
    canceled::Bool=false,
    error=nothing,
    backtrace=nothing,
)::Nothing
    lock(workspace.publish_lock) do
        scan_id == workspace.scan.id || return
        # A scan that ends without reaching :ready must settle cache_state, or busy spins forever.
        if canceled
            workspace.scan.state = :canceled
            workspace.cache_state in (:loading, :writing, :canceling) &&
                (workspace.cache_state = :stale)
        elseif error !== nothing
            workspace.scan.state = :error
            workspace.scan.error =
                "Source scan failed. See the console for full details."
            workspace.cache_state in (:loading, :writing, :canceling) &&
                (workspace.cache_state = :error)
            @error(
                "Source scan job failed",
                source=source_id(workspace.source),
                exception=(error, backtrace),
            )
        else
            workspace.scan.state = cache_hit ? :unchanged : :done
        end
        workspace.scan.cancel_token = nothing
        finish_publish!(workspace)
    end
    return nothing
end

"""
Rebuild the GUI-facing status snapshot when work is or was just running.

Display-only: the engine publishes without it, and idle frames keep reading the cached value.
"""
function refresh_status!(workspace::Workspace)::Nothing
    (workspace_busy(workspace) || workspace.status.busy) &&
        (workspace.status = workspace_status(workspace))
    return nothing
end
