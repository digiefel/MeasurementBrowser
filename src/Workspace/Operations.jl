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

"""Whether the scan or any queued or running work-graph node is still in flight."""
function engine_work_running(workspace::Workspace)::Bool
    return source_scan_running(workspace) || cache_work_running(workspace) ||
        lock(workspace.work.lock) do
            workspace.work.active > 0
        end
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
        workspace.status_dirty[] = true
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
Diff cached source metadata against the open source; update changed cache rows.

Returns items whose work should be invalidated. Cache `read` is the baseline (disk plus buffer
overlay). Rows absent from the baseline are persisted without invalidation.

Reconciles the whole index by default. Pass `collections` (a set of collection path keys) to diff
only those collections and skip the whole-index item pass: the streaming publish path uses this,
because interpret workers already persist each item's metadata and the caller invalidates the
touched collections' work itself. Diffing every collection and item on every publish is what makes
an otherwise linear scan quadratic.
"""
function reconcile_source_metadata_cache!(
    workspace::Workspace;
    collections::Union{Nothing,Vector{String}}=nothing,
    refresh_hierarchy::Bool=false,
)::Vector{ItemRecord}
    old_hierarchy = nothing
    if refresh_hierarchy
        old_hierarchy = workspace.index.hierarchy
        records = sort!(
            collect(values(workspace.index.items));
            by=record -> (record.source_item_id, record.id),
        )
        new_hierarchy = Hierarchy(records, workspace.source, old_hierarchy.skipped_count)
        old_node_analysis = Dict(path => copy(node.analysis) for (path, node) in old_hierarchy.index)
        for (path, analysis) in old_node_analysis
            node = get(new_hierarchy.index, path, nothing)
            node === nothing || merge!(node.analysis, analysis)
        end
        workspace.index.hierarchy = new_hierarchy
        refresh_collection_metadata_keys!(workspace.index)
    end

    cachedb = workspace.cache.db
    if cachedb isa CacheDB
        stale = ItemRecord[]
        baseline_collections = read(cachedb.source_collection_metadata)
        collection_paths = collections === nothing ?
            [collect(path) for path in keys(workspace.index.hierarchy.index)] :
            [collect(collection_path_tuple(key)) for key in unique(collections)]
        for path in collection_paths
            key = collection_path_key(path)
            current = metadata_dict(collection_metadata(workspace.source, path))
            old = get(baseline_collections, key, nothing)
            if old === nothing
                isempty(current) && continue
                edit_source_collection_metadata!(cachedb, key, current)
            elseif !isequal(old, current)
                append!(stale, ItemRecord[
                    record for record in all_items(workspace.index.hierarchy)
                    if length(path) <= length(record.collection) &&
                       record.collection[1:length(path)] == path
                ])
                edit_source_collection_metadata!(cachedb, key, current)
            end
        end

        # Item metadata is reconciled only on a whole-index pass; the scoped publish path relies on
        # interpret workers having persisted each item's row already.
        if collections === nothing
            key_to_id = Dict(row.item_key => row.id for row in values(read(cachedb.items)))
            id_to_key = Dict(id => key for (key, id) in key_to_id)
            baseline_items = read(cachedb.source_item_metadata)
            for record in values(workspace.index.items)
                item_key = get(id_to_key, record.id, nothing)
                item_key === nothing && continue
                current = record.metadata
                old = get(baseline_items, item_key, nothing)
                if old === nothing
                    isempty(current) && continue
                    edit_source_item_metadata!(cachedb, item_key, current)
                elseif !isequal(old, current)
                    push!(stale, record)
                    edit_source_item_metadata!(cachedb, item_key, current)
                end
            end
        end
        return unique(stale)
    end

    old_hierarchy === nothing && return ItemRecord[]
    records = collect(values(workspace.index.items))
    return unique(ItemRecord[
        record for record in records
        if effective_metadata(old_hierarchy, record) !=
           effective_metadata(workspace.index.hierarchy, record)
    ])
end

"""
Bump and re-enqueue the work of a set of records and their collections after their inputs changed.

Shared by the in-session metadata event and the reopen source-metadata diff: the full dependent
subtree is enqueued immediately so live nodes gate delivery while fresh work runs.
"""
function invalidate_records_work!(workspace::Workspace, changed::Vector{ItemRecord})::Nothing
    for record in changed
        delete!(workspace.index.item_metadata, record.id)
        delete!(workspace.index.analysis_errors, record.id)
    end
    enqueue_dependent_subtree!(workspace, changed)
    return nothing
end

"""
Invalidate one record's process/analyze steps and every affected collection step, scheduling their
recomputation only when background processing is enabled.

Invalidation, source removal, and interpret failure all fan out through this helper. The stale cache
result-state is always cleared so the next read (on-demand or background) recomputes; the fresh work
is only enqueued when `background_processing` is on. Without it, processing stays selection-driven and
a live scan interprets the tree without eagerly processing or analyzing every item.
"""
function enqueue_dependent_subtree!(
    workspace::Workspace,
    records::Vector{ItemRecord};
    collection_keys::Vector{String}=affected_collection_keys(records),
    priority::Int=2,
)::Nothing
    background = workspace.background_processing
    for collection_key in unique(collection_keys)
        path = collection_path_tuple(collection_key)
        hierarchy_node = get(workspace.index.hierarchy.index, path, nothing)
        hierarchy_node === nothing && continue
        empty!(hierarchy_node.analysis)
    end
    lock(workspace.work.lock) do
        for record in records
            process_key = WorkKey(ITEM_PROCESS, record.id)
            analyze_key = WorkKey(ITEM_ANALYZE, record.id)
            clear_work_result_state!(workspace, process_key)
            clear_work_result_state!(workspace, analyze_key)
            background || continue
            enqueue_work!(
                workspace, process_key, bump_revision!(workspace.work, process_key);
                priority,
            )
            enqueue_work!(
                workspace, analyze_key, bump_revision!(workspace.work, analyze_key);
                priority,
                dependencies=WorkKey[process_key],
            )
        end
        if background
            enqueue_collection_work!(workspace, collection_keys; supersede=true)
        else
            for collection_key in unique(collection_keys)
                clear_work_result_state!(workspace, WorkKey(COLLECTION_PROCESS, collection_key))
                clear_work_result_state!(workspace, WorkKey(COLLECTION_ANALYZE, collection_key))
            end
        end
    end
    return nothing
end

"""Refresh the authoritative published source snapshot from current index state."""
function refresh_workspace_source!(workspace::Workspace)::Nothing
    workspace.index.source = SourceScan(
        source_id(workspace.source),
        source_label(workspace.source),
        workspace.index.hierarchy,
        ItemFailure[],
    )
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
function ingest_source_changes!(
    workspace::Workspace,
    changes::SourceChanges,
    status::Union{Nothing,ProjectCacheStatus},
)::Bool
    changed = false
    for removed in changes.removals
        old_records, invalidated = remove_source_item_output!(workspace, removed)
        delete_source_item!(workspace.cache.db, removed, old_records)
        delete_collection_metadata!(workspace.cache.db, invalidated)
        enqueue_dependent_subtree!(workspace, old_records; collection_keys=invalidated)
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
    if changes.metadata_changed
        stale = reconcile_source_metadata_cache!(workspace; refresh_hierarchy=true)
        isempty(stale) || invalidate_records_work!(workspace, stale)
        changed = changed || !isempty(stale)
    end
    refresh_workspace_source!(workspace)
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
        job.discovered[] = 0
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
        workspace.status_dirty[] = true
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
                    source_items(
                        workspace.source;
                        on_progress=count -> begin
                            job.discovered[] = count
                            workspace.status_dirty[] = true
                        end,
                    )
                end
            end
            # The cache's `source_items` table is the sole home of previous-session fingerprints;
            # a memory-only cache holds nothing, so every discovered item re-interprets.
            previous = _load_source_item_fingerprints(cachedb)
            current = Dict{String,Any}()
            upserts = AbstractDataSourceItem[]
            seen = Set{String}()
            for item in discovered
                id_value = source_item_id(item)
                push!(seen, id_value)
                current_fingerprint = fingerprint(item)
                current[id_value] = current_fingerprint
                # A `nothing` fingerprint means the source cannot prove the item is unchanged,
                # so it is always re-read.
                if !haskey(previous, id_value) || current_fingerprint === nothing ||
                   !isequal(previous[id_value], current_fingerprint)
                    push!(upserts, item)
                end
            end
            removals = String[id for id in keys(previous) if !(id in seen)]
            stale = count(
                id_value -> haskey(previous, id_value) && !isequal(previous[id_value], current[id_value]),
                keys(current),
            )
            new_items = count(id_value -> !haskey(previous, id_value), keys(current))
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
        if engine_work_running(workspace)
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

"""Replace the workspace item index with one complete hierarchy."""
function replace_item_index!(
    workspace::Workspace,
    hierarchy::Hierarchy,
)::Nothing
    items = Dict{String,ItemRecord}()
    sizehint!(items, length(all_items(hierarchy)))
    items_by_source = Dict{String,Vector{String}}()
    metadata_keys = Set{Symbol}()
    for item in all_items(hierarchy)
        haskey(items, item.id) && error("Duplicate item id generated during scan: $(item.id)")
        items[item.id] = item
        push!(get!(() -> String[], items_by_source, item.source_item_id), item.id)
    end
    for node in values(hierarchy.index)
        union!(metadata_keys, keys(node.metadata))
    end
    item_metadata = Dict{String,Dict{Symbol,Any}}()
    for item in all_items(hierarchy)
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

"""Whether one loaded cache row is `RESULT_READY`."""
function cache_index_ready(
    index::ProjectCacheIndex,
    kind::CacheResultKind,
    entity::AbstractString,
)::Bool
    state = get(index.result_states, CacheResultKey(kind, String(entity)), nothing)
    state !== nothing && CacheResultStatus(state.status) === RESULT_READY
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
        all_items(index.source.hierarchy),
        workspace.source,
        index.source.hierarchy.skipped_count,
    )
    replace_item_index!(workspace, hierarchy)
    # The cache restores the computed metadata layer; source-table diff before seed skips stale work.
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
        hierarchy,
        index.source.analysis_failures,
    )

    workspace.index.analysis_errors = copy(index.analysis_errors)
    stale = reconcile_source_metadata_cache!(workspace)
    isempty(stale) || invalidate_records_work!(workspace, stale)
    for state in values(index.result_states)
        CacheResultStatus(state.status) === RESULT_FAILED || continue
        workspace.index.analysis_errors[state.entity] =
            something(state.message, "cached work failed")
    end
    if workspace.background_processing
        for record in all_items(workspace.index.hierarchy)
            cache_index_ready(index, PROCESSING_RESULT, record.id) ||
                enqueue_processing!(workspace, record)
            cache_index_ready(index, ITEM_ANALYSIS_RESULT, record.id) ||
                enqueue_item_analysis!(workspace, record)
        end
        for (path, hierarchy_node) in workspace.index.hierarchy.index
            isempty(hierarchy_node.items) && continue
            collection_key = collection_path_key(collect(String, path))
            member_kinds = unique(record.kind for record in hierarchy_node.items)
            if any(k -> has_collection_process(workspace.project, k), member_kinds) &&
                    !cache_index_ready(index, COLLECTION_PROCESS_RESULT, collection_key)
                enqueue_collection_work!(workspace, [collection_key]; supersede=false)
            elseif any(k -> has_collection_analysis(workspace.project, k), member_kinds) &&
                    !cache_index_ready(index, COLLECTION_ANALYSIS_RESULT, collection_key)
                enqueue_collection_work!(workspace, [collection_key]; supersede=false)
            end
        end
    end
    return nothing
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
Enqueue process/analyze work for a set of collections, gated on their members' item-analyze nodes.

A collection with a registered `process` recipe queues COLLECTION_PROCESS; its COLLECTION_ANALYZE
then depends on that node, else directly on the member ITEM_ANALYZE nodes. Cached result rows are
cleared so a fresh fold recomputes. `supersede` picks the revision policy — the one axis that
separates the two callers: `false` joins existing work at its current revision (a member finished, or
reopen gap-fill); `true` bumps to a fresh revision to replace work running on stale inputs
(invalidation).
"""
function enqueue_collection_work!(
    workspace::Workspace,
    collection_keys;
    supersede::Bool,
)::Nothing
    revision(key) = supersede ?
        bump_revision!(workspace.work, key) : current_revision(workspace.work, key)
    for collection_key in unique(collection_keys)
        path = collection_path_tuple(collection_key)
        hierarchy_node = get(workspace.index.hierarchy.index, path, nothing)
        hierarchy_node === nothing && continue
        isempty(hierarchy_node.items) && continue
        member_kinds = unique(record.kind for record in hierarchy_node.items)
        member_analyze = WorkKey[
            WorkKey(ITEM_ANALYZE, record.id) for record in hierarchy_node.items]
        has_process = any(k -> has_collection_process(workspace.project, k), member_kinds)
        has_analyze = any(k -> has_collection_analysis(workspace.project, k), member_kinds)
        process_key = WorkKey(COLLECTION_PROCESS, collection_key)
        analyze_key = WorkKey(COLLECTION_ANALYZE, collection_key)
        if !supersede
            process_done = !has_process ||
                cache_work_status(workspace, process_key) === :ready
            analyze_done = !has_analyze ||
                cache_work_status(workspace, analyze_key) === :ready
            already_live = lock(workspace.work.lock) do
                get(workspace.work.nodes, process_key, nothing) !== nothing ||
                get(workspace.work.nodes, analyze_key, nothing) !== nothing
            end
            process_done && analyze_done && !already_live && continue
        end
        if has_process
            supersede && clear_work_result_state!(workspace, process_key)
            enqueue_work!(
                workspace, process_key, revision(process_key);
                priority=1, dependencies=member_analyze,
            )
        end
        has_analyze || continue
        supersede && clear_work_result_state!(workspace, analyze_key)
        enqueue_work!(
            workspace, analyze_key, revision(analyze_key);
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
        enqueue_dependent_subtree!(workspace, records; collection_keys=invalidated)
        publish_metadata_conflicts!(workspace, key.entity, key.kind, result.conflicts)
        reconcile_source_metadata_cache!(workspace; collections=invalidated)
    elseif key.kind === ITEM_PROCESS
        # The worker stored the processed payload before completing, so item analysis becomes
        # runnable the moment this node finishes and late joiners reading the cache find the data.
        record = get(workspace.index.items, key.entity, nothing)
        if record !== nothing
            clear_work_result_state!(workspace, WorkKey(ITEM_ANALYZE, record.id))
            enqueue_item_analysis!(workspace, record)
        end
    elseif key.kind === ITEM_ANALYZE
        record = workspace.index.items[key.entity]
        publish_item_metadata_layer!(workspace, record, metadata_dict(result))
        enqueue_collection_work!(workspace, affected_collection_keys([record]); supersede=false)
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
        publish_metadata_conflicts!(workspace, key.entity, key.kind, dropped)
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
    publish_metadata_conflicts!(workspace, record.id, ITEM_ANALYZE, dropped)
    return nothing
end

function publish_metadata_conflicts!(
    workspace::Workspace,
    entity::AbstractString,
    stage,
    conflicts::Vector{String},
)::Nothing
    isempty(conflicts) && return nothing
    message = join(conflicts, "; ")
    workspace.index.analysis_errors[String(entity)] = message
    @warn(
        "Workspace metadata conflict",
        entity=String(entity),
        stage=stage,
        message=message,
    )
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
        waiters = finish_work_node!(workspace, node)
        for waiter in waiters
            put!(waiter, ProcessingResult(nothing, failure))
        end
        return false
    end
    message = sprint(showerror, failure.ex)
    record = key.kind in (ITEM_PROCESS, ITEM_ANALYZE) ?
        get(workspace.index.items, key.entity, nothing) : nothing
    if key.kind === SOURCE_INTERPRET
        old_records, invalidated = remove_source_item_output!(workspace, key.entity)
        delete_source_item!(workspace.cache.db, key.entity, old_records)
        delete_collection_metadata!(workspace.cache.db, invalidated)
        enqueue_dependent_subtree!(workspace, old_records; collection_keys=invalidated)
        workspace.index.analysis_errors[key.entity] = "interpret_source_item: " * message
        store_source_item_failure!(workspace.cache.db, key.entity, message)
        @error(
            "Source interpretation failed",
            source_item=key.entity,
            stage=key.kind,
            exception=(failure.ex, failure.processed_bt),
        )
    else
        source_item_id_value = record === nothing ? "" : record.source_item_id
        store_result_failure!(
            workspace.cache.db,
            _work_key_cache_kind(key),
            key.entity,
            source_item_id_value,
            message,
        )
        workspace.index.analysis_errors[key.entity] = message
        @warn(
            "Workspace work failed",
            entity=key.entity,
            stage=key.kind,
            exception=(failure.ex, failure.processed_bt),
        )
        key.kind === ITEM_ANALYZE && record !== nothing &&
            enqueue_collection_work!(workspace, affected_collection_keys([record]); supersede=false)
    end
    waiters = finish_work_node!(workspace, node)
    for waiter in waiters
        put!(waiter, ProcessingResult(nothing, failure))
    end
    return true
end

"""Commit one finished work revision the moment its worker completes; stale revisions are ignored."""
function publish_work_completion!(
    workspace::Workspace,
    node::WorkNode,
    result,
)::Nothing
    key = node.key
    lock(workspace.publish_lock) do
        work_node_current(workspace, node) || return
        if result isa CapturedException
            publish_work_failure!(workspace, node, result)
        else
            publish_work_success!(workspace, node, result)
            waiters = finish_work_node!(workspace, node)
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
    workspace.status_dirty[] = true
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
    status::ProjectCacheStatus,
)::Nothing
    lock(workspace.publish_lock) do
        scan_id == workspace.scan.id || return
        ingest_source_changes!(workspace, changes, status)
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
        ingest_source_changes!(workspace, changes, nothing)
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

"""Rebuild the GUI-facing status snapshot after engine publications or discovery progress."""
function refresh_status!(workspace::Workspace)::Nothing
    dirty = Base.Threads.atomic_xchg!(workspace.status_dirty, false)
    dirty && (workspace.status = workspace_status(workspace))
    return nothing
end
