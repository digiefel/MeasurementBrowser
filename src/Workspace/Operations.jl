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
    cancel_job!(workspace.scan)
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
    work_kind_running(workspace, COLLECTION_STATS)

"""Whether any item is waiting for or running processing, or its result is still being written."""
function processing_work_running(workspace::Workspace)::Bool
    return work_kind_running(workspace, PROCESS_ITEM) ||
        work_kind_running(workspace, ITEM_STATS) ||
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

"""Request cancellation of one running workspace operation."""
function cancel_job!(job::WorkspaceJob)::Nothing
    job.cancel_token === nothing && return nothing
    Base.Threads.atomic_xchg!(job.cancel_token, true)
    job.state = :canceling
    return nothing
end

cancel_scan!(workspace::Workspace)::Nothing = (cancel_job!(workspace.scan); nothing)
cancel_analysis!(workspace::Workspace)::Nothing = (cancel_waiting_work!(workspace); nothing)

"""
Cancel in-flight cache work.

Cache writes are now produced by the source scan itself, so cancelling cache work cancels the scan.
"""
cancel_cache!(workspace::Workspace)::Nothing = cancel_scan!(workspace)

"""Return records currently published for one source item."""
function source_item_records(index::WorkspaceIndex, source_item_id_value::String)::Vector{ItemRecord}
    return ItemRecord[
        record for record in values(index.items)
        if record.source_item_id == source_item_id_value
    ]
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

"""Rebuild hierarchy indexes from the currently published records."""
function rebuild_workspace_hierarchy!(workspace::Workspace)::Nothing
    records = sort!(collect(values(workspace.index.items)); by=record -> (record.source_item_id, record.id))
    hierarchy = Hierarchy(records, workspace.source)
    workspace.index.hierarchy = hierarchy
    parameter_keys = Set{Symbol}()
    for node in values(hierarchy.index)
        union!(parameter_keys, keys(node.parameters))
    end
    workspace.index.collection_parameter_keys = sort!(collect(parameter_keys); by=String)
    return nothing
end

"""Reapply source parameters and invalidate only items whose effective inputs changed."""
function apply_parameter_changes!(workspace::Workspace)::Bool
    old_hierarchy = workspace.index.hierarchy
    records = sort!(
        collect(values(workspace.index.items));
        by=record -> (record.source_item_id, record.id),
    )
    old_parameters = Dict(
        record.id => effective_parameters(old_hierarchy, record) for record in records
    )
    new_hierarchy = Hierarchy(records, workspace.source, old_hierarchy.skipped_count)
    changed = ItemRecord[
        record for record in records
        if old_parameters[record.id] != effective_parameters(new_hierarchy, record)
    ]
    old_node_stats = Dict(path => copy(node.stats) for (path, node) in old_hierarchy.index)
    for (path, stats) in old_node_stats
        node = get(new_hierarchy.index, path, nothing)
        node === nothing || merge!(node.stats, stats)
    end
    workspace.index.hierarchy = new_hierarchy
    parameter_keys = Set{Symbol}()
    for node in values(new_hierarchy.index)
        union!(parameter_keys, keys(node.parameters))
    end
    workspace.index.collection_parameter_keys = sort!(collect(parameter_keys); by=String)

    isempty(changed) && return false
    invalidate_item_results!(workspace.cache.db, changed)
    affected = affected_collection_keys(changed)
    for record in changed
        delete!(workspace.index.item_stats, record.id)
        delete!(workspace.index.analysis_errors, record.id)
        mark_work_missing!(
            workspace,
            WorkKey(PROCESS_ITEM, record.id),
            item_revision(workspace, record),
        )
        mark_work_missing!(
            workspace,
            WorkKey(ITEM_STATS, record.id),
            item_revision(workspace, record);
            dependencies=WorkKey[WorkKey(PROCESS_ITEM, record.id)],
        )
        workspace.background_processing && enqueue_processing!(workspace, record)
    end
    delete_collection_stats!(workspace.cache.db, affected)
    for key in affected
        node = get(new_hierarchy.index, collection_path_tuple(key), nothing)
        node === nothing || empty!(node.stats)
    end
    mark_collection_stats_missing!(workspace, affected)
    return true
end

"""Refresh the authoritative published source snapshot from current index state."""
function refresh_workspace_source!(
    workspace::Workspace,
    fingerprints::AbstractDict{String},
)::Nothing
    workspace.index.source = SourceScan(
        source_id(workspace.source),
        source_label(workspace.source),
        Dict{String,Any}(fingerprints),
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
        delete!(workspace.index.item_stats, record.id)
        delete!(workspace.index.analysis_errors, record.id)
    end
    delete!(workspace.index.analysis_errors, source_item_id_value)
    filter!(id -> haskey(workspace.index.items, id), workspace.selection.item_ids)
    rebuild_workspace_hierarchy!(workspace)
    for key in invalidated
        node = get(workspace.index.hierarchy.index, collection_path_tuple(key), nothing)
        node === nothing || empty!(node.stats)
    end
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
    old_records, _ = remove_source_item_output!(workspace, source_item_id_value)
    for record in records
        workspace.index.items[record.id] = record
        isempty(record.stats) || (workspace.index.item_stats[record.id] = copy(record.stats))
    end
    filter!(id -> haskey(workspace.index.items, id), workspace.selection.item_ids)
    rebuild_workspace_hierarchy!(workspace)
    invalidated = unique([old_keys; new_keys])
    for key in invalidated
        node = get(workspace.index.hierarchy.index, collection_path_tuple(key), nothing)
        node === nothing || empty!(node.stats)
    end
    return old_records, invalidated
end

"""Submit one source-change batch to the work graph and remove deleted published output."""
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
        delete_collection_stats!(workspace.cache.db, invalidated)
        mark_collection_stats_missing!(workspace, invalidated)
        lock(workspace.work.lock) do
            delete!(workspace.work.source_items, removed)
        end
        changed = true
    end
    for source_item in changes.upserts
        source_item_id_value = source_item_id(source_item)
        old_records = source_item_records(workspace.index, source_item_id_value)
        delete_source_item!(workspace.cache.db, source_item_id_value, old_records)
        revision = UInt64(hash((source_item_id_value, fingerprint(source_item))))
        lock(workspace.work.lock) do
            workspace.work.source_items[source_item_id_value] = source_item
        end
        enqueue_work!(
            workspace,
            WorkKey(INTERPRET_SOURCE, source_item_id_value),
            revision;
            priority=3,
        )
        changed = true
    end
    changes.parameters_changed && (changed = apply_parameter_changes!(workspace) || changed)
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
    source_scan_running(workspace) && error("A source scan is already running")
    workspace.closed && error("The workspace is closed")

    cachedb = workspace.cache.db
    job = workspace.scan
    job.id += 1
    scan_id = job.id
    job.state = :discovering
    job.error = ""
    workspace.cache_state = :loading
    workspace.cache_error = ""
    workspace.cache.operation = rebuild ? :rebuild : :update
    reset_build_metrics!(workspace.metrics)
    reset_scan_profile!(workspace.project)
    # Detach the previous scan so this one streams into a fresh hierarchy: progressive results update
    # the displayed tree live while any still-finishing analysis from the previous scan keeps reading
    # its own, now-detached hierarchy — no shared mutable state, no race. Errors reset so stale ones
    # from the previous scan never linger.
    workspace.index.source = nothing
    empty!(workspace.index.analysis_errors)
    cancel_token = Base.Threads.Atomic{Bool}(false)
    job.cancel_token = cancel_token
    task = Base.Threads.@spawn begin
        Profiling.@profile_span workspace.profiler :source :scan Profiling.ProfileAttributes(
            source_id=source_id(workspace.source),
        ) begin
        try
            rebuild && clear_cache_index!(cachedb)
            cached = !rebuild && cache_built(cachedb) ? load_cache_index(cachedb) : nothing
            publish_cache_state!(
                workspace, scan_id, cached === nothing ? :missing : :ready, cached)

            write_meta_header!(cachedb)
            discovered = with_cancel(() -> cancel_token[]) do
                source_items(workspace.source)
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
            publish_source_changes!(
                workspace,
                scan_id,
                SourceChanges(upserts, removals; parameters_changed=false),
                current,
                status,
                isempty(upserts) && isempty(removals),
            )
            publish_scan_end!(workspace, scan_id)
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
    if profiler.state === :preparing
        workspace.profile_restart_pending = false
        profiler.state = :idle
        return nothing
    end
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
    parameter_keys = Set{Symbol}()
    for item in hierarchy.all_items
        haskey(items, item.id) && error("Duplicate item id generated during scan: $(item.id)")
        items[item.id] = item
    end
    for node in values(hierarchy.index)
        union!(parameter_keys, keys(node.parameters))
    end
    item_stats = Dict{String,Dict{Symbol,Any}}()
    for item in hierarchy.all_items
        previous = get(workspace.index.item_stats, item.id, nothing)
        if previous !== nothing
            item_stats[item.id] = previous
        elseif !isempty(item.stats)
            item_stats[item.id] = copy(item.stats)
        end
    end
    workspace.index = WorkspaceIndex(
        hierarchy,
        items,
        item_stats,
        sort!(collect(parameter_keys); by=String),
        workspace.index.source,
        workspace.index.analysis_errors,
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
    # parameters so cached and memory-only startup produce the same effective inputs.
    hierarchy = Hierarchy(
        index.source.hierarchy.all_items,
        workspace.source,
        index.source.hierarchy.skipped_count,
    )
    valid_item_stats = Dict{String,Dict{Symbol,Any}}()
    for record in index.source.hierarchy.all_items
        materialized = effective_record(hierarchy, record)
        state = get(
            index.result_states,
            CacheResultKey(ITEM_STATS_RESULT, record.id),
            nothing,
        )
        valid_stats = state !== nothing &&
            CacheResultStatus(state.status) === RESULT_READY &&
            state.input_fingerprint == result_input_fingerprint(materialized)
        valid_stats && haskey(index.item_stats, record.id) &&
            (valid_item_stats[record.id] = copy(index.item_stats[record.id]))
    end
    replace_item_index!(workspace, hierarchy)
    workspace.index.item_stats = valid_item_stats
    for (path, cached_node) in index.source.hierarchy.index
        node = get(hierarchy.index, path, nothing)
        node === nothing && continue
        key = collection_path_key(collect(path))
        state = get(index.result_states, CacheResultKey(COLLECTION_STATS_RESULT, key), nothing)
        state !== nothing &&
            CacheResultStatus(state.status) === RESULT_READY &&
            state.input_fingerprint == collection_input_fingerprint(workspace, path) ||
            continue
        merge!(node.stats, cached_node.stats)
    end
    workspace.index.source = SourceScan(
        index.source.source_id,
        index.source.source_label,
        index.source.source_item_fingerprints,
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
            kind === PROCESSING_RESULT ? PROCESS_ITEM :
            kind === ITEM_STATS_RESULT ? ITEM_STATS : COLLECTION_STATS,
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
    # Cached derived results whose effective inputs changed while the workspace was closed did not
    # seed above; without a fresh interpretation nothing else re-enqueues them.
    if workspace.background_processing
        missing_stats = lock(workspace.work.lock) do
            ItemRecord[
                record for record in workspace.index.hierarchy.all_items
                if !haskey(workspace.work.nodes, WorkKey(ITEM_STATS, record.id))
            ]
        end
        for record in missing_stats
            enqueue_processing!(workspace, record)
            enqueue_item_stats!(workspace, record)
        end
    end
    return nothing
end

"""Seed work graph completion nodes from a loaded cache index."""
function seed_cached_work_nodes!(
    workspace::Workspace,
    index::ProjectCacheIndex,
)::Nothing
    records = Dict(record.id => record for record in workspace.index.hierarchy.all_items)
    for (result_key, state) in index.result_states
        work_state = CacheResultStatus(state.status) === RESULT_READY ? :ready : :failed
        if result_key.kind === PROCESSING_RESULT
            record = get(records, result_key.entity, nothing)
            record === nothing && continue
            materialized = effective_record(workspace.index.hierarchy, record)
            state.input_fingerprint == result_input_fingerprint(materialized) || continue
            seed_work_node!(
                workspace,
                WorkKey(PROCESS_ITEM, result_key.entity),
                item_revision(workspace, record),
                work_state,
            )
        elseif result_key.kind === ITEM_STATS_RESULT
            record = get(records, result_key.entity, nothing)
            record === nothing && continue
            materialized = effective_record(workspace.index.hierarchy, record)
            state.input_fingerprint == result_input_fingerprint(materialized) || continue
            seed_work_node!(
                workspace,
                WorkKey(ITEM_STATS, result_key.entity),
                item_revision(workspace, record),
                work_state,
            )
        elseif result_key.kind === COLLECTION_STATS_RESULT
            path = collection_path_tuple(result_key.entity)
            haskey(workspace.index.hierarchy.index, path) || continue
            state.input_fingerprint == collection_input_fingerprint(workspace, path) || continue
            seed_work_node!(
                workspace,
                WorkKey(COLLECTION_STATS, result_key.entity),
                collection_revision(workspace, path),
                work_state,
            )
        end
    end
    return nothing
end

"""Mark invalidated collection-stat graph nodes stale without exposing cache keys to Workspace."""
function mark_collection_stats_missing!(
    workspace::Workspace,
    collection_keys::Vector{String},
)::Nothing
    for collection_key in unique(collection_keys)
        path = collection_path_tuple(collection_key)
        haskey(workspace.index.hierarchy.index, path) || continue
        mark_work_missing!(
            workspace,
            WorkKey(COLLECTION_STATS, collection_key),
            collection_revision(workspace, path),
        )
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
        notify(workspace.work.condition; all=true)
        waiters
    end
end

"""Queue item stats after processing succeeds."""
function enqueue_item_stats!(workspace::Workspace, record::ItemRecord)::Nothing
    enqueue_work!(
        workspace,
        WorkKey(ITEM_STATS, record.id),
        item_revision(workspace, record);
        priority=2,
        dependencies=WorkKey[WorkKey(PROCESS_ITEM, record.id)],
    )
    return nothing
end

"""Queue collection stats when all current member item-stat nodes are terminal."""
function enqueue_ready_collection_stats!(
    workspace::Workspace,
    collection_keys::Vector{String},
)::Nothing
    for collection_key in collection_keys
        path = collection_path_tuple(collection_key)
        node = get(workspace.index.hierarchy.index, path, nothing)
        node === nothing && continue
        isempty(node.items) && continue
        dependencies = WorkKey[WorkKey(ITEM_STATS, record.id) for record in node.items]
        enqueue_work!(
            workspace,
            WorkKey(COLLECTION_STATS, collection_key),
            collection_revision(workspace, path);
            priority=1,
            dependencies,
        )
    end
    return nothing
end

# In both publish functions, waiter wakeup must happen before any call that can throw, so blocked
# `request_processed_items` callers are never left waiting on a `take!` that will not resolve.
"""Apply one successful work completion to WorkspaceIndex and ProjectCache."""
function publish_work_success!(
    workspace::Workspace,
    node::WorkNode,
    result,
)::Bool
    key = node.key
    if key.kind === INTERPRET_SOURCE
        records = result.records::Vector{ItemRecord}
        old_records, invalidated = publish_source_item_records!(
            workspace, key.entity, records)
        delete_collection_stats!(workspace.cache.db, invalidated)
        mark_collection_stats_missing!(workspace, invalidated)
        store_interpreted_records!(
            workspace.cache.db,
            records,
        )
        for record in records
            mark_work_missing!(
                workspace,
                WorkKey(ITEM_STATS, record.id),
                item_revision(workspace, record);
                dependencies=WorkKey[WorkKey(PROCESS_ITEM, record.id)],
            )
            workspace.background_processing && enqueue_processing!(workspace, record)
        end
        return true
    elseif key.kind === PROCESS_ITEM
        record = workspace.index.items[key.entity]
        materialized_record = effective_record(workspace.index.hierarchy, record)
        store_processed!(workspace.cache.db, materialized_record, result.item)
        waiters = finish_work_node!(workspace, node, :ready)
        loaded = DataItem(materialized_record, item_data(result.item))
        for waiter in waiters
            put!(waiter, ProcessingResult(loaded, nothing))
        end
        enqueue_item_stats!(workspace, record)
        return false
    elseif key.kind === ITEM_STATS
        record = workspace.index.items[key.entity]
        stats = metadata_dict(result)
        workspace.index.item_stats[key.entity] = stats
        store_item_stats!(
            workspace.cache.db,
            effective_record(workspace.index.hierarchy, record),
            stats,
        )
        enqueue_ready_collection_stats!(
            workspace,
            affected_collection_keys([record]),
        )
        return true
    else
        path = collection_path_tuple(key.entity)
        hierarchy_node = get(workspace.index.hierarchy.index, path, nothing)
        hierarchy_node === nothing && return false
        empty!(hierarchy_node.stats)
        merge!(hierarchy_node.stats, metadata_dict(result))
        store_collection_stats!(
            workspace.cache.db,
            key.entity,
            collection_input_fingerprint(workspace, path),
            hierarchy_node.stats,
        )
        return true
    end
end

"""Apply one failed work completion to WorkspaceIndex and ProjectCache."""
function publish_work_failure!(
    workspace::Workspace,
    node::WorkNode,
    failure::CapturedException,
)::Bool
    key = node.key
    message = sprint(showerror, failure.ex)
    if key.kind === INTERPRET_SOURCE
        old_records, invalidated = remove_source_item_output!(workspace, key.entity)
        delete_source_item!(workspace.cache.db, key.entity, old_records)
        delete_collection_stats!(workspace.cache.db, invalidated)
        mark_collection_stats_missing!(workspace, invalidated)
        workspace.index.analysis_errors[key.entity] = "interpret_source_item: " * message
        finish_work_node!(workspace, node, :failed)
        return true
    end
    source_item_id_value = ""
    record = key.kind in (PROCESS_ITEM, ITEM_STATS) ?
        get(workspace.index.items, key.entity, nothing) : nothing
    record === nothing || (source_item_id_value = record.source_item_id)
    input_fingerprint = if record !== nothing
        result_input_fingerprint(effective_record(workspace.index.hierarchy, record))
    elseif key.kind === COLLECTION_STATS
        collection_input_fingerprint(workspace, collection_path_tuple(key.entity))
    else
        ""
    end
    store_result_failure!(
        workspace.cache.db,
        CacheResultKey(
            key.kind === PROCESS_ITEM ? PROCESSING_RESULT :
            key.kind === ITEM_STATS ? ITEM_STATS_RESULT : COLLECTION_STATS_RESULT,
            key.entity,
        ),
        source_item_id_value,
        input_fingerprint,
        message,
    )
    workspace.index.analysis_errors[key.entity] = message
    waiters = finish_work_node!(workspace, node, :failed)
    for waiter in waiters
        put!(waiter, ProcessingResult(nothing, failure))
    end
    if key.kind === ITEM_STATS && record !== nothing
        enqueue_ready_collection_stats!(workspace, affected_collection_keys([record]))
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
            key.kind === PROCESS_ITEM || finish_work_node!(workspace, node, :ready)
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
            replace_item_index!(workspace, Hierarchy(source_id(workspace.source), false))
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
    cache_hit::Bool,
)::Nothing
    lock(workspace.publish_lock) do
        scan_id == workspace.scan.id || return
        apply_source_changes!(workspace, changes, fingerprints, status)
        workspace.scan.state = cache_hit ? :unchanged : :done
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

"""Publish the end of one scan job: finished normally, canceled, or failed."""
function publish_scan_end!(
    workspace::Workspace,
    scan_id::Int;
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
