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
)::Workspace
    workspace = Workspace(
        project,
        source;
        profile_internal,
        profile_cpu,
        profile_output,
        crash_trace,
        rebuild,
    )
    scan_source!(workspace; rebuild)
    return workspace
end

"""
Cancel all work owned by a workspace and wait for it to stop.
"""
function close_workspace!(workspace::Workspace)::Nothing
    workspace.closed && return nothing
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
            workspace.closed = true
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
    workspace.scan.state in (:counting, :discovering, :scanning, :canceling)

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

"""Atomically publish one source item's replacement records after validating duplicate IDs."""
function replace_source_item_output!(
    workspace::Workspace,
    source_item_id_value::String,
    replacement::SourceItemInterpretation,
)::Tuple{Vector{ItemRecord},Vector{String}}
    old_records = source_item_records(workspace.index, source_item_id_value)
    old_ids = Set(record.id for record in old_records)
    replacement_ids = Set{String}()
    for record in replacement.records
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
    new_keys = affected_collection_keys(replacement.records)
    old_records, _ = remove_source_item_output!(workspace, source_item_id_value)
    for record in replacement.records
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
    status::ProjectCacheStatus,
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
    refresh_workspace_source!(workspace, fingerprints)
    workspace.cache.status = status
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
    job.progress = WorkspaceProgress()
    job.error = ""
    workspace.cache_state = :loading
    workspace.cache_error = ""
    workspace.cache.operation = rebuild ? :rebuild : :update
    reset_build_metrics!(workspace.metrics)
    # Detach the previous scan so this one streams into a fresh hierarchy: progressive results update
    # the displayed tree live while any still-finishing analysis from the previous scan keeps reading
    # its own, now-detached hierarchy — no shared mutable state, no race. Errors reset so stale ones
    # from the previous scan never linger.
    workspace.index.source = nothing
    empty!(workspace.index.analysis_errors)
    events = Channel{NamedTuple}(Inf)
    cancel_token = Base.Threads.Atomic{Bool}(false)
    job.events = events
    job.cancel_token = cancel_token
    task = Base.Threads.@spawn begin
        Profiling.@profile_span workspace.profiler :source :scan Profiling.ProfileAttributes(
            source_id=source_id(workspace.source),
        ) begin
        try
            rebuild && clear_cache_index!(cachedb)
            cached = !rebuild && cache_built(cachedb) ? load_cache_index(cachedb) : nothing
            put!(events, (
                kind=:cache_state,
                job_id=scan_id,
                state=cached === nothing ? :missing : :ready,
                index=cached,
            ))

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
                if !haskey(previous, id_value) || !isequal(previous[id_value], current[id_value])
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
            put!(events, (
                kind=:source_changes,
                job_id=scan_id,
                changes=SourceChanges(upserts, removals),
                fingerprints=current,
                status,
                cache_hit=isempty(upserts) && isempty(removals),
            ))
            put!(events, (
                kind=:source_done,
                job_id=scan_id,
            ))
        catch error
            if is_job_cancelled(error)
                put!(events, (kind=:canceled, job_id=scan_id))
            else
                put!(events, (
                    kind=:error,
                    job_id=scan_id,
                    error,
                    backtrace=catch_backtrace(),
                ))
            end
        finally
            close(events)
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
        scan_done=workspace.scan.progress.processed_source_items,
        scan_total=workspace.scan.progress.total_source_items,
        processing_done=completed,
        processing_total=total,
        queue_depth=depth,
    )
end

"""Begin recording immediately when idle, or prepare one clean rebuild after active work stops."""
function start_internal_profile!(workspace::Workspace)::Nothing
    profiler = workspace.profiler
    Profiling.validate_start(profiler)
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

"""Add newly discovered items to the provisional workspace index."""
function append_items!(
    workspace::Workspace,
    items::Vector{ItemRecord},
)::Bool
    isempty(items) && return false
    hierarchy = workspace.index.hierarchy
    item_index = workspace.index.items
    changed = false
    for item in items
        haskey(item_index, item.id) && continue
        insert_item!(hierarchy, item)
        item_index[item.id] = item
        changed = true
    end
    changed || return false
    apply_collection_parameters!(hierarchy, workspace.source)
    parameter_keys = Set{Symbol}()
    for node in values(hierarchy.index)
        union!(parameter_keys, keys(node.parameters))
    end
    workspace.index.collection_parameter_keys = sort!(collect(parameter_keys); by=String)
    return true
end

"""Surface the cached index as an instant, independent first view, before the fresh scan refines it."""
function apply_cache_index!(
    workspace::Workspace,
    index::ProjectCacheIndex,
)::Nothing
    index.identity.source_id == source_id(workspace.source) ||
        error("Loaded cache belongs to $(index.identity.source_id), not $(source_id(workspace.source))")
    workspace.cache.identity = index.identity
    workspace.cache.index = index

    # The cached hierarchy is its own freshly-loaded snapshot (distinct record objects), so the live
    # scan can refine it without touching any hierarchy a previous scan's analysis still reads.
    replace_item_index!(workspace, index.source.hierarchy)
    workspace.index.source = index.source

    cached_items = length(index.source.source_item_fingerprints)
    workspace.cache.status = ProjectCacheStatus(
        cached_items, cached_items, 0, 0, 0, 0, length(index.analysis_errors))
    workspace.index.analysis_errors = copy(index.analysis_errors)
    for state in values(index.result_states)
        CacheResultStatus(state.status) === RESULT_FAILED || continue
        workspace.index.analysis_errors[state.entity] =
            something(state.message, "cached work failed")
    end
    seed_cached_work_nodes!(workspace, index)
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
            seed_work_node!(
                workspace,
                WorkKey(PROCESS_ITEM, result_key.entity),
                item_revision(record),
                work_state,
            )
        elseif result_key.kind === ITEM_STATS_RESULT
            record = get(records, result_key.entity, nothing)
            record === nothing && continue
            seed_work_node!(
                workspace,
                WorkKey(ITEM_STATS, result_key.entity),
                item_revision(record),
                work_state,
            )
        elseif result_key.kind === COLLECTION_STATS_RESULT
            path = collection_path_tuple(result_key.entity)
            haskey(workspace.index.hierarchy.index, path) || continue
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
        item_revision(record);
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

"""Apply one successful work completion to WorkspaceIndex and ProjectCache."""
function publish_work_success!(
    workspace::Workspace,
    node::WorkNode,
    result,
)::Bool
    key = node.key
    if key.kind === INTERPRET_SOURCE
        interpretation = result::SourceItemInterpretation
        old_records, invalidated = replace_source_item_output!(
            workspace, key.entity, interpretation)
        delete_source_item!(workspace.cache.db, key.entity, old_records)
        delete_collection_stats!(workspace.cache.db, invalidated)
        mark_collection_stats_missing!(workspace, invalidated)
        store_interpreted!(
            workspace.cache.db,
            interpretation.records,
            interpretation.interpreted_items,
        )
        for record in interpretation.records
            mark_work_missing!(
                workspace,
                WorkKey(ITEM_STATS, record.id),
                item_revision(record);
                dependencies=WorkKey[WorkKey(PROCESS_ITEM, record.id)],
            )
            enqueue_processing!(workspace, record)
        end
        return true
    elseif key.kind === PROCESS_ITEM
        record = workspace.index.items[key.entity]
        store_processed!(workspace.cache.db, record, result.item)
        waiters = finish_work_node!(workspace, node, :ready)
        loaded = DataItem(record, item_data(result.item))
        for waiter in waiters
            put!(waiter, ProcessingResult(loaded, nothing))
        end
        enqueue_item_stats!(workspace, record)
        return false
    elseif key.kind === ITEM_STATS
        record = workspace.index.items[key.entity]
        stats = metadata_dict(result)
        workspace.index.item_stats[key.entity] = stats
        workspace.index.items[key.entity] = ItemRecord(record; stats)
        rebuild_workspace_hierarchy!(workspace)
        store_item_stats!(workspace.cache.db, workspace.index.items[key.entity], stats)
        enqueue_ready_collection_stats!(
            workspace,
            affected_collection_keys([workspace.index.items[key.entity]]),
        )
        return true
    else
        path = collection_path_tuple(key.entity)
        hierarchy_node = get(workspace.index.hierarchy.index, path, nothing)
        hierarchy_node === nothing && return false
        empty!(hierarchy_node.stats)
        merge!(hierarchy_node.stats, metadata_dict(result))
        store_collection_stats!(workspace.cache.db, key.entity, hierarchy_node.stats)
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
    store_result_failure!(
        workspace.cache.db,
        CacheResultKey(
            key.kind === PROCESS_ITEM ? PROCESSING_RESULT :
            key.kind === ITEM_STATS ? ITEM_STATS_RESULT : COLLECTION_STATS_RESULT,
            key.entity,
        ),
        source_item_id_value,
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

"""Apply graph completion events; stale revisions are ignored."""
function poll_work_events!(workspace::Workspace)::Bool
    changed = false
    events = workspace.work.events
    while isready(events)
        event = take!(events)
        node = current_completion_node(workspace, event.key, event.revision)
        node === nothing && continue
        if event.result isa CapturedException
            changed |= publish_work_failure!(workspace, node, event.result)
        else
            changed |= publish_work_success!(workspace, node, event.result)
            node.key.kind === PROCESS_ITEM || finish_work_node!(workspace, node, :ready)
        end
    end
    return changed
end

"""Apply all completed scan, cache, and work events."""
function poll_workspace!(workspace::Workspace)::Bool
    index_changed = false

    scan_events = workspace.scan.events
    if scan_events !== nothing
        while isready(scan_events)
            event = take!(scan_events)
            event.job_id == workspace.scan.id || continue
            if event.kind == :progress
                workspace.scan.progress = WorkspaceProgress(event.progress)
                workspace.scan.state = event.progress.phase
            elseif event.kind == :cache_state
                if event.index === nothing
                    replace_item_index!(workspace, Hierarchy(source_id(workspace.source), false))
                    empty!(workspace.index.analysis_errors)
                else
                    apply_cache_index!(workspace, event.index)
                end
                index_changed = true
                workspace.cache_state = event.state
                event.state == :missing && workspace.cache.operation != :rebuild &&
                    (workspace.cache.operation = :build)
            elseif event.kind == :source_changes
                index_changed |= apply_source_changes!(
                    workspace, event.changes, event.fingerprints, event.status)
                workspace.scan.state = event.cache_hit ? :unchanged : :done
                workspace.cache_state = :ready
            elseif event.kind == :source_done
                workspace.scan.events = nothing
                workspace.scan.cancel_token = nothing
            elseif event.kind == :canceled
                workspace.scan.state = :canceled
                workspace.scan.events = nothing
                workspace.scan.cancel_token = nothing
            elseif event.kind == :error
                workspace.scan.state = :error
                workspace.scan.error =
                    "Source scan failed. See the console for full details."
                workspace.cache_state == :writing &&
                    (workspace.cache_state = :error)
                @error(
                    "Source scan job failed",
                    source=source_id(workspace.source),
                    current_source_item=workspace.scan.progress.current_source_item,
                    exception=(event.error, event.backtrace),
                )
                workspace.scan.events = nothing
                workspace.scan.cancel_token = nothing
            end
        end
    end

    index_changed |= poll_work_events!(workspace)

    # Refresh the watcher snapshot while work runs (progress advances without an index change) and for
    # the one frame after it stops; idle frames skip the rebuild entirely.
    busy = source_scan_running(workspace) || processing_work_running(workspace) ||
           analysis_work_running(workspace) ||
           cache_work_running(workspace)
    if workspace.profile_restart_pending && !busy
        reset_work_graph!(workspace)
        workspace.profile_restart_pending = false
        workspace.profiler.state = :idle
        Profiling.start!(
            workspace.profiler,
            () -> profile_counter_snapshot(workspace),
        )
        scan_source!(workspace; rebuild=true)
        busy = true
    end
    (index_changed || busy || workspace.status.busy) &&
        (workspace.status = workspace_status(workspace))
    return index_changed
end
