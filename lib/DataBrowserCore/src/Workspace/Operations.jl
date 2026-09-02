"""
Start cache loading and source scanning for one new workspace.
"""
function open_workspace(
    project::AbstractProject,
    source::AbstractDataSource;
    rebuild::Bool=false,
    cache::Bool=true,
    background_processing::Bool=false,
)::Workspace
    opened_source = open_source(source)
    workspace = try
        Workspace(
            project,
            opened_source;
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
            event -> publish_source_event!(workspace, event);
            cancel_token=get_token(workspace.cancel_source),
        )
        return workspace
    catch
        close_workspace!(workspace)
        rethrow()
    end
end

"""
Rebuild `workspace` in place. The project cannot change. Omitted `source`, `cache`, and
`background_processing` stay as they are. Omitted `source` is a `copy` of the current source so the
closed source's watchers are not reused. `rebuild` is a one-shot cache action, not a stored setting.

A failed rebuild leaves the workspace closed.
"""
function modify_workspace!(
    workspace::Workspace;
    source::Union{Nothing,AbstractDataSource}=nothing,
    rebuild::Bool=false,
    cache::Union{Nothing,Bool}=nothing,
    background_processing::Union{Nothing,Bool}=nothing,
)::Nothing
    lock(workspace.lifecycle_lock) do
        workspace.closed && error("Cannot modify a closed workspace")
        next_source = source === nothing ? copy(workspace.source) : source
        next_cache = cache === nothing ? workspace.disk_cache : cache
        next_background =
            background_processing === nothing ? workspace.background_processing :
            background_processing
        opened_source = open_source(next_source)
        _stop_runtime!(workspace)
        try
            _install_runtime!(
                workspace,
                opened_source;
                rebuild,
                cache=next_cache,
                background_processing=next_background,
            )
            scan_source!(workspace; rebuild)
            watch_source(
                opened_source,
                event -> publish_source_event!(workspace, event);
                cancel_token=get_token(workspace.cancel_source),
            )
        catch
            workspace.closed = true
            try
                _stop_runtime!(workspace)
            catch
            end
            rethrow()
        end
    end
    return nothing
end

"""Stop workers, scan, cache, and source without marking the workspace closed."""
function _stop_runtime!(workspace::Workspace)::Nothing
    cancel(workspace.cancel_source)
    cancel_scan!(workspace)
    for task in workspace.background_tasks
        istaskdone(task) || wait(task)
    end
    empty!(workspace.background_tasks)
    stop_work_workers!(workspace)
    try
        close_cache_db!(workspace.cache.db)
    finally
        close_source!(workspace.source)
    end
    return nothing
end

"""Replace runtime fields after `_stop_runtime!` and start cache workers."""
function _install_runtime!(
    workspace::Workspace,
    source::AbstractDataSource;
    rebuild::Bool,
    cache::Bool,
    background_processing::Bool,
)::Nothing
    metrics = BuildMetrics()
    cache_db, disk_error, identity =
        _open_workspace_cache(workspace.project, source, metrics; rebuild, cache)
    workspace.source = source
    workspace.index = WorkspaceIndex(
        CollectionIndex(source_id(source)),
        Dict{String,ItemRecord}(),
        Dict{String,Dict{Symbol,Any}}(),
        Symbol[],
        nothing,
        Dict{Union{String,Int64},String}(),
        Dict{Int64,Vector{String}}(),
    )
    workspace.cache = WorkspaceCache(identity, cache_db, disk_error, nothing, :load)
    workspace.cache_state = :idle
    workspace.cache_error = ""
    workspace.source_error = ""
    workspace.work = WorkDependencyGraph()
    workspace.background_processing = background_processing
    workspace.disk_cache = cache
    workspace.metrics = metrics
    workspace.status = WorkspaceStatus()
    workspace.status_dirty[] = true
    workspace.cancel_source = CancellationTokenSource()
    start_cache!(cache_db)
    start_work_workers!(workspace)
    return nothing
end

"""
Cancel all work owned by a workspace and wait for it to stop.
"""
function close_workspace!(workspace::Workspace)::Nothing
    lock(workspace.lifecycle_lock) do
        workspace.closed && return nothing
        workspace.closed = true
        _stop_runtime!(workspace)
        return nothing
    end
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
        scan = workspace.scan
        scan.state === :discovering || return
        token_source = scan.cancel_token
        token_source !== nothing && cancel(token_source)
        scan.state = :canceling
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
function source_item_records(index::WorkspaceIndex, source_item_key_value::Int64)::Vector{ItemRecord}
    ids = get(index.items_by_source, source_item_key_value, nothing)
    ids === nothing && return ItemRecord[]
    return ItemRecord[index.items[id] for id in ids if haskey(index.items, id)]
end

"""Resolve one source-item surrogate to the stable public source-item id."""
source_item_id(workspace::Workspace, key::Integer)::String =
    DataBrowserCache.source_item_id(workspace.cache.db, key)

"""Resolve one public source-item id to its surrogate, or `nothing` when unknown."""
source_item_key(workspace::Workspace, id_value::AbstractString)::Union{Nothing,Int64} =
    DataBrowserCache.source_item_key(workspace.cache.db, id_value)

"""Return collection keys touched by records and their ancestors."""
function affected_collection_keys(
    collections::CollectionIndex,
    records::Vector{ItemRecord},
)::Vector{Int64}
    keys = Int64[]
    for record in records
        for key in collection_path_keys(collections, record.collection_key)
            push!(keys, key)
        end
    end
    return unique(keys)
end

"""Refresh the collection-metadata key list from the current records."""
function refresh_collection_metadata_keys!(index::WorkspaceIndex)::Nothing
    metadata_keys = Set{Symbol}()
    for collection_record in values(index.collections.records)
        union!(metadata_keys, keys(collection_record.own_metadata))
        union!(metadata_keys, keys(collection_record.analysis))
    end
    index.collection_metadata_keys = sort!(collect(metadata_keys); by=String)
    return nothing
end

"""
    include_collection_metadata_keys!(index, collection_keys) -> Nothing

Add metadata keys from newly registered collection records without rescanning the complete
collection index.
"""
function include_collection_metadata_keys!(index::WorkspaceIndex, collection_keys)::Nothing
    changed = false
    for collection_key in collection_keys
        collection_record = get(index.collections.records, collection_key, nothing)
        collection_record === nothing && continue
        for metadata_key in keys(collection_record.own_metadata)
            metadata_key in index.collection_metadata_keys && continue
            push!(index.collection_metadata_keys, metadata_key)
            changed = true
        end
        for metadata_key in keys(collection_record.analysis)
            metadata_key in index.collection_metadata_keys && continue
            push!(index.collection_metadata_keys, metadata_key)
            changed = true
        end
    end
    changed && sort!(index.collection_metadata_keys; by=String)
    return nothing
end

"""Refresh derived collection metadata after indexed records change."""
function rebuild_workspace_hierarchy!(workspace::Workspace)::Nothing
    refresh_collection_metadata_keys!(workspace.index)
    return nothing
end

"""Rebuild effective collection metadata from the collection values currently in the index."""
function reconcile_source_metadata_cache!(
    workspace::Workspace;
    collections::Union{Nothing,Vector{Int64}}=nothing,
    refresh_hierarchy::Bool=false,
)::Vector{ItemRecord}
    refresh_hierarchy || return ItemRecord[]
    old_collections = @timed_dbg "recon_copy" copy(workspace.index.collections)
    current_keys = collect(keys(workspace.index.collections.records))
    @timed_dbg "recon_resolve" for key in current_keys
        haskey(workspace.index.collections.records, key) || continue
        path = collection_value_path(workspace.index.collections, key)
        isempty(path) && continue
        inputs = collection_inputs(annotate_collection_path(workspace.source, path))
        resolve_collection_path!(
            workspace.index.collections,
            inputs;
            update_existing=true,
        )
    end
    @timed_dbg "recon_keys" refresh_collection_metadata_keys!(workspace.index)
    records = collect(values(workspace.index.items))
    changed = @timed_dbg "recon_diff" ItemRecord[
        record for record in records
        if effective_metadata(old_collections, record) !=
           effective_metadata(workspace.index.collections, record)
    ]
    isempty(records) || @timed_dbg "recon_store" store_collection_index!(
        workspace.cache.db, workspace.index.collections, records)
    return changed
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
    collection_keys::Vector{Int64}=
        affected_collection_keys(workspace.index.collections, records),
    priority::Int=2,
)::Nothing
    background = workspace.background_processing
    for collection_key in unique(collection_keys)
        haskey(workspace.index.collections.records, collection_key) || continue
        clear_collection_analysis!(workspace.index.collections, collection_key)
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
        copy(workspace.index.collections),
        sort!(collect(values(workspace.index.items)); by=record ->
            (record.source_item_key, record.id)),
        ItemFailure[],
    )
    return nothing
end

"""Remove all published output owned by one source item and invalidate affected collections."""
function remove_source_item_output!(
    workspace::Workspace,
    source_item_key_value::Int64,
)::Tuple{Vector{ItemRecord},Vector{Int64}}
    old_records = source_item_records(workspace.index, source_item_key_value)
    invalidated = affected_collection_keys(workspace.index.collections, old_records)
    items = workspace.index.items
    collections = workspace.index.collections
    for record in old_records
        remove_item!(collections, record.id, record.collection_key)
        delete!(items, record.id)
        delete!(workspace.index.item_metadata, record.id)
        delete!(workspace.index.analysis_errors, record.id)
    end
    # Source-level analysis errors stay keyed by the public source-item id.
    delete!(workspace.index.analysis_errors, source_item_id(workspace, source_item_key_value))
    delete!(workspace.index.items_by_source, source_item_key_value)
    analysis_cleared = false
    for collection_key in invalidated
        haskey(collections.records, collection_key) || continue
        analysis_cleared |= !isempty(collections.records[collection_key].analysis)
        clear_collection_analysis!(collections, collection_key)
    end
    (analysis_cleared || any(key -> !haskey(collections.records, key), invalidated)) &&
        refresh_collection_metadata_keys!(workspace.index)
    filter!(id -> haskey(items, id), workspace.selection.item_ids)
    return old_records, invalidated
end

"""Delete persisted collection rows no longer present in the published collection index."""
function delete_pruned_collection_records!(
    workspace::Workspace,
    collection_keys::Vector{Int64},
)::Nothing
    pruned = Int64[
        key for key in unique(collection_keys)
        if !haskey(workspace.index.collections.records, key)
    ]
    delete_collection_records!(workspace.cache.db, pruned)
    return nothing
end

"""Publish one source item's current records after validating duplicate IDs."""
function publish_source_item_records!(
    workspace::Workspace,
    source_item_key_value::Int64,
    records::Vector{ItemRecord},
    collection_paths::Vector{Vector{CollectionInput}},
)::Tuple{Vector{ItemRecord},Vector{Int64}}
    source_ref = source_item_id(workspace, source_item_key_value)
    length(records) == length(collection_paths) || error(
        "Cannot publish source item '$source_ref': records and collection paths differ in length",
    )
    old_records = source_item_records(workspace.index, source_item_key_value)
    old_ids = Set(record.id for record in old_records)
    replacement_ids = Set{String}()
    for record in records
        record.source_item_key == source_item_key_value || error(
            "Cannot publish source item '$source_ref': replacement contains " *
            "record '$(record.id)' from source item key $(record.source_item_key)",
        )
        record.id in replacement_ids && error(
            "Cannot publish source item '$source_ref': duplicate item id '$(record.id)'",
        )
        push!(replacement_ids, record.id)
        haskey(workspace.index.items, record.id) && !(record.id in old_ids) && error(
            "Cannot publish source item '$source_ref': replacement item id " *
            "'$(record.id)' already belongs to another source item",
        )
    end

    collections = workspace.index.collections
    next_collection_key = collections.next_key
    collection_keys = resolve_collection_paths!(collections, collection_paths)
    resolved = ItemRecord[
        ItemRecord(record; collection_key)
        for (record, collection_key) in zip(records, collection_keys)
    ]
    old_keys = affected_collection_keys(workspace.index.collections, old_records)
    items = workspace.index.items
    if isempty(old_records)
        for record in resolved
            append_item!(collections, record.id, record.collection_key)
        end
    else
        old_memberships = Set((record.id, record.collection_key) for record in old_records)
        retained = Set((record.id, record.collection_key) for record in resolved)
        for record in resolved
            (record.id, record.collection_key) in old_memberships ||
                append_item!(collections, record.id, record.collection_key)
        end
        for record in old_records
            (record.id, record.collection_key) in retained ||
                remove_item!(collections, record.id, record.collection_key)
            delete!(items, record.id)
            delete!(workspace.index.item_metadata, record.id)
            delete!(workspace.index.analysis_errors, record.id)
        end
    end
    delete!(workspace.index.analysis_errors, source_ref)
    for record in resolved
        items[record.id] = record
        delete!(workspace.index.item_metadata, record.id)
    end
    if isempty(resolved)
        delete!(workspace.index.items_by_source, source_item_key_value)
    else
        workspace.index.items_by_source[source_item_key_value] =
            String[record.id for record in resolved]
    end
    filter!(id -> haskey(items, id), workspace.selection.item_ids)
    new_keys = affected_collection_keys(collections, resolved)
    invalidated = unique([old_keys; new_keys])
    analysis_cleared = false
    for collection_key in invalidated
        haskey(collections.records, collection_key) || continue
        analysis_cleared |= !isempty(collections.records[collection_key].analysis)
        clear_collection_analysis!(collections, collection_key)
    end
    if analysis_cleared || any(key -> !haskey(collections.records, key), old_keys)
        refresh_collection_metadata_keys!(workspace.index)
    elseif next_collection_key < collections.next_key
        include_collection_metadata_keys!(
            workspace.index,
            next_collection_key:(collections.next_key - 1),
        )
    end
    return resolved, invalidated
end

"""True when the stage ledger already tracks this source item (cached or failed interpret)."""
function _cache_knows_source_item(cachedb::AbstractCacheDB, source_key::Int64)::Bool
    ledger = cachedb.stage_ledger
    return lock(ledger.lock) do
        source_key in ledger.source_items || ("", source_key) in ledger.failures
    end
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
        removed_key = source_item_key(workspace, removed)
        removed_key === nothing && continue
        old_records, invalidated = remove_source_item_output!(workspace, removed_key)
        delete_source_item!(workspace.cache.db, removed_key, old_records)
        delete_collection_metadata!(workspace.cache.db, invalidated)
        delete_pruned_collection_records!(workspace, invalidated)
        enqueue_dependent_subtree!(workspace, old_records; collection_keys=invalidated)
        lock(workspace.work.lock) do
            delete!(workspace.work.source_items, removed_key)
        end
        changed = true
    end
    isempty(changes.upserts) || lock(workspace.work.lock) do
        for source_item in changes.upserts
            source_item_key_value =
                source_item_key!(workspace.cache.db, id(source_item); mint=true)
            old_records = source_item_records(workspace.index, source_item_key_value)
            # Fresh upserts have no published output and no cache row — skip the no-op delete so
            # discovery does not flood the write buffer with BUFFER_DELETE for every new file.
            if !isempty(old_records) ||
               _cache_knows_source_item(workspace.cache.db, source_item_key_value)
                delete_source_item!(workspace.cache.db, source_item_key_value, old_records)
            end
            key = WorkKey(SOURCE_INTERPRET, source_item_key_value)
            workspace.work.source_items[source_item_key_value] = source_item
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
        stale = @timed_dbg "reconcile_metadata" reconcile_source_metadata_cache!(workspace; refresh_hierarchy=true)
        isempty(stale) || @timed_dbg "invalidate_work" invalidate_records_work!(workspace, stale)
        changed = changed || !isempty(stale)
    end
    # Upsert-only batches only enqueue work; the published index is unchanged until interpretation
    # lands or removals/metadata reconcile run. Rebuild SourceScan then, not on every discover batch.
    if !isempty(changes.removals) || changes.metadata_changed
        @timed_dbg "refresh_source" refresh_workspace_source!(workspace)
    end
    status === nothing || (workspace.cache.status = status)
    return changed
end

"""
Start one source scan that progressively populates the workspace index and cache.

The scan first surfaces any already-cached index for an instant first view, then walks the source
while only fingerprinting new or stale items into a queue — so `sources_found` advances at walk
speed. After the walk, a publisher task drains that queue in small batches so `SOURCE_INTERPRET`
starts immediately and overlaps the remaining enqueue (without serializing the walk behind
`publish_lock` / ingest or contending with interpreters on the filesystem). Removals and
source-metadata reconcile run in a final batch after the publisher drains. Interpreted results
stream into the workspace and a bounded cache-write batch that flushes by source-item count, cached
row count, or age. `rebuild=true` first wipes the cache and ignores the prior scan, forcing a full
re-interpretation.
"""
function scan_source!(
    workspace::Workspace;
    rebuild::Bool=false,
)::Nothing
    cachedb = workspace.cache.db
    scan = workspace.scan
    # Scan setup runs under the publish lock so cancel_scan!'s running check can never interleave
    # with the transition into :discovering.
    scan_epoch, cancel_source = lock(workspace.publish_lock) do
        source_scan_running(workspace) && error("A source scan is already running")
        workspace.closed && error("The workspace is closed")
        scan.epoch += 1
        scan.state = :discovering
        scan.error = ""
        scan.discovered[] = 0
        workspace.cache_state = :loading
        workspace.cache_error = ""
        workspace.cache.operation = rebuild ? :rebuild : :update
        reset_build_metrics!(workspace.metrics)
        # Detach the previous scan so this one streams into a fresh hierarchy: progressive results
        # update the displayed tree live while any still-finishing analysis from the previous scan
        # keeps reading its own, now-detached hierarchy — no shared mutable state, no race. Errors
        # reset so stale ones from the previous scan never linger.
        workspace.index.source = nothing
        empty!(workspace.index.analysis_errors)
        cancel_source = CancellationTokenSource(get_token(workspace.cancel_source))
        scan.cancel_token = cancel_source
        workspace.status_dirty[] = true
        (scan.epoch, cancel_source)
    end
    task = Base.Threads.@spawn begin
        scan_token = get_token(cancel_source)
        @timed_dbg "scan_source" begin
            try
                rebuild && clear_cache_index!(cachedb)
                cached = if !rebuild && cache_built(cachedb)
                    try
                        @timed_dbg load_cache_index(cachedb)
                    catch error
                        error isa ProjectCacheDataError || rethrow()
                        @warn(
                            "Generated collection cache is incompatible; rebuilding it",
                            cache=workspace.cache.identity.cache_path,
                            error,
                        )
                        clear_cache_index!(cachedb)
                        nothing
                    end
                else
                    nothing
                end
                publish_cache_state!(
                    workspace, scan_epoch, cached === nothing ? :missing : :ready, cached)

                write_meta_header!(cachedb)
                # Fingerprints must be loaded before discovery so upserts can stream during the walk.
                # The cache's `source_items` table is the sole home of previous-session fingerprints;
                # a memory-only cache holds nothing, so every discovered item re-interprets.
                previous = _load_source_item_fingerprints(cachedb)
                current = Dict{String,Any}()
                seen = Set{String}()
                pending_upserts = Ref{Vector{AbstractDataSourceItem}}(AbstractDataSourceItem[])
                any_upserts = Ref(false)
                streamed = Ref(false)
                # Walk only fingerprints into an unbounded channel — do not publish yet. Starting
                # SOURCE_INTERPRET during the walk contends on the same filesystem (and publish_lock)
                # and made discovery unusably slow on large trees. After the walk, a publisher task
                # drains batches so interpretation overlaps the remaining enqueue, not the walk.
                upsert_batch = 32
                upsert_channel = Channel{Vector{AbstractDataSourceItem}}(Inf)
                flush_pending_upserts! = () -> begin
                    isempty(pending_upserts[]) && return
                    batch = pending_upserts[]
                    pending_upserts[] = AbstractDataSourceItem[]
                    put!(upsert_channel, batch)
                end
                consider_source_item! = item -> begin
                    is_cancellation_requested(scan_token) &&
                        throw(OperationCanceledException(scan_token))
                    id_value = id(item)
                    push!(seen, id_value)
                    current_fingerprint = fingerprint(item)
                    current[id_value] = current_fingerprint
                    # A `nothing` fingerprint means the source cannot prove the item is unchanged,
                    # so it is always re-read.
                    if !haskey(previous, id_value) || current_fingerprint === nothing ||
                       !isequal(previous[id_value], current_fingerprint)
                        push!(pending_upserts[], item)
                        length(pending_upserts[]) >= upsert_batch && flush_pending_upserts!()
                    end
                end
                discovered = @timed_dbg "source_walk" source_items(
                    workspace.source;
                    cancel_token=scan_token,
                    on_progress=count -> begin
                        scan.discovered[] = count
                        workspace.status_dirty[] = true
                    end,
                    on_item=item -> begin
                        streamed[] = true
                        consider_source_item!(item)
                    end,
                )
                is_cancellation_requested(scan_token) &&
                    throw(OperationCanceledException(scan_token))
                # Sources that ignore `on_item` still queue batches after discovery returns.
                if !streamed[]
                    for item in discovered
                        consider_source_item!(item)
                    end
                end
                flush_pending_upserts!()
                publisher = Threads.@spawn begin
                    for batch in upsert_channel
                        any_upserts[] = true
                        publish_source_changes!(
                            workspace,
                            scan_epoch,
                            SourceChanges(batch, String[]; metadata_changed=false),
                        )
                    end
                end
                close(upsert_channel)
                wait(publisher)
                is_cancellation_requested(scan_token) &&
                    throw(OperationCanceledException(scan_token))
                removals = String[id for id in keys(previous) if !(id in seen)]
                stale = count(
                    id_value -> haskey(previous, id_value) &&
                        !isequal(previous[id_value], current[id_value]),
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
                is_cancellation_requested(scan_token) &&
                    throw(OperationCanceledException(scan_token))
                # Final batch: removals + deferred metadata reconcile (not on every upsert batch).
                @timed_dbg publish_source_changes!(
                    workspace,
                    scan_epoch,
                    SourceChanges(AbstractDataSourceItem[], removals; metadata_changed=true),
                    status,
                )
                publish_scan_end!(
                    workspace, scan_epoch;
                    cache_hit=!any_upserts[] && isempty(removals),
                )
            catch error
                if error isa OperationCanceledException
                    publish_scan_end!(workspace, scan_epoch; canceled=true)
                else
                    publish_scan_end!(
                        workspace, scan_epoch; error, backtrace=catch_backtrace())
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

"""Return the stable IDs of all items currently published in a workspace."""
query_items(workspace::Workspace)::Vector{String} = sort!(collect(keys(workspace.index.items)))

"""Resize this workspace's DuckDB cache buffer-pool limit (MiB) live; suits a GUI control."""
set_cache_memory_limit!(workspace::Workspace, mib::Integer)::Int =
    set_cache_memory_limit!(workspace.cache.db, mib)

"""Replace the workspace item index with one complete collection record index."""
function replace_item_index!(
    workspace::Workspace,
    collections::CollectionIndex,
    records::Vector{ItemRecord},
)::Nothing
    items = Dict{String,ItemRecord}()
    sizehint!(items, length(records))
    items_by_source = Dict{Int64,Vector{String}}()
    metadata_keys = Set{Symbol}()
    for item in records
        haskey(items, item.id) && error("Duplicate item id generated during scan: $(item.id)")
        items[item.id] = item
        push!(get!(() -> String[], items_by_source, item.source_item_key), item.id)
    end
    for collection_record in values(collections.records)
        union!(metadata_keys, keys(collection_record.own_metadata))
        union!(metadata_keys, keys(collection_record.analysis))
    end
    item_metadata = Dict{String,Dict{Symbol,Any}}()
    for item in records
        previous = get(workspace.index.item_metadata, item.id, nothing)
        previous === nothing || (item_metadata[item.id] = previous)
    end
    workspace.index = WorkspaceIndex(
        collections,
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
    entity::Union{String,Int64},
)::Bool
    key_entity = kind in (COLLECTION_PROCESS_RESULT, COLLECTION_ANALYSIS_RESULT) ?
        entity::Int64 : entity::String
    state = get(index.result_states, CacheResultKey(kind, key_entity), nothing)
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

    collections = copy(index.source.collections)
    replace_item_index!(workspace, collections, copy(index.source.items))
    # The cache restores the computed metadata layer; source-table diff before seed skips stale work.
    workspace.index.item_metadata = Dict{String,Dict{Symbol,Any}}(
        id => copy(dict) for (id, dict) in index.item_metadata)
    workspace.index.source = SourceScan(
        index.source.source_id,
        index.source.source_label,
        collections,
        copy(index.source.items),
        index.source.analysis_failures,
    )

    workspace.index.analysis_errors = Dict{Union{String,Int64},String}(index.analysis_errors)
    for state in values(index.result_states)
        CacheResultStatus(state.status) === RESULT_FAILED || continue
        workspace.index.analysis_errors[state.entity] =
            something(state.message, "cached work failed")
    end
    if workspace.background_processing
        for record in values(workspace.index.items)
            cache_index_ready(index, PROCESSING_RESULT, record.id) ||
                enqueue_processing!(workspace, record)
            cache_index_ready(index, ITEM_ANALYSIS_RESULT, record.id) ||
                enqueue_item_analysis!(workspace, record)
        end
        for collection_key_value in keys(workspace.index.collections.records)
            member_ids = collection_item_ids(workspace.index.collections, collection_key_value)
            isempty(member_ids) && continue
            members = ItemRecord[
                workspace.index.items[id] for id in member_ids
                if haskey(workspace.index.items, id)
            ]
            collection_key = collection_key_value
            if !cache_index_ready(index, COLLECTION_PROCESS_RESULT, collection_key) ||
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
        haskey(workspace.index.collections.records, collection_key) || continue
        members = ItemRecord[
            workspace.index.items[id]
            for id in collection_item_ids(workspace.index.collections, collection_key)
            if haskey(workspace.index.items, id)
        ]
        isempty(members) && continue
        # Collection process reads processed payloads. Do not queue until every current member
        # has one; this function runs again when the next member finishes.
        all(record -> cache_work_status(
                workspace, WorkKey(ITEM_PROCESS, record.id)) === :ready, members) || continue
        member_analyze = WorkKey[
            WorkKey(ITEM_ANALYZE, record.id) for record in members]
        process_key = WorkKey(COLLECTION_PROCESS, collection_key)
        analyze_key = WorkKey(COLLECTION_ANALYZE, collection_key)
        if !supersede
            process_done = cache_work_status(workspace, process_key) === :ready
            analyze_done = cache_work_status(workspace, analyze_key) === :ready
            already_live = lock(workspace.work.lock) do
                get(workspace.work.nodes, process_key, nothing) !== nothing ||
                get(workspace.work.nodes, analyze_key, nothing) !== nothing
            end
            process_done && analyze_done && !already_live && continue
        end
        supersede && clear_work_result_state!(workspace, process_key)
        enqueue_work!(
            workspace, process_key, revision(process_key);
            priority=1, dependencies=member_analyze,
        )
        supersede && clear_work_result_state!(workspace, analyze_key)
        enqueue_work!(
            workspace, analyze_key, revision(analyze_key);
            priority=1,
            dependencies=WorkKey[WorkKey(COLLECTION_PROCESS, collection_key)],
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
        interpretation = result.interpretation::SourceItemInterpretation
        records = interpretation.records
        source_ref = source_item_id(workspace, key.entity::Int64)
        resolved, invalidated = @timed_dbg publish_source_item_records!(
                workspace, key.entity::Int64, records, interpretation.collection_paths)
        conflicts = store_interpreted!(
            workspace.cache.db,
            result.source_item,
            interpretation.source_item_label,
            resolved,
            interpretation.interpreted_items,
        )
        store_collection_index!(workspace.cache.db, workspace.index.collections, resolved)
        delete_collection_metadata!(workspace.cache.db, invalidated)
        delete_pruned_collection_records!(workspace, invalidated)
        enqueue_dependent_subtree!(workspace, resolved; collection_keys=invalidated)
        publish_metadata_conflicts!(workspace, source_ref, key.kind, conflicts)
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
        enqueue_collection_work!(
            workspace,
            affected_collection_keys(workspace.index.collections, [record]);
            supersede=false,
        )
    elseif key.kind === COLLECTION_PROCESS
        # Rewritten member payloads were persisted worker-side before publication.
        nothing
    else
        collection_key = key.entity::Int64
        haskey(workspace.index.collections.records, collection_key) || return nothing
        analysis = metadata_dict(result)
        set_collection_analysis!(workspace.index.collections, collection_key, analysis)
        dropped = store_collection_metadata!(
            workspace.cache.db,
            collection_key,
            analysis,
        )
        publish_metadata_conflicts!(workspace, collection_key, key.kind, dropped)
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
        workspace.cache.db,
        record,
        delivered_metadata(workspace, record, workspace.index.collections),
    )
    publish_metadata_conflicts!(workspace, record.id, ITEM_ANALYZE, dropped)
    return nothing
end

function publish_metadata_conflicts!(
    workspace::Workspace,
    entity::Union{String,Int64},
    stage,
    conflicts::Vector{String},
)::Nothing
    isempty(conflicts) && return nothing
    message = join(conflicts, "; ")
    workspace.index.analysis_errors[entity] = message
    @warn(
        "Workspace metadata conflict",
        entity,
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
    if failure.ex isa OperationCanceledException
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
        source_key = key.entity::Int64
        source_ref = source_item_id(workspace, source_key)
        old_records, invalidated = remove_source_item_output!(workspace, source_key)
        delete_source_item!(workspace.cache.db, source_key, old_records)
        delete_collection_metadata!(workspace.cache.db, invalidated)
        delete_pruned_collection_records!(workspace, invalidated)
        enqueue_dependent_subtree!(workspace, old_records; collection_keys=invalidated)
        workspace.index.analysis_errors[source_ref] = "interpret_source_item: " * message
        store_source_item_failure!(workspace.cache.db, source_key, message)
        @error(
            "Source interpretation failed",
            source_item=source_ref,
            stage=key.kind,
            exception=(failure.ex, failure.processed_bt),
        )
    else
        source_item_key_value = record === nothing ? Int64(0) : record.source_item_key
        cache_entity = key.kind in (COLLECTION_PROCESS, COLLECTION_ANALYZE) ?
            key.entity::Int64 : key.entity::String
        store_result_failure!(
            workspace.cache.db,
            _work_key_cache_kind(key),
            cache_entity,
            source_item_key_value,
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
            enqueue_collection_work!(
                workspace,
                affected_collection_keys(workspace.index.collections, [record]);
                supersede=false,
            )
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
    scan_epoch::Int,
    state::Symbol,
    index::Union{Nothing,ProjectCacheIndex},
)::Nothing
    lock(workspace.publish_lock) do
        scan_epoch == workspace.scan.epoch || return
        if index === nothing
            replace_item_index!(
                workspace,
                CollectionIndex(source_id(workspace.source)),
                ItemRecord[],
            )
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
    scan_epoch::Int,
    changes::SourceChanges,
    status::Union{Nothing,ProjectCacheStatus}=nothing,
)::Nothing
    lock(workspace.publish_lock) do
        scan_epoch == workspace.scan.epoch || return
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
        is_cancellation_requested(get_token(workspace.cancel_source)) && return
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
        is_cancellation_requested(get_token(workspace.cancel_source)) && return
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
    scan_epoch::Int;
    cache_hit::Bool=false,
    canceled::Bool=false,
    error=nothing,
    backtrace=nothing,
)::Nothing
    lock(workspace.publish_lock) do
        scan_epoch == workspace.scan.epoch || return
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
