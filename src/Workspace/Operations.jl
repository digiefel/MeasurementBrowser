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
)::Workspace
    workspace = Workspace(
        project,
        source;
        profile_internal,
        profile_cpu,
        profile_output,
        crash_trace,
    )
    scan_source!(workspace)
    return workspace
end

"""
Cancel all work owned by a workspace and wait for it to stop.
"""
function close_workspace!(workspace::Workspace)::Nothing
    workspace.closed && return nothing
    cancel_job!(workspace.scan)
    cancel_job!(workspace.analysis)
    for task in workspace.background_tasks
        istaskdone(task) || wait(task)
    end
    stop_processing_workers!(workspace)
    # Workers are gone, so no more deposits can arrive; drain everything still staged and stop the
    # flusher before the database handle closes.
    stop_cache_buffer!(workspace.cache.db.buffer)
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
    workspace.analysis.state in (:analyzing, :canceling)

"""Whether any item is waiting for or running processing, or its result is still being written."""
function processing_work_running(workspace::Workspace)::Bool
    jobs_active = lock(workspace.processing.lock) do
        !isempty(workspace.processing.jobs)
    end
    return jobs_active || buffer_has_pending_writes(workspace.cache.db.buffer)
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
cancel_analysis!(workspace::Workspace)::Nothing = (cancel_job!(workspace.analysis); nothing)

"""
Cancel in-flight cache work.

Cache writes are now produced by the source scan itself, so cancelling cache work cancels the scan.
"""
cancel_cache!(workspace::Workspace)::Nothing = cancel_scan!(workspace)

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
    analysis_work_running(workspace) && cancel_analysis!(workspace)

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
            cached = if rebuild
                nothing
            else
                cache_built = with_reader(cachedb) do connection
                    only(DBInterface.execute(
                        connection,
                        "SELECT count(*) > 0 AS built FROM meta WHERE key = 'schema_version'",
                    )).built
                end
                cache_built ? load_cache_index(cachedb) : nothing
            end
            # The UI owns `cached` (it displays it and streams progressive items into it). The scan
            # reuses unchanged source items by record object, then rewrites their effective parameters
            # while building its result — so it must work on its OWN deep copy, or it would mutate the
            # very records the render thread is mutating, racing on each record's parameter dict. The
            # copy is taken here, before `cached` is handed to the UI, so nothing else touches it yet.
            reuse_index = cached === nothing ? nothing : deepcopy(cached)
            put!(events, (
                kind=:cache_state,
                job_id=scan_id,
                state=cached === nothing ? :missing : :ready,
                index=cached,
            ))

            # Establish cache identity before the scan transaction writes its replacement contents.
            write_scan_identity!(cachedb)
            wrote = Ref(false)
            written_ids = Set{String}()
            kept_ids = Set{String}()

            source = with_cancel(() -> cancel_token[]) do
                scan_source(
                    workspace.project,
                    workspace.source;
                    cached_source=reuse_index === nothing ? nothing : reuse_index.source,
                    on_progress=(progress) -> put!(events, (
                        kind=:progress,
                        job_id=scan_id,
                        progress,
                    )),
                    on_items=(items) -> begin
                        records = ItemRecord[ItemRecord(item) for item in items]
                        put!(events, (kind=:items, job_id=scan_id, items=records))
                    end,
                    on_source_item=(interpretation::SourceItemInterpretation) -> begin
                        records = interpretation.records
                        # Hand the buffer every interpreted item, cacheable or not — it serves them
                        # from memory and the flush persists only the cacheable ones. Discarding the
                        # non-cacheable items here was what forced a redundant whole-file re-read.
                        data = Any[item for item in interpretation.interpreted_items]
                        wrote[] || put!(events, (
                            kind=:cache_writing,
                            job_id=scan_id,
                        ))
                        wrote[] = true
                        # Deposit into the buffer (never blocks on the database) and enqueue the
                        # data-less records. Processing reads the interpreted data straight back from
                        # the buffer's staged store, before it is even durable.
                        store_interpreted!(workspace.cache.db, records, data)
                        for record in records
                            push!(written_ids, record.source_item_id)
                            enqueue_processing!(workspace, record)
                        end
                    end,
                    on_kept_source_item=(source_item_id::String) ->
                        push!(kept_ids, source_item_id),
                    # The opening move of an update: drop source items the source no longer contains,
                    # while the buffers are still empty so no append can race the delete.
                    on_discovered=(present_ids) ->
                        retain_source_items!(cachedb, present_ids),
                    on_failure=(failure::ItemFailure) ->
                        put!(events, (kind=:failure, job_id=scan_id, failure)),
                    profiler=workspace.profiler,
                )
            end
            # Drain the per-item deposits before the readiness probe below queries the cache for
            # already-processed items.
            wait_cache_flushed!(workspace.cache.db.buffer)
            # Unchanged source items keep their existing rows: treat them as already written so the
            # collection-level append below records no redundant bare row for them.
            union!(written_ids, kept_ids)
            cache_hit = reuse_index !== nothing && source === reuse_index.source
            cache_hit || store_scan_collection_data!(cachedb, source, written_ids)
            # A fresh build already enqueued every record from `on_source_item`, so the readiness
            # probe below would query DuckDB for thousands of rows only to learn nothing is cached.
            # It earns its keep only on an incremental reopen, where unchanged source items were kept
            # (never enqueued during the scan) and may still need processing.
            if reuse_index !== nothing
                records = source.hierarchy.all_items
                processed_ids = cached_item_data_ids(cachedb, records; stage=:processed)
                for record in records
                    record.id in processed_ids || enqueue_processing!(workspace, record)
                end
            end
            put!(events, (
                kind=:source,
                job_id=scan_id,
                source,
                cache_hit,
                status=_post_write_status(workspace.cache.identity, source, reuse_index),
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
    completed, total, jobs = lock(workspace.processing.lock) do
        (
            workspace.processing.completed,
            workspace.processing.total,
            length(workspace.processing.jobs),
        )
    end
    depth = jobs + buffer_pending_counts(workspace.cache.db.buffer).items
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
        cancel_waiting_processing!(workspace)
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

"""Start background collection stats for the completed source scan."""
function start_analysis!(workspace::Workspace)::Nothing
    analysis_work_running(workspace) && cancel_analysis!(workspace)
    source = workspace.index.source
    source isa SourceScan || return nothing

    job = workspace.analysis
    job.id += 1
    analysis_id = job.id
    total = count(node -> !isempty(node.items), values(source.hierarchy.index))
    job.state = :analyzing
    job.progress = WorkspaceProgress(phase=:analyzing, total_source_items=total)
    job.error = ""
    events = Channel{NamedTuple}(Inf)
    cancel_token = Base.Threads.Atomic{Bool}(false)
    job.events = events
    job.cancel_token = cancel_token

    cachedb = workspace.cache.db
    task = Base.Threads.@spawn begin
        try
            collection_node_stats = Dict{Tuple{Vararg{String}},Dict{Symbol,Any}}()
            no_item_stats = Dict{String,Dict{Symbol,Any}}()
            failures = Dict{String,String}()
            processed = 0

            with_cancel(() -> cancel_token[]) do
                for (path, node) in source.hierarchy.index
                    check_cancel()
                    isempty(node.items) && continue
                    try
                        items = AbstractDataItem[
                            DataItem(
                                ItemRecord(record; stats=get(
                                    workspace.index.item_stats, record.id, record.stats)),
                                nothing,
                            )
                            for record in node.items
                        ]
                        computed = Profiling.@profile_span workspace.profiler :project :collection_stats Profiling.ProfileAttributes(
                            items=Int64(length(items)),
                        ) begin
                            collection_stats(
                                workspace.project,
                                workspace.source,
                                collect(path),
                                items,
                            )
                        end
                        isempty(computed) || (collection_node_stats[path] = computed)
                    catch error
                        is_job_cancelled(error) && rethrow()
                        failures[first(node.items).id] =
                            "collection_stats: " * sprint(showerror, error)
                    end
                    processed += 1
                    put!(events, (
                        kind=:progress,
                        job_id=analysis_id,
                        progress=(
                            phase=:analyzing,
                            total_source_items=total,
                            processed_source_items=processed,
                            loaded_items=length(node.items),
                            skipped_source_items=0,
                            current_source_item=first(node.items).source_item_id,
                        ),
                    ))
                end
                if !isempty(collection_node_stats)
                    write_started = time_ns()
                    store_stats!(cachedb, no_item_stats, collection_node_stats)
                    record_cache_phase!(
                        workspace.metrics.stats_write_ns,
                        workspace.metrics.stats_writes,
                        write_started,
                    )
                end
            end

            put!(events, (
                kind=:analysis,
                job_id=analysis_id,
                source_id=source.source_id,
                item_stats=Dict{String,Dict{Symbol,Any}}(),
                collection_stats=collection_node_stats,
                failures,
            ))
        catch error
            if is_job_cancelled(error)
                put!(events, (kind=:canceled, job_id=analysis_id))
            else
                put!(events, (
                    kind=:error,
                    job_id=analysis_id,
                    error,
                    backtrace=catch_backtrace(),
                ))
            end
        finally
            close(events)
        end
    end
    track_task!(workspace, task)
    return nothing
end

"""Merge completed background stats into the live hierarchy."""
function apply_analysis!(
    workspace::Workspace,
    source_id_value::String,
    item_stats::Dict{String,Dict{Symbol,Any}},
    collection_node_stats::Dict{Tuple{Vararg{String}},Dict{Symbol,Any}},
    failures::Dict{String,String},
)::Bool
    source_id_value == source_id(workspace.source) ||
        error("Analysis belongs to $(source_id_value), not $(source_id(workspace.source))")
    changed = false
    for (id, computed) in item_stats
        haskey(workspace.index.items, id) || continue
        workspace.index.item_stats[id] = metadata_dict(computed)
        changed = true
    end
    for (path, computed) in collection_node_stats
        node = get(workspace.index.hierarchy.index, path, nothing)
        node === nothing && continue
        merge!(node.stats, metadata_dict(computed))
        changed = true
    end
    for (id, message) in failures
        workspace.index.analysis_errors[id] = message
    end
    source = workspace.index.source
    if source isa SourceScan
        for (id, message) in failures
            record = get(workspace.index.items, id, nothing)
            push!(source.analysis_failures, ItemFailure(
                record === nothing ? "" : record.source_item_id,
                id,
                message,
            ))
        end
    end
    return changed || !isempty(failures)
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

    cached_items = length(index.source.source_item_fingerprints)
    workspace.cache.status = ProjectCacheStatus(
        cached_items, cached_items, 0, 0, 0, 0, length(index.analysis_errors))
    workspace.index.analysis_errors = copy(index.analysis_errors)
    return nothing
end

"""Apply the authoritative completed source scan (already mirrored into the cache) to the workspace."""
function apply_source_scan!(
    workspace::Workspace,
    source::SourceScan,
    status::ProjectCacheStatus;
    analyze::Bool=true,
)::Nothing
    source.source_id == source_id(workspace.source) ||
        error("Source scan belongs to $(source.source_id), not $(source_id(workspace.source))")
    workspace.index.source = source
    replace_item_index!(workspace, source.hierarchy)
    identity = workspace.cache.identity
    source_errors = ProjectCacheIndex(identity, source).analysis_errors
    item_ids = keys(workspace.index.items)
    processing_errors = Dict(
        id => message
        for (id, message) in workspace.index.analysis_errors
        if id in item_ids
    )
    merge!(source_errors, processing_errors)
    workspace.index.analysis_errors = source_errors
    workspace.cache.status = status
    analyze && (workspace.analysis.state = :pending)
    return nothing
end

"""Apply all completed scan and cache events."""
function poll_workspace!(workspace::Workspace)::Bool
    index_changed = false

    scan_events = workspace.scan.events
    if scan_events !== nothing
        pending_items = ItemRecord[]
        while isready(scan_events)
            event = take!(scan_events)
            event.job_id == workspace.scan.id || continue
            if event.kind == :progress
                workspace.scan.progress = WorkspaceProgress(event.progress)
                workspace.scan.state = event.progress.phase
            elseif event.kind == :failure
                # Source-item failures stream in as they happen, so a long scan shows current errors
                # instead of revealing them all at the end. apply_source_scan! later replaces this with
                # the authoritative set.
                failure = event.failure
                previous = get(workspace.index.analysis_errors, failure.source_item_id, "")
                workspace.index.analysis_errors[failure.source_item_id] =
                    isempty(previous) ? failure.message : previous * "\n" * failure.message
                index_changed = true
            elseif event.kind == :items
                # Progressive first-view items only matter until a complete scan is applied. Once
                # `index.source` is set, a background analysis task iterates that hierarchy, so
                # mutating it here (insert_item!/apply_collection_parameters!) would race that task
                # and corrupt the heap. The next :source replaces the whole index anyway.
                workspace.index.source === nothing && append!(pending_items, event.items)
            elseif event.kind == :cache_state
                if event.index === nothing
                    # Rebuild, or first build with no cache: clear to a fresh empty tree so the scan
                    # streams the new hierarchy in from scratch.
                    replace_item_index!(workspace, Hierarchy(source_id(workspace.source), false))
                    empty!(workspace.index.analysis_errors)
                else
                    # The cache found on disk, surfaced as a complete, independent snapshot before fresh
                    # results refine it.
                    apply_cache_index!(workspace, event.index)
                end
                index_changed = true
                workspace.cache_state = event.state
                event.state == :missing && workspace.cache.operation != :rebuild &&
                    (workspace.cache.operation = :build)
            elseif event.kind == :cache_writing
                workspace.cache_state = :writing
            elseif event.kind == :source
                # No progressive append here: apply_source_scan! rebuilds the index wholesale from
                # event.source, which already contains every record accumulated in pending_items.
                empty!(pending_items)
                apply_source_scan!(
                    workspace, event.source, event.status; analyze=!event.cache_hit)
                index_changed = true
                workspace.scan.state = event.cache_hit ? :unchanged : :done
                workspace.cache_state = :ready
                workspace.scan.events = nothing
                workspace.scan.cancel_token = nothing
            elseif event.kind == :canceled
                workspace.scan.state = :canceled
                workspace.cache_state =
                    workspace.cache_state == :writing ? :canceled : workspace.cache_state
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
        workspace.index.source === nothing &&
            (index_changed |= append_items!(workspace, pending_items))
    end

    processing_events = workspace.processing.events
    while isready(processing_events)
        event = take!(processing_events)
        if event.kind == :stats
            workspace.index.item_stats[event.item_id] = metadata_dict(event.stats)
            delete!(workspace.index.analysis_errors, event.item_id)
        elseif event.kind == :failure
            workspace.index.analysis_errors[event.item_id] = event.message
        else
            error("Unknown processing event kind $(event.kind)")
        end
        index_changed = true
    end

    analysis_events = workspace.analysis.events
    if analysis_events !== nothing
        while isready(analysis_events)
            event = take!(analysis_events)
            event.job_id == workspace.analysis.id || continue
            if event.kind == :progress
                workspace.analysis.progress = WorkspaceProgress(event.progress)
            elseif event.kind == :analysis
                index_changed |= apply_analysis!(
                    workspace,
                    event.source_id,
                    event.item_stats,
                    event.collection_stats,
                    event.failures,
                )
                workspace.analysis.state = :done
                workspace.analysis.events = nothing
                workspace.analysis.cancel_token = nothing
            elseif event.kind == :canceled
                workspace.analysis.state = :canceled
                workspace.analysis.events = nothing
                workspace.analysis.cancel_token = nothing
            elseif event.kind == :error
                workspace.analysis.state = :error
                workspace.analysis.error =
                    "Analysis failed. See the console for full details."
                @error(
                    "Analysis job failed",
                    source=source_id(workspace.source),
                    exception=(event.error, event.backtrace),
                )
                workspace.analysis.events = nothing
                workspace.analysis.cancel_token = nothing
            end
        end
    end

    !workspace.profile_restart_pending && workspace.analysis.state == :pending &&
        !processing_work_running(workspace) &&
        start_analysis!(workspace)

    # Refresh the watcher snapshot while work runs (progress advances without an index change) and for
    # the one frame after it stops; idle frames skip the rebuild entirely.
    busy = source_scan_running(workspace) || processing_work_running(workspace) ||
           analysis_work_running(workspace) ||
           cache_work_running(workspace)
    if workspace.profile_restart_pending && !busy
        reset_processing_queue!(workspace)
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
