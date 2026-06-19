"""
Start cache loading and source scanning for one new workspace.
"""
function open_workspace(project::Project, source::AbstractDataSource)::Workspace
    workspace = Workspace(project, source)
    load_cache!(workspace)
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
    cancel_job!(workspace.cache_job)
    for task in workspace.background_tasks
        istaskdone(task) || wait(task)
    end
    close_source!(workspace.source)
    workspace.closed = true
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
    workspace.cache_job.state in (:loading, :writing, :canceling)

analysis_work_running(workspace::Workspace)::Bool =
    workspace.analysis.state in (:analyzing, :canceling)

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
cancel_cache!(workspace::Workspace)::Nothing = (cancel_job!(workspace.cache_job); nothing)

"""
Start one source scan that progressively populates the workspace index.
"""
function scan_source!(workspace::Workspace)::Nothing
    source_scan_running(workspace) && error("A source scan is already running")
    workspace.closed && error("The workspace is closed")

    job = workspace.scan
    job.id += 1
    scan_id = job.id
    job.state = :discovering
    job.progress = WorkspaceProgress()
    job.error = ""
    events = Channel{NamedTuple}(Inf)
    cancel_token = Base.Threads.Atomic{Bool}(false)
    job.events = events
    job.cancel_token = cancel_token

    task = Base.Threads.@spawn begin
        try
            cached = try
                load_project_cache(workspace.cache.identity)
            catch
                nothing
            end
            source = with_cancel(() -> cancel_token[]) do
                scan_source(
                    workspace.project,
                    workspace.source;
                    cached_source=cached === nothing ? nothing : cached.source,
                    on_progress=(progress) -> put!(events, (
                        kind=:progress,
                        job_id=scan_id,
                        progress,
                    )),
                    on_items=(items) -> put!(events, (
                        kind=:items,
                        job_id=scan_id,
                        items=[ItemRecord(item) for item in items],
                    )),
                )
            end
            cache_hit = cached !== nothing && source === cached.source
            put!(events, (kind=:source, job_id=scan_id, source, cache_hit))
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
    track_task!(workspace, task)
    return nothing
end

"""
Start reading the compact cache index for a workspace.
"""
function load_cache!(workspace::Workspace)::Nothing
    cache_work_running(workspace) && error("Cache work is already running")
    workspace.closed && error("The workspace is closed")

    identity = workspace.cache.identity
    job = workspace.cache_job
    job.id += 1
    cache_id = job.id
    job.state = :loading
    job.progress = WorkspaceProgress()
    job.error = ""
    workspace.cache.index = nothing
    workspace.cache.status = nothing
    workspace.cache.operation = :load
    events = Channel{NamedTuple}(Inf)
    cancel_token = Base.Threads.Atomic{Bool}(false)
    job.events = events
    job.cancel_token = cancel_token

    task = Base.Threads.@spawn begin
        try
            index = with_cancel(() -> cancel_token[]) do
                load_project_cache(
                    identity;
                    on_progress=(progress) -> put!(events, (
                        kind=:progress,
                        job_id=cache_id,
                        progress,
                    )),
                )
            end
            put!(events, (kind=:cache, job_id=cache_id, index))
        catch error
            if is_job_cancelled(error)
                put!(events, (kind=:canceled, job_id=cache_id))
            elseif error isa ProjectCacheError
                put!(events, (
                    kind=:missing,
                    job_id=cache_id,
                    message=sprint(showerror, error),
                ))
            else
                put!(events, (
                    kind=:error,
                    job_id=cache_id,
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

"""Update the cache from the completed source scan already owned by the workspace."""
function update_cache!(
    workspace::Workspace;
    rebuild::Bool=false,
)::Nothing
    cache_work_running(workspace) && error("Cache work is already running")
    workspace.closed && error("The workspace is closed")
    source = workspace.index.source
    source isa SourceScan || error("Complete the source scan before updating the cache")

    identity = workspace.cache.identity
    previous_state = workspace.cache_job.state
    job = workspace.cache_job
    job.id += 1
    cache_id = job.id
    job.state = :writing
    job.progress = WorkspaceProgress()
    job.error = ""
    workspace.cache.operation =
        rebuild ? (previous_state == :missing ? :build : :rebuild) : :update
    events = Channel{NamedTuple}(Inf)
    cancel_token = Base.Threads.Atomic{Bool}(false)
    job.events = events
    job.cancel_token = cancel_token

    task = Base.Threads.@spawn begin
        try
            index = with_cancel(() -> cancel_token[]) do
                write_project_cache!(
                    identity,
                    source;
                    replace=rebuild,
                    on_progress=(progress) -> put!(events, (
                        kind=:progress,
                        job_id=cache_id,
                        progress,
                    )),
                )
            end
            put!(events, (kind=:cache, job_id=cache_id, index))
        catch error
            if is_job_cancelled(error)
                put!(events, (kind=:canceled, job_id=cache_id))
            else
                put!(events, (
                    kind=:error,
                    job_id=cache_id,
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
    workspace.index = WorkspaceIndex(
        hierarchy,
        items,
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
    item_index = copy(workspace.index.items)
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
    workspace.index.items = item_index
    workspace.index.collection_parameter_keys = sort!(collect(parameter_keys); by=String)
    return true
end

function _effective_loaded_item(record::ItemRecord, item::AbstractDataItem)::AbstractDataItem
    item isa DataItem || return item
    return DataItem(
        item.id,
        item.label,
        item.kind,
        item.collection,
        record.parameters,
        item.stats,
        item.data,
    )
end

function _load_for_analysis(
    workspace::Workspace,
    record::ItemRecord,
    loaded::Dict{String,AbstractDataItem},
)::AbstractDataItem
    cached = get(loaded, record.id, nothing)
    cached !== nothing && return cached
    item = load_data_item(workspace.project, workspace.source, record.source_item_id, record.id)
    item = _effective_loaded_item(record, item)
    loaded[record.id] = item
    return item
end

"""Start background per-item and collection stats for the completed source scan."""
function start_analysis!(workspace::Workspace)::Nothing
    analysis_work_running(workspace) && cancel_analysis!(workspace)
    source = workspace.index.source
    source isa SourceScan || return nothing

    job = workspace.analysis
    job.id += 1
    analysis_id = job.id
    total = length(source.hierarchy.all_items)
    job.state = :analyzing
    job.progress = WorkspaceProgress(phase=:analyzing, total_source_items=total)
    job.error = ""
    events = Channel{NamedTuple}(Inf)
    cancel_token = Base.Threads.Atomic{Bool}(false)
    job.events = events
    job.cancel_token = cancel_token

    task = Base.Threads.@spawn begin
        try
            item_stats = Dict{String,Dict{Symbol,Any}}()
            collection_node_stats = Dict{Tuple{Vararg{String}},Dict{Symbol,Any}}()
            failures = Dict{String,String}()
            loaded = Dict{String,AbstractDataItem}()
            processed = 0

            with_cancel(() -> cancel_token[]) do
                for record in source.hierarchy.all_items
                    check_cancel()
                    try
                        item = _load_for_analysis(workspace, record, loaded)
                        computed = compute_item_stats(workspace.project, workspace.source, item)
                        isempty(computed) || (item_stats[record.id] = computed)
                    catch error
                        is_job_cancelled(error) && rethrow()
                        failures[record.id] = "stats: " * sprint(showerror, error)
                    end
                    processed += 1
                    put!(events, (
                        kind=:progress,
                        job_id=analysis_id,
                        progress=(
                            phase=:analyzing,
                            total_source_items=total,
                            processed_source_items=processed,
                            loaded_items=processed,
                            skipped_source_items=0,
                            current_source_item=record.source_item_id,
                        ),
                    ))
                end

                for (path, node) in source.hierarchy.index
                    check_cancel()
                    isempty(node.items) && continue
                    try
                        items = AbstractDataItem[
                            _load_for_analysis(workspace, record, loaded)
                            for record in node.items
                        ]
                        computed = collection_stats(
                            workspace.project,
                            workspace.source,
                            collect(path),
                            items,
                        )
                        isempty(computed) || (collection_node_stats[path] = computed)
                    catch error
                        is_job_cancelled(error) && rethrow()
                        failures[first(node.items).id] =
                            "collection_stats: " * sprint(showerror, error)
                    end
                end
            end
            put!(events, (
                kind=:analysis,
                job_id=analysis_id,
                source_id=source.source_id,
                item_stats,
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
        record = get(workspace.index.items, id, nothing)
        record === nothing && continue
        merge!(record.stats, computed)
        changed = true
    end
    for (path, computed) in collection_node_stats
        node = get(workspace.index.hierarchy.index, path, nothing)
        node === nothing && continue
        merge!(node.stats, computed)
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

"""Apply one loaded cache index without replacing items already streamed by the source scan."""
function apply_cache_index!(
    workspace::Workspace,
    index::ProjectCacheIndex,
)::Bool
    index.identity.source_id == source_id(workspace.source) ||
        error("Loaded cache belongs to $(index.identity.source_id), not $(source_id(workspace.source))")
    workspace.cache.identity = index.identity
    workspace.cache.index = index

    index_changed = false
    source = workspace.index.source
    if source === nothing
        if isempty(workspace.index.items)
            replace_item_index!(workspace, index.source.hierarchy)
            index_changed = true
        else
            index_changed = append_items!(workspace, index.source.hierarchy.all_items)
        end
    end

    if source isa SourceScan
        workspace.cache.status = cache_status(index, source)
    else
        cached_items = length(index.source.source_item_fingerprints)
        workspace.cache.status = ProjectCacheStatus(
            cached_items,
            cached_items,
            0,
            0,
            0,
            0,
            length(index.analysis_errors),
        )
        workspace.index.analysis_errors = copy(index.analysis_errors)
    end
    workspace.cache_job.state = :ready
    workspace.cache_job.error = ""
    return index_changed
end

"""Apply the authoritative completed source scan to the workspace."""
function apply_source_scan!(
    workspace::Workspace,
    source::SourceScan;
    analyze::Bool=true,
)::Nothing
    source.source_id == source_id(workspace.source) ||
        error("Source scan belongs to $(source.source_id), not $(source_id(workspace.source))")
    workspace.index.source = source
    replace_item_index!(workspace, source.hierarchy)
    identity = workspace.cache.identity
    workspace.index.analysis_errors = ProjectCacheIndex(identity, source).analysis_errors

    cache_state = workspace.cache_job.state
    if cache_state == :ready
        workspace.cache.status = cache_status(workspace.cache.index, source)
    elseif !isfile(identity.cache_path) || cache_state == :missing
        source_item_count = length(source.source_item_fingerprints)
        workspace.cache.status = ProjectCacheStatus(
            source_item_count,
            0,
            0,
            0,
            source_item_count,
            0,
            0,
        )
    end
    analyze && start_analysis!(workspace)
    return nothing
end

"""Start cache repair when the completed source scan differs from the loaded cache."""
function repair_cache_if_needed!(workspace::Workspace)::Nothing
    cache_work_running(workspace) && return nothing
    analysis_work_running(workspace) && return nothing
    workspace.index.source isa SourceScan || return nothing
    if workspace.cache_job.state == :missing
        update_cache!(workspace; rebuild=true)
        return nothing
    end
    status = workspace.cache.status
    status isa ProjectCacheStatus || return nothing
    if status.stale_source_items > 0 ||
       status.new_source_items > 0 ||
       status.deleted_source_items > 0
        update_cache!(workspace)
    end
    return nothing
end

"""Apply all completed scan and cache events."""
function poll_workspace!(workspace::Workspace)::Bool
    index_changed = false
    cache_events = workspace.cache_job.events
    if cache_events !== nothing
        while isready(cache_events)
            event = take!(cache_events)
            event.job_id == workspace.cache_job.id || continue
            if event.kind == :progress
                workspace.cache_job.progress = WorkspaceProgress(event.progress)
            elseif event.kind == :cache
                index_changed |= apply_cache_index!(workspace, event.index)
                workspace.cache_job.events = nothing
                workspace.cache_job.cancel_token = nothing
            elseif event.kind == :missing
                workspace.cache_job.state = :missing
                workspace.cache_job.error = event.message
                workspace.cache_job.events = nothing
                workspace.cache_job.cancel_token = nothing
            elseif event.kind == :canceled
                workspace.cache_job.state = :canceled
                workspace.cache_job.events = nothing
                workspace.cache_job.cancel_token = nothing
            elseif event.kind == :error
                workspace.cache_job.state = :error
                workspace.cache_job.error =
                    "Cache job failed. See the console for full details."
                @error(
                    "Cache job failed",
                    cache_path=workspace.cache.identity.cache_path,
                    operation=workspace.cache.operation,
                    exception=(event.error, event.backtrace),
                )
                workspace.cache_job.events = nothing
                workspace.cache_job.cancel_token = nothing
            end
        end
    end

    scan_events = workspace.scan.events
    if scan_events !== nothing
        pending_items = ItemRecord[]
        while isready(scan_events)
            event = take!(scan_events)
            event.job_id == workspace.scan.id || continue
            if event.kind == :progress
                workspace.scan.progress = WorkspaceProgress(event.progress)
                workspace.scan.state = event.progress.phase
            elseif event.kind == :items
                append!(pending_items, event.items)
            elseif event.kind == :source
                apply_source_scan!(workspace, event.source; analyze=!event.cache_hit)
                index_changed = true
                workspace.scan.state = event.cache_hit ? :unchanged : :done
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
        index_changed |= append_items!(workspace, pending_items)
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

    repair_cache_if_needed!(workspace)
    return index_changed
end
