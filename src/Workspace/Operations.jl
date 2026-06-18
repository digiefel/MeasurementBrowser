"""
Start cache loading and source scanning for one new workspace.
"""
function open_workspace(project::Project, source::AbstractDataSource)::Workspace
    workspace = Workspace(project, source)
    load_cache!(workspace)
    scan_source!(workspace)
    return workspace
end

open_workspace(project::Project, root_path::AbstractString)::Workspace =
    open_workspace(project, RegisteredProjectSource(project, root_path))

"""
Cancel all work owned by a workspace and wait for it to stop.
"""
function close_workspace!(workspace::Workspace)::Nothing
    workspace.closed && return nothing
    cancel_job!(workspace.scan)
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
    workspace.selection.item_keys = [item_record_key(item) for item in items]
    return nothing
end

"""Replace the selection with exact workspace item keys."""
function select_items!(
    workspace::Workspace,
    item_keys::AbstractVector{<:AbstractString},
)::Nothing
    keys = String[item_key for item_key in item_keys]
    missing_keys = String[item_key for item_key in keys if !haskey(workspace.index.items, item_key)]
    isempty(missing_keys) || error(
        "Cannot select $(length(missing_keys)) item key(s) that are not in this workspace: " *
        join(missing_keys, ", "),
    )
    workspace.selection.item_keys = keys
    return nothing
end

"""Replace the selection with data items whose `item_id` is unambiguous in the workspace."""
function select_items!(
    workspace::Workspace,
    items::AbstractVector{<:AbstractDataItem},
)::Nothing
    keys = String[]
    for item in items
        id = item_id(item)
        # ponytail: linear scan per requested item; add an item_id index if bulk selection gets slow.
        matches = [record for record in values(workspace.index.items) if record.item_id == id]
        isempty(matches) && error(
            "Cannot select item_id '$id': no indexed item with that id exists in this workspace",
        )
        length(matches) == 1 || error(
            "Cannot select item_id '$id' because it matches $(length(matches)) indexed items; " *
            "select by ItemRecord or exact workspace item key instead",
        )
        push!(keys, item_record_key(only(matches)))
    end

    workspace.selection.item_keys = keys
    return nothing
end

"""Clear the selection, or fail clearly for unsupported item selectors."""
function select_items!(
    workspace::Workspace,
    items::AbstractVector,
)::Nothing
    if isempty(items)
        empty!(workspace.selection.item_keys)
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
        "exact workspace item keys, or AbstractDataItem values",
    )
end

source_scan_running(workspace::Workspace)::Bool =
    workspace.scan.state in (:counting, :discovering, :scanning, :analyzing, :canceling)

cache_work_running(workspace::Workspace)::Bool =
    workspace.cache_job.state in (:loading, :writing, :canceling)

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
    metadata_keys = Set{Symbol}()
    for item in hierarchy.all_items
        key = item_record_key(item)
        haskey(items, key) && error("Duplicate item key generated during scan: $(key)")
        items[key] = item
        union!(metadata_keys, keys(item.collection_metadata))
    end
    workspace.index = WorkspaceIndex(
        hierarchy,
        items,
        sort!(collect(metadata_keys); by=String),
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
    metadata_keys = Set(workspace.index.collection_metadata_keys)
    changed = false
    for item in items
        key = item_record_key(item)
        haskey(item_index, key) && continue
        insert_item!(hierarchy, item)
        item_index[key] = item
        union!(metadata_keys, keys(item.collection_metadata))
        changed = true
    end
    changed || return false
    workspace.index.collection_metadata_keys = sort!(collect(metadata_keys); by=String)
    return true
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
    source::SourceScan,
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
    return nothing
end

"""Start cache repair when the completed source scan differs from the loaded cache."""
function repair_cache_if_needed!(workspace::Workspace)::Nothing
    cache_work_running(workspace) && return nothing
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
                apply_source_scan!(workspace, event.source)
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

    repair_cache_if_needed!(workspace)
    return index_changed
end
