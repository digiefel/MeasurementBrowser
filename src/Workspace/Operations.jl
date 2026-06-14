"""
Start cache loading and source scanning for one new workspace.
"""
function open_workspace(
    project::AbstractProject,
    root_path::AbstractString,
)::Workspace
    workspace = Workspace(project, root_path)
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
    cancel_job!(workspace.cache_job)
    for task in workspace.background_tasks
        istaskdone(task) || wait(task)
    end
    workspace.closed = true
    return nothing
end

"""
Replace the selected measurements with the supplied measurement identities.
"""
function select_measurements!(
    workspace::Workspace,
    measurements::Vector{MeasurementInfo},
)::Nothing
    workspace.selection.measurement_ids =
        [measurement.unique_id for measurement in measurements]
    return nothing
end

source_scan_running(workspace::Workspace)::Bool =
    workspace.scan.state in (:counting, :discovering, :scanning, :analyzing, :canceling)

cache_work_running(workspace::Workspace)::Bool =
    workspace.cache_job.state in (:loading, :writing, :canceling)

"""
Keep a background task alive until workspace shutdown.
"""
function track_task!(workspace::Workspace, task::Task)::Task
    filter!(!istaskdone, workspace.background_tasks)
    push!(workspace.background_tasks, task)
    return task
end

"""
Request cancellation of one running workspace operation.
"""
function cancel_job!(job::WorkspaceJob)::Nothing
    job.cancel_token === nothing && return nothing
    Base.Threads.atomic_xchg!(job.cancel_token, true)
    job.state = :canceling
    return nothing
end

"""
Request cancellation of the current source scan.
"""
function cancel_scan!(workspace::Workspace)::Nothing
    cancel_job!(workspace.scan)
    return nothing
end

"""
Request cancellation of the current cache operation.
"""
function cancel_cache!(workspace::Workspace)::Nothing
    cancel_job!(workspace.cache_job)
    return nothing
end

"""
Start one filesystem scan that progressively populates the workspace index.
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
            source = with_cancel(() -> cancel_token[]) do
                scan_source(
                    workspace.root_path;
                    project=workspace.project,
                    on_progress=(progress) -> put!(events, (
                        kind=:progress,
                        job_id=scan_id,
                        progress,
                    )),
                    # Stream immutable snapshots for incremental display. The scan keeps mutating the
                    # original measurements' stats on worker threads, so the rendering thread must
                    # never see those objects mid-mutation (that race corrupts the heap). The complete
                    # originals replace these snapshots atomically when the :source event arrives.
                    on_measurements=(measurements) -> put!(events, (
                        kind=:measurements,
                        job_id=scan_id,
                        measurements=[MeasurementInfo(m) for m in measurements],
                    )),
                )
            end
            put!(events, (kind=:source, job_id=scan_id, source))
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

"""
Update the cache from the completed source scan already owned by the workspace.
"""
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

"""
Replace the workspace measurement index with one complete hierarchy.
"""
function replace_measurement_index!(
    workspace::Workspace,
    hierarchy::MeasurementHierarchy,
)::Nothing
    measurements = Dict{String,MeasurementInfo}()
    sizehint!(measurements, length(hierarchy.all_measurements))
    metadata_keys = Set{Symbol}()
    for measurement in hierarchy.all_measurements
        haskey(measurements, measurement.unique_id) &&
            error("Duplicate measurement id generated during scan: $(measurement.unique_id)")
        measurements[measurement.unique_id] = measurement
        union!(metadata_keys, keys(measurement.device_info.parameters))
    end
    workspace.index = WorkspaceIndex(
        hierarchy,
        measurements,
        sort!(collect(metadata_keys); by=String),
        workspace.index.source,
        workspace.index.analysis_errors,
    )
    return nothing
end

"""
Add newly discovered measurements to the provisional workspace index.

The completed source scan later replaces this hierarchy with the final sorted hierarchy.
"""
function append_measurements!(
    workspace::Workspace,
    measurements::Vector{MeasurementInfo},
)::Bool
    isempty(measurements) && return false
    hierarchy = workspace.index.hierarchy
    measurement_index = workspace.index.measurements
    metadata_keys = Set(workspace.index.device_metadata_keys)
    changed = false
    for measurement in measurements
        haskey(measurement_index, measurement.unique_id) && continue
        insert_measurement!(hierarchy, measurement)
        measurement_index[measurement.unique_id] = measurement
        union!(metadata_keys, keys(measurement.device_info.parameters))
        changed = true
    end
    changed || return false
    workspace.index.device_metadata_keys = sort!(collect(metadata_keys); by=String)
    return true
end

"""
Apply one loaded cache index without replacing measurements already streamed by the source scan.
"""
function apply_cache_index!(
    workspace::Workspace,
    index::ProjectCacheIndex,
)::Bool
    index.identity.root_path == workspace.root_path ||
        error("Loaded cache belongs to $(index.identity.root_path), not $(workspace.root_path)")
    workspace.cache.identity = index.identity
    workspace.cache.index = index

    index_changed = false
    source = workspace.index.source
    if source === nothing
        if isempty(workspace.index.measurements)
            replace_measurement_index!(workspace, index.source.hierarchy)
            index_changed = true
        else
            index_changed =
                append_measurements!(workspace, index.source.hierarchy.all_measurements)
        end
    end

    if source isa SourceScan
        workspace.cache.status = cache_status(index, source)
    else
        cached_files = length(index.files)
        workspace.cache.status = ProjectCacheStatus(
            cached_files,
            cached_files,
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

"""
Apply the authoritative completed source scan to the workspace.
"""
function apply_source_scan!(
    workspace::Workspace,
    source::SourceScan,
)::Nothing
    source.root_path == workspace.root_path ||
        error("Source scan belongs to $(source.root_path), not $(workspace.root_path)")
    workspace.index.source = source
    replace_measurement_index!(workspace, source.hierarchy)
    identity = workspace.cache.identity
    workspace.index.analysis_errors = ProjectCacheIndex(identity, source).analysis_errors

    cache_state = workspace.cache_job.state
    if cache_state == :ready
        workspace.cache.status = cache_status(workspace.cache.index, source)
    elseif !isfile(identity.cache_path) || cache_state == :missing
        workspace.cache.status = ProjectCacheStatus(
            length(source.files),
            0,
            0,
            0,
            length(source.files),
            0,
            0,
        )
    end
    return nothing
end

"""
Start cache repair when the completed source scan differs from the loaded cache.
"""
function repair_cache_if_needed!(workspace::Workspace)::Nothing
    cache_work_running(workspace) && return nothing
    workspace.index.source isa SourceScan || return nothing
    if workspace.cache_job.state == :missing
        update_cache!(workspace; rebuild=true)
        return nothing
    end
    status = workspace.cache.status
    status isa ProjectCacheStatus || return nothing
    if status.stale_files > 0 || status.new_files > 0 || status.deleted_files > 0
        update_cache!(workspace)
    end
    return nothing
end

"""
Apply all completed scan and cache events.

Returns `true` when the measurement index changed and browser references must be refreshed.
"""
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
        pending_measurements = MeasurementInfo[]
        while isready(scan_events)
            event = take!(scan_events)
            event.job_id == workspace.scan.id || continue
            if event.kind == :progress
                workspace.scan.progress = WorkspaceProgress(event.progress)
                workspace.scan.state = event.progress.phase
            elseif event.kind == :measurements
                append!(pending_measurements, event.measurements)
            elseif event.kind == :source
                apply_source_scan!(workspace, event.source)
                index_changed = true
                workspace.scan.state = :done
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
                    project=project_name(workspace.project),
                    root=workspace.root_path,
                    current_file=workspace.scan.progress.current_path,
                    exception=(event.error, event.backtrace),
                )
                workspace.scan.events = nothing
                workspace.scan.cancel_token = nothing
            end
        end
        index_changed |= append_measurements!(workspace, pending_measurements)
    end

    repair_cache_if_needed!(workspace)
    return index_changed
end
