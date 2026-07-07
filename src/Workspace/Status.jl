"""Fold a workspace's job and cache state into the watcher-facing [`WorkspaceStatus`](@ref)."""
function workspace_status(workspace::Workspace)::WorkspaceStatus
    # Never wait on publish_lock here: a source-item publish or cache flush can hold it for 100ms+,
    # and this runs every GUI frame. If a publish is mid-flight, reuse the last snapshot's errors —
    # they change rarely, so a one-frame delay is invisible and keeps the frame from stalling.
    errors = if trylock(workspace.publish_lock)
        try
            lock(workspace.work.lock) do
                sort!(Pair{String,String}[k => v for (k, v) in workspace.index.analysis_errors]; by=first)
            end
        finally
            unlock(workspace.publish_lock)
        end
    else
        workspace.status.errors
    end
    scan = workspace.scan

    # Active work: one live line, determinate bar when totals are known.
    if scan.state == :discovering
        label = workspace.cache.operation === :rebuild ? "Rebuilding" :
                workspace.cache.operation === :build ? "Building" : "Scanning"
        return WorkspaceStatus(
            :busy, label, "Finding source items…", true, nothing, errors)
    elseif scan.state == :canceling
        return WorkspaceStatus(:busy, "Canceling", "Canceling…", true, nothing, errors)
    end
    completed, total, active = work_counts(workspace)
    if active > 0
        progress = total > 0 ? Float32(clamp(completed / total, 0, 1)) : nothing
        return WorkspaceStatus(
            :busy,
            "Caching",
            @sprintf("Processing %d/%d tasks", completed, total),
            true,
            progress,
            errors,
        )
    end

    # Settled: a single merged summary line, no bar.
    workspace.cache.identity === nothing && return WorkspaceStatus(
        :none, "No Project", "Open a project folder to build a cache.", false, nothing, errors)
    scan.state == :error && return WorkspaceStatus(
        :error, "Scan Error", scan.error, false, nothing, errors)
    isempty(workspace.source_error) || return WorkspaceStatus(
        :error, "Source Error", workspace.source_error, false, nothing, errors)
    workspace.cache_state == :error && return WorkspaceStatus(
        :error, "Cache Error", workspace.cache_error, false, nothing, errors)
    workspace.cache_state == :missing && return WorkspaceStatus(
        :missing, "Cache Missing", "No cache has been built for this project yet.",
        false, nothing, errors)
    scan.state == :canceled && return WorkspaceStatus(
        :none, "Canceled", "The last scan was canceled.", false, nothing, errors)

    status = workspace.cache.status
    status isa ProjectCacheStatus || return WorkspaceStatus(
        :none, "Idle", "No cache loaded.", false, nothing, errors)
    item_count = length(workspace.index.items)
    if !isempty(errors) || status.error_source_items > 0
        failed = max(length(errors), status.error_source_items)
        return WorkspaceStatus(:error, "Errors",
            @sprintf("%d source item(s) failed · %d items cached", failed, item_count),
            false, nothing, errors)
    elseif !(workspace.index.source isa SourceScan)
        return WorkspaceStatus(:fresh, "Loaded",
            @sprintf("%d source items cached (source not checked)", status.cached_source_items),
            false, nothing, errors)
    end
    # A settled successful scan has ingested every difference it found, so the cache mirrors the
    # source; the scan's stale/new/deleted counts describe the work it did, not pending work.
    changed = status.stale_source_items + status.new_source_items + status.deleted_source_items
    detail = changed > 0 ?
        @sprintf("%d source items · %d items cached · last scan: %d stale · %d new · %d deleted",
            status.total_source_items, item_count, status.stale_source_items,
            status.new_source_items, status.deleted_source_items) :
        @sprintf("%d source items · %d items cached", status.total_source_items, item_count)
    return WorkspaceStatus(:fresh, "Fresh", detail, false, nothing, errors)
end
