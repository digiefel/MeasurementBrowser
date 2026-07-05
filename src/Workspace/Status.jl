"""Fold a workspace's job and cache state into the watcher-facing [`WorkspaceStatus`](@ref)."""
function workspace_status(workspace::Workspace)::WorkspaceStatus
    errors = sort!(collect(workspace.index.analysis_errors); by=first)
    scan = workspace.scan

    # Active work: one live line, determinate bar when totals are known.
    if scan.state == :discovering
        label = workspace.cache.operation === :rebuild ? "Rebuilding" :
                workspace.cache.operation === :build ? "Building" : "Scanning"
        return WorkspaceStatus(
            :busy, label, "Finding source items…", true, nothing, errors)
    elseif scan.state == :canceling
        return WorkspaceStatus(:busy, "Canceling", "Canceling…", true, nothing, errors)
    elseif processing_work_running(workspace)
        completed, total, _active = work_counts(workspace)
        progress = total > 0 ? Float32(clamp(completed / total, 0, 1)) : nothing
        return WorkspaceStatus(
            :busy,
            "Caching",
            @sprintf("Processing %d/%d items", completed, total),
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
