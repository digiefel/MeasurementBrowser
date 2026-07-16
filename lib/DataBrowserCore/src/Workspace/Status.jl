"""Fold a workspace's job and cache state into the watcher-facing [`WorkspaceStatus`](@ref)."""
function workspace_status(workspace::Workspace)::WorkspaceStatus
    # Never wait on publish_lock here: a source-item publish or cache flush can hold it for 100ms+,
    # and this runs every GUI frame. If a publish is mid-flight, reuse the last snapshot's errors —
    # they change rarely, so a one-frame delay is invisible and keeps the frame from stalling.
    errors = if trylock(workspace.publish_lock)
        try
            lock(workspace.work.lock) do
                sort!(
                    Pair{String,String}[string(k) => v for (k, v) in workspace.index.analysis_errors];
                    by=first,
                )
            end
        finally
            unlock(workspace.publish_lock)
        end
    else
        workspace.status.errors
    end
    scan = workspace.scan
    counts = workspace_stage_counts(workspace)
    # During an incremental update, keep the last settled verdict on the chip instead of flashing
    # colors. A rebuild or first build voids any verdict — the cache was just discarded or has
    # never existed — and without one nothing is known fresh yet, so the chip says Building.
    settled = workspace.cache.operation === :update &&
        workspace.status.label in ("Fresh", "Loaded", "Errors")
    calm_level = settled ? workspace.status.level : :busy
    calm_label = settled ? workspace.status.label : "Building"

    # Active work: one live line, determinate bar when totals are known.
    if scan.state == :discovering
        noun = source_item_noun(workspace.source)
        return WorkspaceStatus(
            calm_level,
            calm_label,
            @sprintf("Finding %s… (%d found)", noun, counts.sources_found),
            true,
            nothing,
            counts,
            errors,
        )
    elseif scan.state == :canceling
        return WorkspaceStatus(
            calm_level, "Canceling", "Canceling…", true, nothing, counts, errors)
    end
    completed, total, active = work_counts(workspace)
    if active > 0
        progress = total > 0 ? Float32(clamp(completed / total, 0, 1)) : nothing
        return WorkspaceStatus(
            calm_level,
            calm_label,
            @sprintf("Processing %d/%d tasks", completed, total),
            true,
            progress,
            counts,
            errors,
        )
    end

    # Settled: a single merged summary line, no bar.
    workspace.cache.identity === nothing && return WorkspaceStatus(
        :none, "No Project", "Open a project folder to build a cache.",
        false, nothing, counts, errors)
    scan.state == :error && return WorkspaceStatus(
        :error, "Scan Error", scan.error, false, nothing, counts, errors)
    isempty(workspace.source_error) || return WorkspaceStatus(
        :error, "Source Error", workspace.source_error, false, nothing, counts, errors)
    workspace.cache_state == :error && return WorkspaceStatus(
        :error, "Cache Error", workspace.cache_error, false, nothing, counts, errors)
    workspace.cache_state == :missing && return WorkspaceStatus(
        :missing, "Cache Missing", "No cache has been built for this project yet.",
        false, nothing, counts, errors)
    scan.state == :canceled && return WorkspaceStatus(
        :none, "Canceled", "The last scan was canceled.", false, nothing, counts, errors)

    status = workspace.cache.status
    status isa ProjectCacheStatus || return WorkspaceStatus(
        :none, "Idle", "No cache loaded.", false, nothing, counts, errors)
    if !isempty(errors) || status.error_source_items > 0
        return WorkspaceStatus(
            :error, "Errors", "Some work failed.", false, nothing, counts, errors)
    elseif !(workspace.index.source isa SourceScan)
        return WorkspaceStatus(:fresh, "Loaded",
            "Cache loaded; source not checked.", false, nothing, counts, errors)
    end
    # A settled successful scan has ingested every difference it found, so the cache mirrors the
    # source; the scan's stale/new/deleted counts describe the work it did, not pending work.
    return WorkspaceStatus(:fresh, "Fresh", "Cache is current.", false, nothing, counts, errors)
end

function workspace_stage_counts(workspace::Workspace)::WorkspaceStageCounts
    cache_counts = cache_stage_summary(workspace.cache.db)
    status = workspace.cache.status
    found = status isa ProjectCacheStatus ? status.total_source_items :
        cache_counts.cached_sources + cache_counts.failed_interpret
    workspace.scan.state === :discovering && (found = max(found, workspace.scan.discovered[]))
    done = cache_counts.cached_sources + cache_counts.failed_interpret
    return WorkspaceStageCounts(
        source_item_noun(workspace.source),
        found,
        max(found - done, 0),
        cache_counts,
    )
end
