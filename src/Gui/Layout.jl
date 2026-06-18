using Printf
using Statistics: mean
import GLFW
import CImGui as ig
using NativeFileDialog: pick_folder

using ..Projects:
    DEFAULT_PROJECT,
    PROJECTS,
    project_description,
    project_name,
    scan_profile_summary,
    source_label
using ..ItemIndex: SourceScan
using ..Cache:
    ProjectCacheIdentity,
    ProjectCacheStatus
import ..Workspace
using ..Workspace:
    WorkspaceProgress,
    cache_work_running,
    cancel_cache!,
    cancel_scan!,
    poll_workspace!,
    scan_source!,
    source_scan_running,
    update_cache!

"""Describe the current cache operation for toolbar and progress rendering."""
function _cache_activity_model(state::BrowserState)::NamedTuple
    workspace = state.workspace
    job_state = workspace isa Workspace.Workspace ? workspace.cache_job.state : :idle
    progress = workspace isa Workspace.Workspace ?
        workspace.cache_job.progress :
        WorkspaceProgress()
    total = progress.total_source_items
    processed = progress.processed_source_items
    loaded = progress.loaded_items
    fraction = total > 0 ? Float32(clamp(processed / total, 0, 1)) : 0.0f0

    if job_state == :loading
        progress_text = total > 0 ?
            @sprintf("Read %d/%d cached source items, loaded %d items", processed, total, loaded) :
            @sprintf("Loaded %d items", loaded)
        return (
            title="Cache: Loading",
            detail="Reading cached items from the HDF5 file.",
            progress=progress_text,
            fraction,
            cancel_label="Cancel Reload",
        )
    elseif job_state == :writing
        phase = progress.phase
        progress_text = phase == :cache_finalize ?
            "Finalizing cache metadata and indexes" :
            total > 0 ?
            @sprintf("Processed %d/%d cache entries, cached %d items", processed, total, loaded) :
            @sprintf("Cached %d items", loaded)
        operation = workspace.cache.operation
        title = phase == :cache_finalize ? "Cache: Finalizing" :
            operation == :build ? "Cache: Building" :
            operation == :rebuild ? "Cache: Rebuilding" : "Cache: Updating"
        return (
            title,
            detail=phase == :cache_finalize ?
                "Finalizing the HDF5 cache." :
                "Writing source items into the HDF5 cache.",
            progress=progress_text,
            fraction,
            cancel_label="Cancel Cache Build",
        )
    elseif job_state == :canceling
        return (
            title="Canceling",
            detail="Waiting for the active background job to stop.",
            progress="Canceling...",
            fraction,
            cancel_label="Cancel",
        )
    elseif job_state == :canceled
        return (
            title="Canceled",
            detail="The last background job was canceled.",
            progress="Canceled",
            fraction,
            cancel_label="Cancel",
        )
    elseif job_state == :error
        return (
            title="Error",
            detail=workspace.cache_job.error,
            progress="Failed",
            fraction,
            cancel_label="Cancel",
        )
    elseif job_state == :done
        return (
            title="Ready",
            detail="No background cache or source job is running.",
            progress="Complete",
            fraction=1.0f0,
            cancel_label="Cancel",
        )
    end
    return (
        title="Idle",
        detail="No background cache or source job is running.",
        progress="Idle",
        fraction=0.0f0,
        cancel_label="Cancel",
    )
end

"""Describe source-scan progress independently from cache progress."""
function _source_rescan_progress_model(
    state::BrowserState,
)::Union{Nothing,NamedTuple}
    workspace = state.workspace
    workspace isa Workspace.Workspace || return nothing
    source_state = workspace.scan.state
    source_state in (:counting, :discovering, :scanning,
                     :canceling, :canceled, :done, :unchanged, :error) || return nothing

    progress = workspace.scan.progress
    total = progress.total_source_items
    processed = progress.processed_source_items
    loaded = progress.loaded_items
    skipped = progress.skipped_source_items
    fraction = total > 0 ? Float32(clamp(processed / total, 0, 1)) : 0.0f0

    # A single status line; the bar appears only while files are actively being read.
    if source_state == :error
        return (text=workspace.scan.error, fraction=0.0f0, show_bar=false)
    elseif source_state == :canceling
        return (text="Canceling scan...", fraction=fraction, show_bar=false)
    elseif source_state == :canceled
        return (text="Scan canceled", fraction=0.0f0, show_bar=false)
    elseif source_state in (:counting, :discovering)
        text = processed > 0 ? "Finding source items... $processed found" : "Finding source items..."
        return (text=text, fraction=0.0f0, show_bar=false)
    elseif source_state == :scanning
        text = total > 0 ?
            @sprintf("Reading %d/%d source items, %d items", processed, total, loaded) :
            "Reading source items..."
        return (text=text, fraction=fraction, show_bar=total > 0)
    elseif workspace.analysis.state == :analyzing
        analysis = workspace.analysis.progress
        total = analysis.total_source_items
        processed = analysis.processed_source_items
        loaded = analysis.loaded_items
        fraction = total > 0 ? Float32(clamp(processed / total, 0, 1)) : 0.0f0
        text = total > 0 ?
            @sprintf("Analyzing %d/%d, %d items", processed, total, loaded) :
            @sprintf("Analyzing %d items...", loaded)
        return (text=text, fraction=fraction, show_bar=total > 0)
    elseif workspace.analysis.state == :error
        return (text=workspace.analysis.error, fraction=0.0f0, show_bar=false)
    elseif source_state == :unchanged
        return (
            text=@sprintf("Up to date: %d items (%d source items unchanged)", loaded, total),
            fraction=0.0f0,
            show_bar=false,
        )
    else # :done
        text = skipped > 0 ?
            @sprintf("Scanned %d source items, %d items (%d skipped)", total, loaded, skipped) :
            @sprintf("Scanned %d source items, %d items", total, loaded)
        return (text=text, fraction=0.0f0, show_bar=false)
    end
end

"""Describe the cache button from the workspace cache state."""
function _cache_toolbar_model(state::BrowserState)::NamedTuple
    workspace = state.workspace
    identity = workspace isa Workspace.Workspace ? workspace.cache.identity : nothing
    status = workspace isa Workspace.Workspace ? workspace.cache.status : nothing
    cache_state = workspace isa Workspace.Workspace ?
        workspace.cache_job.state :
        :idle
    activity = _cache_activity_model(state)

    if identity === nothing
        return (
            label="Cache: No Project",
            color=(0.28, 0.28, 0.28, 1.0),
            detail="Open a project folder before building a cache.",
        )
    elseif cache_state == :writing
        return (
            label=activity.title,
            color=(0.18, 0.42, 0.78, 1.0),
            detail=activity.detail,
        )
    elseif cache_state == :loading
        return (
            label="Cache: Loading",
            color=(0.18, 0.42, 0.78, 1.0),
            detail="Reading cached items from the HDF5 file.",
        )
    elseif cache_state == :missing
        message = workspace.cache_job.error
        return (
            label="Cache: Missing",
            color=(0.82, 0.48, 0.12, 1.0),
            detail=isempty(message) ? "No HDF5 cache has been built for this project." : message,
        )
    elseif cache_state == :error
        return (
            label="Cache: Error",
            color=(0.72, 0.18, 0.18, 1.0),
            detail=workspace.cache_job.error,
        )
    elseif cache_state == :canceled
        return (
            label="Cache: Canceled",
            color=(0.45, 0.45, 0.45, 1.0),
            detail="Cache update was canceled.",
        )
    elseif status isa ProjectCacheStatus
        if !(workspace.index.source isa SourceScan)
            return (
                label="Cache: Loaded",
                color=(0.20, 0.58, 0.30, 1.0),
                detail="$(status.cached_source_items) source items cached; source not checked",
            )
        end
        if status.error_source_items > 0
            return (
                label="Cache: Errors",
                color=(0.72, 0.18, 0.18, 1.0),
                detail="$(status.error_source_items) cached source item(s) failed to transform",
            )
        end
        has_changes = status.stale_source_items > 0 ||
            status.new_source_items > 0 ||
            status.deleted_source_items > 0
        if has_changes
            return (
                label="Cache: Stale",
                color=(0.78, 0.58, 0.12, 1.0),
                detail="$(status.stale_source_items) stale, $(status.new_source_items) new, $(status.deleted_source_items) deleted",
            )
        end
        return (
            label="Cache: Fresh",
            color=(0.20, 0.58, 0.30, 1.0),
            detail="$(status.fresh_source_items) source items cached",
        )
    end

    return (
        label="Cache: Idle",
        color=(0.36, 0.36, 0.36, 1.0),
        detail="No cache loaded.",
    )
end

"""Render the cache status button and its detailed control popup."""
function _render_cache_toolbar_button!(state::BrowserState)::Nothing
    model = _cache_toolbar_model(state)
    color = model.color
    ig.PushStyleColor(ig.ImGuiCol_Button, model.color)
    ig.PushStyleColor(
        ig.ImGuiCol_ButtonHovered,
        (
            min(color[1] + 0.10, 1.0),
            min(color[2] + 0.10, 1.0),
            min(color[3] + 0.10, 1.0),
            color[4],
        ),
    )
    ig.PushStyleColor(
        ig.ImGuiCol_ButtonActive,
        (
            max(color[1] - 0.08, 0.0),
            max(color[2] - 0.08, 0.0),
            max(color[3] - 0.08, 0.0),
            color[4],
        ),
    )
    if ig.Button(model.label)
        ig.OpenPopup("cache_toolbar_popup")
    end
    ig.PopStyleColor()
    ig.PopStyleColor()
    ig.PopStyleColor()
    if ig.BeginItemTooltip()
        ig.TextUnformatted(model.detail)
        ig.EndTooltip()
    end
    ig.SetNextWindowSize((960, 560), ig.ImGuiCond_Always)
    if ig.BeginPopup("cache_toolbar_popup")
        _render_cache_controls!(state; compact=false)
        ig.EndPopup()
    end
    return nothing
end

"""Render diagnostic timing, allocation, memory, and OpenGL information."""
function render_perf_window(state::BrowserState)::Nothing
    state.show_performance_window || return nothing
    performance = state.performance

    if ig.Begin("Performance")
        raw_io = ig.GetIO()
        fps = unsafe_load(raw_io.Framerate)
        if fps > 0
            ig.Text(
                "FPS: $(round(fps; digits=1))  Frame: " *
                "$(round(1000 / fps; digits=2)) ms"
            )
        end

        if !isempty(performance.gl_info)
            gi = performance.gl_info
            for k in (:vendor, :renderer, :version)
                haskey(gi, k) && ig.Text("GL $(k): $(gi[k])")
            end
        end

        ig.Text("Tree nodes rendered: $(performance.node_count)")
        rows_visible = performance.item_rows_visible
        rows_rendered = performance.item_rows_rendered
        ig.Text("Items visible/rendered: $rows_visible / $rows_rendered")

        mem = _memory_snapshot()
        if mem.vmrss_kb !== nothing
            performance.memory_start_rss_kb === nothing &&
                (performance.memory_start_rss_kb = mem.vmrss_kb)
            start_rss = something(performance.memory_start_rss_kb)
            peak_rss = max(
                something(performance.memory_peak_rss_kb, mem.vmrss_kb),
                mem.vmrss_kb,
            )
            performance.memory_peak_rss_kb = peak_rss

            read_bytes = mem.read_bytes === nothing ? 0 : mem.read_bytes
            performance.memory_start_read_bytes === nothing &&
                (performance.memory_start_read_bytes = read_bytes)
            start_read = something(performance.memory_start_read_bytes)

            ig.Separator()
            ig.Text("Process memory")
            ig.BulletText(
                "RSS=$( _fmt_kb(mem.vmrss_kb) )  Δstart=$( _fmt_kb(mem.vmrss_kb - start_rss) )  peak=$( _fmt_kb(peak_rss) )",
            )
            ig.BulletText("Anon=$( _fmt_kb(mem.rssanon_kb) )  VSize=$( _fmt_kb(mem.vmsize_kb) )  VPeak=$( _fmt_kb(mem.vmpeak_kb) )")
            ig.BulletText("GC live=$( _fmt_bytes(mem.gc_live_bytes) )  maxrss=$( _fmt_bytes(mem.maxrss_bytes) )")
            ig.BulletText("read_bytes Δstart=$( _fmt_bytes(read_bytes - start_read) )")
        end

        timings = performance.timings
        allocs = performance.allocations

        plot_timing_keys = haskey(timings, :plot_draw) ? [:plot_draw] : Symbol[]
        other_timing_keys = sort([key for key in keys(timings) if key != :plot_draw])
        for (title, keys_to_show) in (("Plot drawing", plot_timing_keys), ("Frame/UI timings", other_timing_keys))
            isempty(keys_to_show) && continue
            ig.Separator()
            ig.Text(title)
            for key in keys_to_show
                samples = timings[key]
                isempty(samples) && continue
                bytes = get(allocs, key, Int[])
                last_ms = round(samples[end]; digits=2)
                mean_ms = round(mean(samples); digits=2)
                last_alloc = isempty(bytes) ? 0.0 : round(bytes[end] / 1024; digits=1)
                mean_alloc = isempty(bytes) ? 0.0 : round(mean(bytes) / 1024; digits=1)
                msg = @sprintf "%s: n=%d  last=%.2f ms  mean=%.2f ms  alloc_last=%.1f KB  alloc_mean=%.1f KB" String(key) length(samples) last_ms mean_ms last_alloc mean_alloc

                ig.BulletText(msg)
            end
        end

        workspace = state.workspace
        scan_rows = workspace isa Workspace.Workspace ?
            scan_profile_summary(workspace.project) :
            NamedTuple[]
        if ig.CollapsingHeader("Scan timings (per kind)", ig.ImGuiTreeNodeFlags_DefaultOpen)
            if isempty(scan_rows)
                ig.TextDisabled("No scan timing yet (cache hit or no scan run)")
            else
                total_read = sum(row.read_seconds for row in scan_rows)
                total_stats = sum(row.stats_seconds for row in scan_rows)
                source_items = sum(row.source_items for row in scan_rows)
                meas = sum(row.items for row in scan_rows)
                ig.Text(@sprintf(
                    "Last scan: %d source items, %d items", source_items, meas,
                ))
                # Reads and stats run across worker threads, so these are summed CPU-time, not
                # wall-clock; the total can exceed elapsed time by up to the thread count.
                ig.TextDisabled(@sprintf(
                    "Work summed across threads: read %.1fs, stats %.1fs", total_read, total_stats,
                ))
                for row in scan_rows
                    ig.BulletText(@sprintf(
                        "%s: %d source items read %.2fs  ·  %d items, stats %.2fs",
                        String(row.kind), row.source_items, row.read_seconds,
                        row.items, row.stats_seconds,
                    ))
                end
            end
        end

        if ig.Button("Clear timings")
            empty!(performance.timings)
            empty!(performance.allocations)
        end
    end
    ig.End()
    return nothing
end

"""Render project, cache, annotation, workflow, and debug controls."""
function render_menu_bar(state::BrowserState)::Nothing
    workspace = state.workspace
    if ig.BeginMenuBar()
        if ig.BeginMenu("Project")
            ig.TextDisabled(
                workspace isa Workspace.Workspace ?
                "Active: $(source_label(workspace.source))" :
                "No project loaded",
            )
            ig.Separator()

            if ig.MenuItem("Open Folder...")
                path = pick_folder()
                if !isnothing(path) && !isempty(path)
                    @info "Selected path: $path"
                    _open_project_path!(state, path)
                end
            end

            recents = state.recent_projects
            if ig.BeginMenu("Recent Projects")
                if isempty(recents)
                    ig.TextDisabled("No recent projects")
                else
                    for (idx, entry) in enumerate(recents)
                        path = entry.path
                        isempty(path) && continue
                        pref = entry.project_preference
                        label = "$(basename(path)) [$pref]###recent_project_$idx"
                        if ig.MenuItem(label)
                            _open_project_path!(state, path)
                        end
                        if ig.BeginItemTooltip()
                            ig.TextUnformatted(path)
                            ig.EndTooltip()
                        end
                    end
                end
                ig.EndMenu()
            end

            has_root = workspace isa Workspace.Workspace
            source_rescan_running =
                has_root && source_scan_running(workspace)
            rescan_label = source_rescan_running ? "Cancel Rescan" : "Rescan"
            can_rescan = has_root
            !can_rescan && ig.BeginDisabled()
            if ig.MenuItem(rescan_label)
                if source_rescan_running
                    cancel_scan!(workspace)
                else
                    @info "Rescanning source: $(source_label(workspace.source))"
                    scan_source!(workspace)
                end
            end
            !can_rescan && ig.EndDisabled()

            ig.Separator()
            if ig.MenuItem("Project Settings", C_NULL, state.show_project_window)
                state.show_project_window = !state.show_project_window
            end
            ig.EndMenu()
        end
        _render_cache_toolbar_button!(state)
        show_bad_effective = _show_bad_effective(state)
        bad_visibility_toggle_enabled = _tag_state_ready(state) ||
                                        isempty(state.tag_state_error)
        !bad_visibility_toggle_enabled && ig.BeginDisabled()
        if ig.MenuItem("Show Bad", C_NULL, show_bad_effective)
            state.show_bad = !state.show_bad
        end
        if !bad_visibility_toggle_enabled
            ig.EndDisabled()
            if ig.BeginItemTooltip()
                ig.TextUnformatted("Fix tags.txt and rescan before hiding bad items")
                ig.EndTooltip()
            end
        end
        _render_table_inspector_menu!(state)
        if ig.BeginMenu("Debug")
            if ig.MenuItem(
                "Performance Window",
                C_NULL,
                state.show_performance_window,
            )
                state.show_performance_window = !state.show_performance_window
            end
            plots = state.plots
            if ig.MenuItem("Debug Plot Mode", C_NULL, plots.debug)
                plots.debug = !plots.debug
                clear_plot_views!(plots)
            end
            ig.EndMenu()
        end
        ig.EndMenuBar()
    end
    return nothing
end

"""
Draw source and cache progress, differences, errors, and controls.
"""
function _render_cache_controls!(
    state::BrowserState;
    compact::Bool,
)::Nothing
    workspace = state.workspace
    identity = workspace isa Workspace.Workspace ? workspace.cache.identity : nothing
    status = workspace isa Workspace.Workspace ? workspace.cache.status : nothing
    cache_state = workspace isa Workspace.Workspace ?
        workspace.cache_job.state :
        :idle
    model = _cache_toolbar_model(state)
    running = workspace isa Workspace.Workspace &&
        cache_work_running(workspace)
    activity = _cache_activity_model(state)

    if compact
        ig.Text("Cache")
    else
        ig.TextColored(model.color, model.label)
    end
    ig.TextWrapped(model.detail)
    if running
        if cache_state in (:loading, :writing, :canceling)
            ig.TextDisabled(activity.title)
            ig.TextDisabled(activity.progress)
            ig.ProgressBar(activity.fraction, (-1, 0))
        end
    end

    source_progress = _source_rescan_progress_model(state)
    if source_progress !== nothing
        ig.Separator()
        ig.TextDisabled(source_progress.text)
        source_progress.show_bar &&
            ig.ProgressBar(source_progress.fraction, (-1, 0))
    end

    if identity isa ProjectCacheIdentity
        ig.Separator()
        ig.Text("Source: $(identity.source_label)")
        ig.TextWrapped("Source ID: $(identity.source_id)")
        ig.TextWrapped("File: $(identity.cache_path)")
    end

    if status isa ProjectCacheStatus
        ig.Separator()
        ig.Text("Source Items")
        ig.Text("Total: $(status.total_source_items)")
        ig.Text("Cached: $(status.cached_source_items)")
        ig.Text("Fresh: $(status.fresh_source_items)")
        ig.Text("Stale: $(status.stale_source_items)")
        ig.Text("New: $(status.new_source_items)")
        ig.Text("Deleted: $(status.deleted_source_items)")
        status.error_source_items > 0 && ig.TextColored(
            (1.0, 0.35, 0.35, 1.0),
            "Errors: $(status.error_source_items)",
        )
    elseif cache_state == :error
        message = workspace.cache_job.error
        !isempty(message) && ig.TextWrapped(message)
    end

    cache_errors = workspace isa Workspace.Workspace ?
        sort!(collect(workspace.index.analysis_errors); by=first) :
        Pair{String,String}[]
    if !isempty(cache_errors)
        ig.Separator()
        ig.TextColored((1.0, 0.35, 0.35, 1.0), "Source Item Errors")
        ig.TextDisabled("Select a source item to show its items")
        for (index, source_item_error) in enumerate(cache_errors)
            index > 20 && break
            source_item_id, message = source_item_error
            ig.PushID(source_item_id)
            if ig.TextLink(basename(source_item_id)) && select_source_item!(state, source_item_id)
                state.tree_filter = ""
                state.item_filter = ""
                state.reset_project_filters = true
            end
            if ig.BeginItemTooltip()
                ig.TextUnformatted(source_item_id)
                ig.EndTooltip()
            end
            ig.TextWrapped(first(split(message, '\n'; limit=2)))
            ig.PopID()
            index < min(length(cache_errors), 20) && ig.Separator()
        end
        length(cache_errors) > 20 && ig.TextDisabled("$(length(cache_errors) - 20) more file errors")
    end

    ig.Separator()
    has_identity = identity isa ProjectCacheIdentity
    can_update = has_identity &&
        status isa ProjectCacheStatus &&
        workspace.index.source isa SourceScan &&
        cache_state != :missing &&
        cache_state != :error &&
        (status.stale_source_items > 0 ||
         status.new_source_items > 0 ||
         status.deleted_source_items > 0)
    can_build = has_identity &&
        workspace.index.source isa SourceScan &&
        cache_state != :error
    can_scan = workspace isa Workspace.Workspace
    scan_running = can_scan && source_scan_running(workspace)

    !can_scan && ig.BeginDisabled()
    scan_label = scan_running ? "Cancel Scan" : "Scan Source"
    if ig.Button(scan_label, (-1, 0))
        if scan_running
            cancel_scan!(workspace)
        else
            scan_source!(workspace)
        end
    end
    !can_scan && ig.EndDisabled()

    (!can_update || running) && ig.BeginDisabled()
    if ig.Button("Update Cache", (-1, 0))
        update_cache!(workspace)
    end
    (!can_update || running) && ig.EndDisabled()

    (!can_build || running) && ig.BeginDisabled()
    build_label = cache_state == :missing ? "Build Cache" : "Rebuild Cache"
    if ig.Button(build_label, (-1, 0))
        update_cache!(workspace; rebuild=true)
    end
    (!can_build || running) && ig.EndDisabled()

    if running
        if ig.Button(activity.cancel_label, (-1, 0))
            cancel_cache!(workspace)
        end
    end
    return nothing
end

"""Render project selection and project-specific settings."""
function render_project_window(state::BrowserState)::Nothing
    state.show_project_window || return nothing

    open_ref = Ref(true)
    if ig.Begin("Project Settings###project_window", open_ref, ig.ImGuiWindowFlags_AlwaysAutoResize)
        workspace = state.workspace
        if workspace isa Workspace.Workspace
            ig.Text("Active: $(source_label(workspace.source))")
        else
            ig.TextDisabled("No folder loaded yet")
        end

        ig.Separator()

        if state.project_locked
            ig.TextDisabled("Project supplied by the Julia caller")
        else
            pref = state.project_preference
            changed = false

            default_project = something(DEFAULT_PROJECT[])
            default_label = "Default ($(project_name(default_project)))"

            if ig.RadioButton(default_label, pref == "auto")
                state.project_preference = "auto"
                changed = true
            end
            ig.SameLine()
            _helpmarker("Use the default project without trying to infer the project from the files in the folder.")

            for p in PROJECTS
                pn = project_name(p)
                if ig.RadioButton(pn, pref == pn)
                    state.project_preference = pn
                    changed = true
                end
                ig.SameLine()
                _helpmarker(project_description(p))
            end

            if changed
                current_root = (
                    workspace isa Workspace.Workspace &&
                    hasproperty(workspace.source, :root_path)
                ) ? workspace.source.root_path : ""
                _persist_preferences!(
                    state;
                    path=isempty(current_root) ? nothing : current_root,
                )
                if workspace isa Workspace.Workspace
                    @info(
                        "Project preference changed to '$(state.project_preference)' - " *
                        "reloading cache",
                    )
                    _open_project_path!(state, current_root)
                end
            end
        end
    end
    open_ref[] || (state.show_project_window = false)
    ig.End()
    return nothing
end

"""Create the initial hierarchy, information, and plot docking layout."""
function _setup_docking_layout!(dockspace_id::UInt32)::Nothing
    vp = ig.GetMainViewport()
    sz = unsafe_load(vp.Size)

    ig.DockBuilderRemoveNode(dockspace_id)
    ig.DockBuilderAddNode(dockspace_id, Int(ig.ImGuiDockNodeFlags_DockSpace))
    ig.DockBuilderSetNodeSize(dockspace_id, sz)

    # Vertical split: left column 2/5, right column 3/5
    left  = Ref{UInt32}(0)
    right = Ref{UInt32}(0)
    ig.DockBuilderSplitNode(dockspace_id, ig.ImGuiDir_Left, 2/5, left, right)

    # Left column: top 3/4 hierarchy, bottom 1/4 info
    top_left    = Ref{UInt32}(0)
    bottom_left = Ref{UInt32}(0)
    ig.DockBuilderSplitNode(left[], ig.ImGuiDir_Up, 3/4, top_left, bottom_left)

    ig.DockBuilderDockWindow("Hierarchy",         top_left[])
    ig.DockBuilderDockWindow("Plot Area",          right[])
    ig.DockBuilderDockWindow("Information Panel",  bottom_left[])

    ig.DockBuilderFinish(dockspace_id)
    return nothing
end

"""Return whether the native GLFW window requested shutdown."""
function _window_close_requested()::Bool
    window = ig.current_window()
    window === nothing && return false
    return GLFW.WindowShouldClose(window)
end

"""Create the ImGui context with the docking/viewport configuration the browser needs."""
function _init_browser_context!()
    ig.set_backend(:GlfwOpenGL3)
    ctx = ig.CreateContext()
    io = ig.GetIO()
    io.ConfigFlags = unsafe_load(io.ConfigFlags) | ig.ImGuiConfigFlags_DockingEnable
    io.ConfigFlags = unsafe_load(io.ConfigFlags) | ig.ImGuiConfigFlags_ViewportsEnable
    io.ConfigFlags = unsafe_load(io.ConfigFlags) | ig.ImGuiConfigFlags_NavEnableKeyboard
    io.IniFilename = Ptr{Cchar}(C_NULL)   # disable layout persistence
    ig.StyleColorsDark()
    return ctx
end

"""
On macOS, promote the process to a regular foreground app so its window gets a Dock icon and shows
up in Cmd-Tab. When the renderloop is driven from the REPL the process can otherwise stay an
accessory app with no Dock presence. No-op off macOS; never throws, since a failure here must not
take down the render loop.
"""
function _promote_to_foreground_app()
    Sys.isapple() || return nothing
    try
        objc_class(name) = @ccall objc_getClass(name::Cstring)::Ptr{Cvoid}
        objc_sel(name) = @ccall sel_registerName(name::Cstring)::Ptr{Cvoid}
        nsapplication = objc_class("NSApplication")
        nsapplication == C_NULL && return nothing
        app = @ccall objc_msgSend(
            nsapplication::Ptr{Cvoid}, objc_sel("sharedApplication")::Ptr{Cvoid},
        )::Ptr{Cvoid}
        app == C_NULL && return nothing
        # NSApplicationActivationPolicyRegular = 0
        @ccall objc_msgSend(
            app::Ptr{Cvoid}, objc_sel("setActivationPolicy:")::Ptr{Cvoid}, 0::Clong,
        )::Bool
        @ccall objc_msgSend(
            app::Ptr{Cvoid}, objc_sel("activateIgnoringOtherApps:")::Ptr{Cvoid}, true::Bool,
        )::Cvoid
    catch error
        @warn "Could not promote to a macOS foreground app" exception=error
    end
    return nothing
end

"""
Run the browser render loop for a prepared state. With `wait=false` the loop runs as a background
task pinned to thread 1 (required for GLFW) and the call returns that task, leaving the REPL live.
"""
function _run_browser(state::BrowserState, ctx; engine, spawn, wait::Bool)
    first_frame   = Ref(true)
    setup_layout  = Ref(true)
    return ig.render(
        ctx;
        engine,
        spawn,
        wait,
        window_size=(1920, 1080),
        window_title="Measurement Browser",
        opengl_version=v"3.3",
        wait_events=false,
        on_exit=() -> begin
            _shutdown_background_jobs!(state)
            _save_project_view_if_changed!(state)
            _print_perf_summary(state)
        end,
    ) do
        if _window_close_requested()
            _shutdown_background_jobs!(state)
            return :imgui_exit_loop
        end
        state.performance.frame += 1
        workspace = state.workspace
        workspace isa Workspace.Workspace && poll_workspace!(workspace)
        if first_frame[]
            state.performance.gl_info = _gl_info()
            _promote_to_foreground_app()
            first_frame[] = false
        end
        dockspace_id = ig.DockSpaceOverViewport(0, ig.GetMainViewport())
        if setup_layout[]
            setup_layout[] = false
            _setup_docking_layout!(dockspace_id)
        end
        if !(workspace isa Workspace.Workspace && source_scan_running(workspace))
            _time!(state, :plot_warmup) do
                ensure_plot_runtime_warmed!(state)
            end
        end
        render_selection_window(state)
        render_project_window(state)
        _time!(state, :info) do
            render_info_window(state)
        end
        _time!(state, :table_inspector) do
            render_table_inspector_window(state)
        end
        _time!(state, :plot) do
            render_plot_window(state)
        end
        _time!(state, :extra_plots) do
            render_additional_plot_windows(state)
        end
        _time!(state, :perf_window) do
            render_perf_window(state)
        end
        _save_project_view_if_changed!(state)
        # Show metadata guidance modal if needed
        render_collection_metadata_modal(state)
    end
end

"""
Open the browser on an already-opened workspace (see `open_workspace`). With `wait=false` the REPL
stays interactive while the window runs; the returned task completes when the window closes.
"""
function open_browser(
    workspace::Workspace.Workspace;
    engine::Any=nothing,
    spawn::Int=1,
    wait::Bool=true,
)
    prefs = _load_prefs()
    state = BrowserState(
        project_locked=true,
        project_preference=project_name(workspace.project),
        recent_projects=_parse_recent_projects(prefs),
    )
    ctx = _init_browser_context!()
    _attach_workspace!(state, workspace)
    return _run_browser(state, ctx; engine, spawn, wait)
end
