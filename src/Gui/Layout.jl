function _cache_activity_model(ui_state)
    state = get(ui_state, :cache_state, :idle)
    progress = get(ui_state, :cache_progress, _new_scan_progress())
    total = get(progress, :total_csv, 0)
    processed = get(progress, :processed_csv, 0)
    loaded = get(progress, :loaded_measurements, 0)
    fraction = total > 0 ? Float32(clamp(processed / total, 0, 1)) : 0.0f0

    if state == :loading
        progress_text = total > 0 ?
            @sprintf("Read %d/%d cached files, loaded %d measurements", processed, total, loaded) :
            @sprintf("Loaded %d measurements", loaded)
        return (
            title="Cache: Loading",
            detail="Reading cached measurements from the HDF5 file.",
            progress=progress_text,
            fraction,
            cancel_label="Cancel Reload",
        )
    elseif state == :writing
        phase = get(progress, :phase, :idle)
        progress_text = phase == :cache_finalize ?
            "Finalizing cache metadata and indexes" :
            total > 0 ?
            @sprintf("Processed %d/%d cache entries, cached %d measurements", processed, total, loaded) :
            @sprintf("Cached %d measurements", loaded)
        operation = get(ui_state, :cache_operation, :update)
        title = phase == :cache_finalize ? "Cache: Finalizing" :
            operation == :build ? "Cache: Building" :
            operation == :rebuild ? "Cache: Rebuilding" : "Cache: Updating"
        return (
            title,
            detail=phase == :cache_finalize ?
                "Finalizing the HDF5 cache." :
                "Writing source measurements into the HDF5 cache.",
            progress=progress_text,
            fraction,
            cancel_label="Cancel Cache Build",
        )
    elseif state == :canceling
        return (
            title="Canceling",
            detail="Waiting for the active background job to stop.",
            progress="Canceling...",
            fraction,
            cancel_label="Cancel",
        )
    elseif state == :canceled
        return (
            title="Canceled",
            detail="The last background job was canceled.",
            progress="Canceled",
            fraction,
            cancel_label="Cancel",
        )
    elseif state == :error
        return (
            title="Error",
            detail=get(ui_state, :cache_error, "Cache job failed."),
            progress="Failed",
            fraction,
            cancel_label="Cancel",
        )
    elseif state == :done
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

function _progress_fraction(progress)
    total = get(progress, :total_csv, 0)
    processed = get(progress, :processed_csv, 0)
    total <= 0 && return 0.0f0
    return Float32(clamp(processed / total, 0, 1))
end

function _render_progress_indicator!(ui_state, item)
    get(item, :show_bar, true) || return
    ig.ProgressBar(item.fraction, (-1, 0))
end

function _source_rescan_progress_model(ui_state)
    source_state = get(ui_state, :source_scan_state, :idle)
    if !(source_state in (:counting, :discovering, :scanning, :canceling, :done, :canceled, :error))
        return nothing
    end
    if source_state == :error
        return (
            title="Rescan: Error",
            progress=get(ui_state, :source_scan_error, "Source scan failed."),
            fraction=0.0f0,
            show_bar=false,
        )
    end

    progress = get(ui_state, :source_scan_progress, _new_scan_progress())
    total = get(progress, :total_csv, 0)
    processed = get(progress, :processed_csv, 0)
    loaded = get(progress, :loaded_measurements, 0)
    skipped = get(progress, :skipped_csv, 0)
    text = if source_state == :canceling
        "Canceling source scan..."
    elseif source_state == :canceled
        "Source scan canceled"
    elseif source_state == :done
        total > 0 ?
            "Source scan complete: scanned $processed/$total CSV files" :
            "Source scan complete"
    elseif source_state in (:counting, :discovering)
        processed > 0 ? "Finding source CSV files: $processed found" : "Finding source CSV files"
    else
        total > 0 ?
            @sprintf("Scanning %d/%d source files, loaded %d measurements, skipped %d", processed, total, loaded, skipped) :
            "Scanning source files"
    end
    title = if source_state == :done
        "Source: Complete"
    elseif source_state == :canceled
        "Source: Canceled"
    elseif source_state == :canceling
        "Source: Canceling"
    elseif source_state in (:counting, :discovering)
        "Source: Finding Files"
    else
        "Source: Scanning"
    end
    return (
        title,
        progress=text,
        fraction=source_state == :done ? 1.0f0 : _progress_fraction(progress),
        show_bar=total > 0 || source_state == :done,
    )
end

function _cache_progress_models(ui_state)
    models = NamedTuple[]
    cache_state = get(ui_state, :cache_state, :idle)
    if cache_state in (:loading, :writing, :canceling)
        activity = _cache_activity_model(ui_state)
        push!(models, (
            title=activity.title,
            progress=activity.progress,
            fraction=activity.fraction,
            show_bar=true,
        ))
    end
    return models
end

function _source_progress_models(ui_state)
    models = NamedTuple[]
    rescan_model = _source_rescan_progress_model(ui_state)
    if rescan_model !== nothing
        push!(models, rescan_model)
    end

    return models
end

function _cache_toolbar_model(ui_state)
    identity = get(ui_state, :cache_identity, nothing)
    status = get(ui_state, :cache_status, nothing)
    cache_state = get(ui_state, :cache_state, :idle)
    activity = _cache_activity_model(ui_state)

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
            detail="Reading cached measurements from the HDF5 file.",
        )
    elseif cache_state == :missing
        return (
            label="Cache: Missing",
            color=(0.82, 0.48, 0.12, 1.0),
            detail="No HDF5 cache has been built for this project.",
        )
    elseif cache_state == :error
        return (
            label="Cache: Error",
            color=(0.72, 0.18, 0.18, 1.0),
            detail=get(ui_state, :cache_error, "Cache error"),
        )
    elseif cache_state == :canceled
        return (
            label="Cache: Canceled",
            color=(0.45, 0.45, 0.45, 1.0),
            detail="Cache update was canceled.",
        )
    elseif status isa ProjectCacheStatus
        if !get(ui_state, :cache_source_checked, true)
            return (
                label="Cache: Loaded",
                color=(0.20, 0.58, 0.30, 1.0),
                detail="$(status.cached_files) files cached; source not checked",
            )
        end
        if status.error_files > 0
            return (
                label="Cache: Errors",
                color=(0.72, 0.18, 0.18, 1.0),
                detail="$(status.error_files) cached file(s) failed to transform",
            )
        end
        has_changes = status.stale_files > 0 || status.new_files > 0 || status.deleted_files > 0
        if has_changes
            return (
                label="Cache: Stale",
                color=(0.78, 0.58, 0.12, 1.0),
                detail="$(status.stale_files) stale, $(status.new_files) new, $(status.deleted_files) deleted",
            )
        end
        return (
            label="Cache: Fresh",
            color=(0.20, 0.58, 0.30, 1.0),
            detail="$(status.fresh_files) files cached",
        )
    end

    return (
        label="Cache: Idle",
        color=(0.36, 0.36, 0.36, 1.0),
        detail="No cache loaded.",
    )
end

function _source_scan_matches_cache(ui_state, identity::ProjectCacheIdentity)
    source = get(ui_state, :source_scan_result, nothing)
    return source isa SourceScan &&
        source.root_path == identity.root_path &&
        project_name(source.project) == identity.project_name
end

function _cache_build_label(ui_state, identity::ProjectCacheIdentity, cache_state::Symbol)
    has_source = _source_scan_matches_cache(ui_state, identity)
    verb = cache_state == :missing ? "Build" : "Rebuild"
    return has_source ? "$(verb) Cache" : "Scan and $(verb) Cache"
end

function _cache_button_hover_color(color)
    return (
        min(color[1] + 0.10, 1.0),
        min(color[2] + 0.10, 1.0),
        min(color[3] + 0.10, 1.0),
        color[4],
    )
end

function _cache_button_active_color(color)
    return (
        max(color[1] - 0.08, 0.0),
        max(color[2] - 0.08, 0.0),
        max(color[3] - 0.08, 0.0),
        color[4],
    )
end

function _render_cache_toolbar_button!(ui_state)
    model = _cache_toolbar_model(ui_state)
    ig.PushStyleColor(ig.ImGuiCol_Button, model.color)
    ig.PushStyleColor(ig.ImGuiCol_ButtonHovered, _cache_button_hover_color(model.color))
    ig.PushStyleColor(ig.ImGuiCol_ButtonActive, _cache_button_active_color(model.color))
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
    _render_cache_toolbar_popup!(ui_state)
end

function _render_cache_toolbar_popup!(ui_state)
    ig.SetNextWindowSize((960, 560), ig.ImGuiCond_Always)
    if ig.BeginPopup("cache_toolbar_popup")
        _render_cache_controls!(ui_state; compact=false)
        ig.EndPopup()
    end
end

function render_perf_window(ui_state)
    if !get(ui_state, :show_performance_window, false)
        return
    end

    if ig.Begin("Performance")
        raw_io = ig.GetIO()
        fps = unsafe_load(raw_io.Framerate)
        if fps > 0
            ig.Text(
                "FPS: $(round(fps; digits=1))  Frame: " *
                "$(round(1000 / fps; digits=2)) ms"
            )
        end

        if haskey(ui_state, :_gl_info)
            gi = ui_state[:_gl_info]
            for k in (:vendor, :renderer, :version)
                haskey(gi, k) && ig.Text("GL $(k): $(gi[k])")
            end
        end

        if haskey(ui_state, :_node_count)
            ig.Text("Tree nodes rendered: $(ui_state[:_node_count])")
        end
        rows_visible = get(ui_state, :_measurement_rows_visible, 0)
        rows_rendered = get(ui_state, :_measurement_rows_rendered, 0)
        ig.Text("Measurements visible/rendered: $rows_visible / $rows_rendered")

        mem = _memory_snapshot()
        if mem.vmrss_kb !== nothing
            start_rss = get!(ui_state, :_mem_start_rss_kb, mem.vmrss_kb)
            peak_rss = max(get(ui_state, :_mem_peak_rss_kb, mem.vmrss_kb), mem.vmrss_kb)
            ui_state[:_mem_peak_rss_kb] = peak_rss

            read_bytes = mem.read_bytes === nothing ? 0 : mem.read_bytes
            start_read = get!(ui_state, :_mem_start_read_bytes, read_bytes)

            ig.Separator()
            ig.Text("Process memory")
            ig.BulletText(
                "RSS=$( _fmt_kb(mem.vmrss_kb) )  Δstart=$( _fmt_kb(mem.vmrss_kb - start_rss) )  peak=$( _fmt_kb(peak_rss) )",
            )
            ig.BulletText("Anon=$( _fmt_kb(mem.rssanon_kb) )  VSize=$( _fmt_kb(mem.vmsize_kb) )  VPeak=$( _fmt_kb(mem.vmpeak_kb) )")
            ig.BulletText("GC live=$( _fmt_bytes(mem.gc_live_bytes) )  maxrss=$( _fmt_bytes(mem.maxrss_bytes) )")
            ig.BulletText("read_bytes Δstart=$( _fmt_bytes(read_bytes - start_read) )")
        end

        timings = get(ui_state, :_timings,
            Dict{Symbol,Vector{Float64}}())
        allocs = get(ui_state, :_allocs,
            Dict{Symbol,Vector{Int}}())

        for (k, v) in timings
            isempty(v) && continue
            a = get(allocs, k, Int[])
            last_ms = round(v[end]; digits=2)
            mean_ms = round(mean(v); digits=2)
            last_alloc = isempty(a) ? 0.0 : round(a[end] / 1024; digits=1)
            mean_alloc = isempty(a) ? 0.0 : round(mean(a) / 1024; digits=1)
            msg = @sprintf "%s: last=%.2f ms  mean=%.2f ms  alloc_last=%.1f KB  alloc_mean=%.1f KB" String(k) last_ms mean_ms last_alloc mean_alloc

            ig.BulletText(msg)
        end

        profiles = get(ui_state, :figure_script_job_profiles, Any[])
        if ig.CollapsingHeader("Figure-Script Jobs", ig.ImGuiTreeNodeFlags_DefaultOpen)
            if isempty(profiles)
                ig.TextDisabled("No figure-script jobs profiled yet")
            else
                latest = profiles[end]
                latest_profile = latest.profile
                ig.Text(
                    "Last: $(latest.operation) $(repr(latest_profile.group_name))  " *
                    "$(round(latest_profile.total_ms; digits=1)) ms  " *
                    "selected=$(latest_profile.selected_count) / total=$(latest_profile.measurement_count)",
                )
                if !isempty(latest_profile.sections)
                    ig.Text("Sections")
                    for section in latest_profile.sections
                        ig.BulletText(
                            "$(section.key): calls=$(section.calls)  " *
                            "time=$(round(section.duration_ms; digits=2)) ms  " *
                            "alloc=$(round(section.alloc_bytes / 1024; digits=1)) KB",
                        )
                    end
                end
                if !isempty(latest_profile.counters)
                    ig.Text("Counters")
                    for (key, value) in sort!(collect(latest_profile.counters); by=entry -> String(first(entry)))
                        ig.BulletText("$(key): $(value)")
                    end
                end

                if length(profiles) > 1
                    ig.Text("History")
                    for entry in Iterators.reverse(profiles[max(1, end - 5):end-1])
                        profile = entry.profile
                        ig.BulletText(
                            "$(entry.operation) $(repr(profile.group_name)): " *
                            "$(round(profile.total_ms; digits=1)) ms  " *
                            "selected=$(profile.selected_count)",
                        )
                    end
                end
            end
        end

        if ig.Button("Clear timings")
            empty!(get!(() -> Dict{Symbol,Vector{Float64}}(), ui_state, :_timings))
            empty!(get!(() -> Dict{Symbol,Vector{Int}}(), ui_state, :_allocs))
            empty!(get!(ui_state, :figure_script_job_profiles) do
                Any[]
            end)
        end
    end
    ig.End()
end

function render_menu_bar(ui_state)
    if ig.BeginMenuBar()
        if ig.BeginMenu("Project")
            ig.TextDisabled(_project_status_text(ui_state))
            ig.Separator()

            if ig.MenuItem("Open Folder...")
                path = pick_folder()
                if !isnothing(path) && !isempty(path)
                    @info "Selected path: $path"
                    _open_project_path!(ui_state, path)
                end
            end

            recents = get(ui_state, :recent_projects, Dict{String,String}[])
            if ig.BeginMenu("Recent Projects")
                if isempty(recents)
                    ig.TextDisabled("No recent projects")
                else
                    for (idx, entry) in enumerate(recents)
                        path = get(entry, "path", "")
                        isempty(path) && continue
                        pref = get(entry, "project_preference", "auto")
                        label = "$(basename(path)) [$pref]###recent_project_$idx"
                        if ig.MenuItem(label)
                            _open_project_path!(ui_state, path)
                        end
                        if ig.BeginItemTooltip()
                            ig.TextUnformatted(path)
                            ig.EndTooltip()
                        end
                    end
                end
                ig.EndMenu()
            end

            has_root = haskey(ui_state, :root_path) && !isempty(ui_state[:root_path])
            source_rescan_running = _source_scan_running(ui_state)
            rescan_label = source_rescan_running ? "Cancel Rescan" : "Rescan"
            can_rescan = has_root
            !can_rescan && ig.BeginDisabled()
            if ig.MenuItem(rescan_label)
                if source_rescan_running
                    _request_source_scan_cancel!(ui_state)
                else
                    proj = haskey(ui_state, :project) ?
                        ui_state[:project] :
                        _project_for_preference(get(ui_state, :project_preference, "auto"))
                    @info "Rescanning path: $(ui_state[:root_path])"
                    _launch_source_scan_job!(ui_state, ui_state[:root_path], proj)
                end
            end
            !can_rescan && ig.EndDisabled()

            ig.Separator()
            if ig.MenuItem("Project Settings", C_NULL, get(ui_state, :show_project_window, false))
                ui_state[:show_project_window] = !get(ui_state, :show_project_window, false)
            end
            ig.EndMenu()
        end
        _render_cache_toolbar_button!(ui_state)
        show_bad_effective = _show_bad_effective(ui_state)
        bad_visibility_toggle_enabled = _tag_state_ready(ui_state) ||
                                        isempty(get(ui_state, :tag_state_error, ""))
        !bad_visibility_toggle_enabled && ig.BeginDisabled()
        if ig.MenuItem("Show Bad", C_NULL, show_bad_effective)
            ui_state[:show_bad] = !get(ui_state, :show_bad, true)
            _apply_visible_selection!(ui_state)
        end
        if !bad_visibility_toggle_enabled
            ig.EndDisabled()
            if ig.BeginItemTooltip()
                ig.TextUnformatted("Fix tags.txt and rescan before hiding bad items")
                ig.EndTooltip()
            end
        end
        has_measurement_selection = !isempty(get(ui_state, :selected_measurements, MeasurementInfo[]))
        can_open_figure_scripts = has_measurement_selection ||
                                  !isempty(_figure_script_groups(ui_state)) ||
                                  get(ui_state, :show_figure_script_window, false)
        !can_open_figure_scripts && ig.BeginDisabled()
        if ig.MenuItem("Figure Script", C_NULL, get(ui_state, :show_figure_script_window, false))
            ui_state[:show_figure_script_window] = !get(ui_state, :show_figure_script_window, false)
            ui_state[:figure_script_root_path] = get(ui_state, :root_path, "")
        end
        if !can_open_figure_scripts
            ig.EndDisabled()
            if ig.BeginItemTooltip()
                ig.TextUnformatted("Select one or more measurements to build a figure script")
                ig.EndTooltip()
            end
        end
        if ig.BeginMenu("Debug")
            if ig.MenuItem("Performance Window", C_NULL, get(ui_state, :show_performance_window, false))
                ui_state[:show_performance_window] = !get(ui_state, :show_performance_window, false)
            end
            if ig.MenuItem("Debug Plot Mode", C_NULL, get(ui_state, :debug_plot_mode, false))
                ui_state[:debug_plot_mode] = !get(ui_state, :debug_plot_mode, false)
                # Invalidate cached figures when toggled
                _clear_plot_jobs!(ui_state)
                delete!(ui_state, :plot_figure)
                delete!(ui_state, :_last_plot_key)
            end
            ig.EndMenu()
        end
        ig.EndMenuBar()
    end
end

function _render_cache_controls!(ui_state; compact::Bool)
    identity = get(ui_state, :cache_identity, nothing)
    status = get(ui_state, :cache_status, nothing)
    cache_state = get(ui_state, :cache_state, :idle)
    model = _cache_toolbar_model(ui_state)
    running = _cache_action_blocked(ui_state)
    activity = _cache_activity_model(ui_state)

    if compact
        ig.Text("Cache")
    else
        ig.TextColored(model.color, model.label)
    end
    ig.TextWrapped(model.detail)
    if running
        cache_progress = _cache_progress_models(ui_state)
        if !isempty(cache_progress)
            for item in cache_progress
                ig.TextDisabled(item.title)
                ig.TextDisabled(item.progress)
                _render_progress_indicator!(ui_state, item)
            end
        end
    end

    source_progress = _source_progress_models(ui_state)
    if !isempty(source_progress)
        for item in source_progress
            ig.Separator()
            ig.TextDisabled(item.title)
            ig.TextDisabled(item.progress)
            _render_progress_indicator!(ui_state, item)
        end
    end

    if identity isa ProjectCacheIdentity
        ig.Separator()
        ig.Text("Project: $(identity.project_name)")
        ig.TextWrapped("Root: $(identity.root_path)")
        ig.TextWrapped("File: $(identity.cache_path)")
    end

    if status isa ProjectCacheStatus
        ig.Separator()
        ig.Text("Files")
        ig.Text("Source files: $(status.total_files)")
        ig.Text("Cached: $(status.cached_files)")
        ig.Text("Fresh: $(status.fresh_files)")
        ig.Text("Stale: $(status.stale_files)")
        ig.Text("New: $(status.new_files)")
        ig.Text("Deleted: $(status.deleted_files)")
        status.error_files > 0 && ig.TextColored(
            (1.0, 0.35, 0.35, 1.0),
            "Errors: $(status.error_files)",
        )
    elseif cache_state == :error
        message = get(ui_state, :cache_error, "")
        !isempty(message) && ig.TextWrapped(message)
    end

    cache_errors = get(ui_state, :cache_errors, ProjectCacheFileError[])
    if !isempty(cache_errors)
        ig.Separator()
        ig.TextColored((1.0, 0.35, 0.35, 1.0), "File Errors")
        for (index, file_error) in enumerate(cache_errors)
            index > 20 && break
            ig.TextWrapped(file_error.path)
            ig.TextColored((1.0, 0.55, 0.35, 1.0), file_error.message)
            index < min(length(cache_errors), 20) && ig.Separator()
        end
        length(cache_errors) > 20 && ig.TextDisabled("$(length(cache_errors) - 20) more file errors")
    end

    ig.Separator()
    has_identity = identity isa ProjectCacheIdentity
    can_update = has_identity &&
        status isa ProjectCacheStatus &&
        get(ui_state, :cache_source_checked, false) &&
        cache_state != :missing &&
        cache_state != :error &&
        (status.stale_files > 0 || status.new_files > 0 || status.deleted_files > 0)
    can_build = has_identity && cache_state != :error
    can_scan = haskey(ui_state, :root_path) && !isempty(get(ui_state, :root_path, ""))
    can_press_scan = can_scan && !_cache_action_blocked(ui_state)

    !can_press_scan && ig.BeginDisabled()
    scan_label = _source_scan_running(ui_state) ? "Cancel Scan" : "Scan Source"
    if ig.Button(scan_label, (-1, 0))
        if _source_scan_running(ui_state)
            _request_source_scan_cancel!(ui_state)
        else
            proj = haskey(ui_state, :project) ?
                ui_state[:project] :
                _project_for_preference(get(ui_state, :project_preference, "auto"))
            _launch_source_scan_job!(ui_state, ui_state[:root_path], proj)
        end
    end
    !can_press_scan && ig.EndDisabled()

    (!can_update || running) && ig.BeginDisabled()
    if ig.Button("Update Cache", (-1, 0))
        _queue_cache_update!(ui_state)
    end
    (!can_update || running) && ig.EndDisabled()

    (!can_build || running) && ig.BeginDisabled()
    build_label = has_identity ? _cache_build_label(ui_state, identity, cache_state) : "Build Cache"
    if ig.Button(build_label, (-1, 0))
        _queue_cache_update!(ui_state; full_rebuild=true)
    end
    (!can_build || running) && ig.EndDisabled()

    if running
        if ig.Button(activity.cancel_label, (-1, 0))
            _request_cache_cancel!(ui_state)
        end
    end
end

function render_project_window(ui_state)
    get(ui_state, :show_project_window, false) || return

    open_ref = Ref(true)
    if ig.Begin("Project Settings###project_window", open_ref, ig.ImGuiWindowFlags_AlwaysAutoResize)
        if haskey(ui_state, :project)
            proj = ui_state[:project]
            ig.Text("Active: $(project_name(proj))")
        else
            ig.TextDisabled("No folder loaded yet")
        end

        ig.Separator()

        pref = get(ui_state, :project_preference, "auto")
        changed = false

        default_project = something(_default_project[])
        default_label = "Default ($(project_name(default_project)))"

        if ig.RadioButton(default_label, pref == "auto")
            ui_state[:project_preference] = "auto"
            changed = true
        end
        ig.SameLine()
        _helpmarker("Use the default project without trying to infer the project from the files in the folder.")

        for p in KNOWN_PROJECTS
            pn = project_name(p)
            if ig.RadioButton(pn, pref == pn)
                ui_state[:project_preference] = pn
                changed = true
            end
            ig.SameLine()
            _helpmarker(project_description(p))
        end

        if changed
            current_root = get(ui_state, :root_path, "")
            _persist_preferences!(ui_state; path=isempty(current_root) ? nothing : current_root)
            if !isempty(current_root)
                @info "Project preference changed to '$(ui_state[:project_preference])' - reloading cache"
                _open_project_path!(ui_state, current_root)
            end
        end
    end
    open_ref[] || (ui_state[:show_project_window] = false)
    ig.End()
end

function _setup_docking_layout!(dockspace_id)
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

    # Right column: top 3/4 plot, bottom 1/4 combined plots
    top_right    = Ref{UInt32}(0)
    bottom_right = Ref{UInt32}(0)
    ig.DockBuilderSplitNode(right[], ig.ImGuiDir_Up, 3/4, top_right, bottom_right)

    ig.DockBuilderDockWindow("Hierarchy",         top_left[])
    ig.DockBuilderDockWindow("Plot Area",          top_right[])
    ig.DockBuilderDockWindow("Information Panel",  bottom_left[])
    ig.DockBuilderDockWindow("Combined Plots",     bottom_right[])

    ig.DockBuilderFinish(dockspace_id)
end

function create_window_and_run_loop(root_path::Union{Nothing,String}=nothing; engine=nothing, spawn=1)
    ig.set_backend(:GlfwOpenGL3)
    ui_state = Dict{Symbol,Any}()
    _init_scan_state!(ui_state)
    _init_cache_state!(ui_state)
    _init_tag_state!(ui_state)
    _init_figure_script_state!(ui_state)
    _init_plot_state!(ui_state)
    ui_state[:_frame] = 0
    ctx = ig.CreateContext()
    io = ig.GetIO()
    io.ConfigFlags = unsafe_load(io.ConfigFlags) | ig.ImGuiConfigFlags_DockingEnable
    io.ConfigFlags = unsafe_load(io.ConfigFlags) | ig.ImGuiConfigFlags_ViewportsEnable
    io.ConfigFlags = unsafe_load(io.ConfigFlags) | ig.ImGuiConfigFlags_NavEnableKeyboard
    io.IniFilename = Ptr{Cchar}(C_NULL)   # disable layout persistence
    ig.StyleColorsDark()

    # Load persisted preferences
    prefs = _load_prefs()
    ui_state[:project_preference] = get(prefs, "project", "auto")
    ui_state[:recent_projects] = _parse_recent_projects(prefs)

    if root_path !== nothing && root_path != ""
        _open_project_path!(ui_state, root_path)
    end
    first_frame   = Ref(true)
    setup_layout  = Ref(true)
    ig.render(
        ctx;
        engine,
        window_size=(1920, 1080),
        window_title="Measurement Browser",
        opengl_version=v"3.3",
        spawn,
        wait_events=false,
        on_exit=() -> _print_perf_summary(ui_state),
    ) do
        ui_state[:_frame] += 1
        _poll_cache_events!(ui_state)
        _poll_source_scan_events!(ui_state)
        _poll_plot_events!(ui_state)
        _poll_figure_script_job_events!(ui_state)
        if first_frame[] && !haskey(ui_state, :_gl_info)
            ui_state[:_gl_info] = _collect_gl_info!()
            first_frame[] = false
        end
        dockspace_id = ig.DockSpaceOverViewport(0, ig.GetMainViewport())
        if setup_layout[]
            setup_layout[] = false
            _setup_docking_layout!(dockspace_id)
        end
        if !_scan_running(ui_state)
            _time!(ui_state, :plot_warmup) do
                _ensure_plot_runtime_warmed!(ui_state)
            end
        end
        render_selection_window(ui_state)
        render_project_window(ui_state)
        _time!(ui_state, :info) do
            render_info_window(ui_state)
        end
        _time!(ui_state, :figure_scripts) do
            render_figure_script_window(ui_state)
        end
        _time!(ui_state, :plot) do
            render_plot_window(ui_state)
        end
        _time!(ui_state, :extra_plots) do
            render_additional_plot_windows(ui_state)
        end
        _time!(ui_state, :perf_window) do
            render_perf_window(ui_state)
        end
        _time!(ui_state, :combined_plots) do
            render_combined_plots_window(ui_state)
        end
        # Show metadata guidance modal if needed
        render_device_info_modal(ui_state)
    end
end

function start_browser(root_path::Union{Nothing,String}=nothing)
    create_window_and_run_loop(root_path)
end
