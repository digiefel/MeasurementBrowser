function render_info_window(ui_state)
    proj = ui_state[:project]
    selected_devices = get(ui_state, :selected_devices, HierarchyNode[])
    selected_measurements = get(ui_state, :selected_measurements, MeasurementInfo[])
    if ig.Begin("Information Panel")
        flags = ig.ImGuiTableFlags_Borders | ig.ImGuiTableFlags_RowBg | ig.ImGuiTableFlags_ScrollY
        ig.BeginTable("info_cols", 2, flags)
        ig.TableSetupColumn("Device")
        ig.TableSetupColumn("Measurement")
        ig.TableHeadersRow()
        ig.TableNextRow()
        ig.TableNextColumn()

        if length(selected_devices) == 1
            meas_vec = selected_devices[1].measurements
            sel_name = join(get(ui_state, :selected_path, [""]), "/")
            ig.Text("Location: $sel_name")
            ig.Separator()
            stats = begin
                try
                    get_measurements_stats(meas_vec, proj)
                catch err
                    @warn "Failed to compute stats" error = err
                    Dict{Symbol,Any}()
                end
            end
            if !isempty(stats)
                ig.Text("Stats")
                ig.BulletText("Total: $(stats[:total_measurements])")
                ig.BulletText("Types: $(join(stats[:measurement_types], ", "))")
                if haskey(stats, :first_measurement)
                    ig.BulletText("First: $(stats[:first_measurement]) ")
                    ig.BulletText("Last:  $(stats[:last_measurement]) ")
                end
            else
                ig.TextDisabled("No stats available")
            end
            ig.Separator()
            # Device-level metadata
            if !isempty(meas_vec)
                dev_meta = first(meas_vec).device_info.parameters
                if !isempty(dev_meta)
                    ig.Text("Device metadata")
                    for (k, v) in dev_meta
                        ig.BulletText("$(k): $(v)")
                    end
                else
                    ig.TextDisabled("No metadata parameters found")
                end
            end
        elseif isempty(selected_devices)
            ig.TextDisabled("Select a device to see details")
        else
            ig.TextDisabled("Select a single device to see details")
        end

        ig.TableNextColumn()
        if length(selected_measurements) == 1
            m = selected_measurements[1]
            ig.Text("Title: $(m.clean_title)")
            ig.Separator()
            ig.BulletText("Type: $(kind_label(proj, m.measurement_kind))")
            ig.BulletText("Timestamp: $(m.timestamp)")
            ig.BulletText("Filename:")
            ig.SameLine()
            ig.TextLinkOpenURL(m.filename, m.filepath)
            ig.Separator()
            if !isempty(m.parameters)
                ig.Text("Parameters")
                for (k, v) in m.parameters
                    ig.BulletText("$(k) = $(v)")
                end
            else
                ig.TextDisabled("No parameters extracted")
            end

            # ---- data-derived statistics (PUND measurements only) ---------
            if m.measurement_kind == :pund
                ig.Separator()
                ig.Text("Statistics")
                cache = get(ui_state, :computed_stats_cache, Dict{Tuple{String,Int},Dict{Symbol,Any}}())
                cache_key = (m.filepath, get(m.parameters, :fatigue_cycle, 0))
                stats = get!(cache, cache_key) do
                    cached = _compute_cached_pund_stats(ui_state, m)
                    cached === nothing ?
                        compute_pund_stats(m.filepath, m.parameters, m.device_info.parameters) :
                        cached
                end
                if !isempty(stats)
                    order = [:voltage_max_V, :voltage_baseline_V, :voltage_min_V, :frequency_kHz, :Pr_max_uCcm2]
                    for key in order
                        haskey(stats, key) || continue
                        v = stats[key]
                        label = replace(string(key), "_" => " ")
                        if v === nothing
                            ig.BulletText("$label: —")
                        else
                            ig.BulletText("$label = $v")
                        end
                    end
                else
                    ig.TextDisabled("Could not compute statistics")
                end
            end

        elseif isempty(selected_measurements)
            ig.TextDisabled("Select a measurement to view details")
        else
            ig.TextDisabled("Select a single measurement to view details")
        end
        ig.EndTable()
    end
    ig.End()
end

function _compute_cached_pund_stats(ui_state, measurement::MeasurementInfo)
    identity = get(ui_state, :cache_identity, nothing)
    identity isa ProjectCacheIdentity || return nothing
    try
        analyzed = _measurement_group_for_cached_plot(identity, measurement)
        return compute_pund_stats_from_analyzed_plot(analyzed, measurement.device_info.parameters)
    catch err
        @debug "Could not compute PUND statistics from cache" exception=(err, catch_backtrace())
        return nothing
    end
end

function _figure_script_output_path(ui_state)
    output_directory = _buffer_string(ui_state[:figure_script_output_dir_buffer])
    script_name = _buffer_string(ui_state[:figure_script_name_buffer])
    try
        return figure_script_path(output_directory, script_name)
    catch err
        err isa FigureScriptValidationError || rethrow()
        return nothing
    end
end

function _write_figure_script_from_ui!(ui_state; overwrite::Bool=false)
    root_path = get(ui_state, :root_path, "")
    isempty(root_path) && throw(FigureScriptValidationError("Open a project folder before generating a figure script"))
    project = get(ui_state, :project, nothing)
    project isa AbstractProject || error("Figure script generation requires an active project")
    path = write_figure_script(
        _buffer_string(ui_state[:figure_script_output_dir_buffer]),
        root_path,
        project,
        _buffer_string(ui_state[:figure_script_name_buffer]),
        copy(_figure_script_groups(ui_state)),
        _current_scan_measurements(ui_state);
        overwrite=overwrite,
    )
    ui_state[:figure_script_overwrite_confirm] = ""
    _set_figure_script_status!(ui_state, "Wrote $(basename(path))")
    return path
end

function _render_figure_script_group_tooltip(proj, preview_measurements::Vector{MeasurementInfo})
    ig.BeginItemTooltip() || return
    for measurement in preview_measurements[1:min(6, end)]
        ig.BulletText("$(display_label(proj, measurement))")
    end
    length(preview_measurements) > 6 && ig.TextDisabled("...")
    ig.EndTooltip()
end

function render_figure_script_window(ui_state)
    get(ui_state, :show_figure_script_window, false) || return

    proj = ui_state[:project]
    selected_measurements = _selected_measurements_in_panel_order(ui_state)
    selected_count = length(selected_measurements)
    job_running = _figure_script_job_running(ui_state)
    groups = _figure_script_groups(ui_state)
    group_matches = isempty(groups) ? Dict{String,Vector{MeasurementInfo}}() : _ensure_figure_script_group_matches!(ui_state)
    output_path = _figure_script_output_path(ui_state)
    overwrite_path = get(ui_state, :figure_script_overwrite_confirm, "")
    if output_path === nothing || output_path != overwrite_path
        ui_state[:figure_script_overwrite_confirm] = ""
        overwrite_path = ""
    end

    open_ref = Ref(true)
    ig.SetNextWindowSize((520, 430), ig.ImGuiCond_FirstUseEver)
    if ig.Begin("Figure Script", open_ref, ig.ImGuiWindowFlags_NoDocking)
        ig.TextDisabled("$selected_count measurements selected")
        ig.Separator()

        ig.Text("Output Directory")
        output_dir_changed = ig.InputText(
            "##figure_script_output_dir",
            ui_state[:figure_script_output_dir_buffer],
            length(ui_state[:figure_script_output_dir_buffer]),
        )
        output_dir_changed && _persist_current_project_preferences!(ui_state)
        ig.SameLine()
        if ig.Button("Choose...")
            selected_dir = pick_folder()
            if !isnothing(selected_dir) && !isempty(selected_dir)
                _set_buffer_string!(ui_state[:figure_script_output_dir_buffer], _normalize_project_path(selected_dir))
                _persist_current_project_preferences!(ui_state)
            end
        end
        ig.SameLine()
        _helpmarker("Scripts are written to this directory.")

        ig.Spacing()
        ig.Text("Script Name")
        ig.InputText("##figure_script_name", ui_state[:figure_script_name_buffer], length(ui_state[:figure_script_name_buffer]))
        ig.TextDisabled(output_path === nothing ? "<choose output directory and enter script name>" : output_path)

        ig.Spacing()
        ig.Text("Group Name")
        ig.InputText("##figure_group_name", ui_state[:figure_script_group_name_buffer], length(ui_state[:figure_script_group_name_buffer]))

        ig.Spacing()
        ig.Text("Groups")
        ig.SameLine()
        _helpmarker("Each group becomes one entry in the generated script: data[\"group_name\"]::Vector{MeasurementBrowser.FigureMeasurement}")
        if ig.BeginChild("figure_script_groups", (0, 190), true)
            if isempty(groups)
                ig.TextDisabled("No groups yet")
            else
                for (index, group) in enumerate(groups)
                    preview_measurements = get(group_matches, group.name, MeasurementInfo[])
                    match_count = length(preview_measurements)
                    label = "$(group.name) ($(match_count))###figure_group_$index"
                    if ig.Selectable(label, get(ui_state, :figure_script_selected_group, 0) == index)
                        _set_selected_figure_script_group!(ui_state, index)
                    end
                    ig.IsItemHovered() && _render_figure_script_group_tooltip(proj, preview_measurements)
                end
            end
        end
        ig.EndChild()

        group_selected = _selected_figure_script_group(ui_state) !== nothing
        can_apply_selection = selected_count > 0

        (!can_apply_selection || job_running) && ig.BeginDisabled()
        if ig.Button("New Group From Selection")
            _clear_figure_script_messages!(ui_state)
            try
                _create_figure_script_group_from_selection!(ui_state)
            catch err
                if err isa FigureScriptValidationError
                    _set_figure_script_error!(ui_state, err)
                else
                    rethrow()
                end
            end
        end
        (!can_apply_selection || job_running) && ig.EndDisabled()

        (!group_selected || job_running) && ig.BeginDisabled()
        ig.SameLine()
        if ig.Button("Rename Group")
            _clear_figure_script_messages!(ui_state)
            try
                _rename_selected_figure_script_group!(ui_state)
            catch err
                if err isa FigureScriptValidationError
                    _set_figure_script_error!(ui_state, err)
                else
                    rethrow()
                end
            end
        end
        ig.SameLine()
        if ig.Button("Delete Group")
            _clear_figure_script_messages!(ui_state)
            try
                _delete_selected_figure_script_group!(ui_state)
            catch err
                if err isa FigureScriptValidationError
                    _set_figure_script_error!(ui_state, err)
                else
                    rethrow()
                end
            end
        end
        (!group_selected || job_running) && ig.EndDisabled()

        (!group_selected || !can_apply_selection || job_running) && ig.BeginDisabled()
        if group_selected || can_apply_selection
            ig.SameLine()
        end
        if ig.Button("Add Selection")
            _clear_figure_script_messages!(ui_state)
            try
                _add_selection_to_figure_script_group!(ui_state)
            catch err
                if err isa FigureScriptValidationError
                    _set_figure_script_error!(ui_state, err)
                else
                    rethrow()
                end
            end
        end
        ig.SameLine()
        if ig.Button("Remove Selection")
            _clear_figure_script_messages!(ui_state)
            try
                _remove_selection_from_figure_script_group!(ui_state)
            catch err
                if err isa FigureScriptValidationError
                    _set_figure_script_error!(ui_state, err)
                else
                    rethrow()
                end
            end
        end
        (!group_selected || !can_apply_selection || job_running) && ig.EndDisabled()

        ig.Separator()
        job_running && ig.BeginDisabled()
        if ig.Button("Generate Script")
            _clear_figure_script_messages!(ui_state)
            try
                _write_figure_script_from_ui!(ui_state)
            catch err
                if err isa FigureScriptExistsError
                    ui_state[:figure_script_overwrite_confirm] = err.path
                    _set_figure_script_error!(ui_state, err)
                elseif err isa FigureScriptValidationError || err isa FigureScriptIOError
                    _set_figure_script_error!(ui_state, err)
                else
                    rethrow()
                end
            end
        end
        job_running && ig.EndDisabled()

        if !isempty(overwrite_path)
            job_running && ig.BeginDisabled()
            ig.SameLine()
            if ig.Button("Overwrite Existing")
                _clear_figure_script_messages!(ui_state)
                try
                    _write_figure_script_from_ui!(ui_state; overwrite=true)
                catch err
                    if err isa FigureScriptValidationError || err isa FigureScriptIOError
                        _set_figure_script_error!(ui_state, err)
                    else
                        rethrow()
                    end
                end
            end
            ig.SameLine()
            if ig.Button("Cancel Overwrite")
                ui_state[:figure_script_overwrite_confirm] = ""
            end
            job_running && ig.EndDisabled()
        end

        error_message = get(ui_state, :figure_script_error, "")
        !isempty(error_message) && ig.TextColored((1.0, 0.45, 0.45, 1.0), error_message)
        status_message = get(ui_state, :figure_script_status, "")
        !isempty(status_message) && ig.TextColored((0.45, 0.85, 0.55, 1.0), status_message)
        job_running && ig.TextDisabled("Figure-script worker is running...")
    end
    open_ref[] || (ui_state[:show_figure_script_window] = false)
    ig.End()
end

# ------------------------------------------------------------------
# Modal for missing device metadata (shown each scan when missing)
# ------------------------------------------------------------------
function render_device_info_modal(ui_state)
    # Reset dismissal when root path changes
    current_root = get(ui_state, :root_path, "")
    if get(ui_state, :_modal_last_root_path, "") != current_root
        ui_state[:_modal_last_root_path] = current_root
        ui_state[:dev_info_modal] = true
    end
    # always center
    center = ig.ImVec2(0.5, 0.5)
    @c ig.ImGuiViewport_GetCenter(&center, ig.GetMainViewport())
    ig.SetNextWindowPos(center, ig.ImGuiCond_Always, (0.5, 0.5))

    # Show modal if: missing metadata and user hasn't dismissed it this scan
    if get(ui_state, :dev_info_modal, true) && !get(ui_state, :has_device_metadata, true)
        ig.OpenPopup("Device Metadata Missing")
    end

    opened = get(ui_state, :dev_info_modal, true)

    if @c ig.BeginPopupModal("Device Metadata Missing", &opened, ig.ImGuiWindowFlags_AlwaysAutoResize)
        ig.Text("No device metadata file (device_info.txt) was found.")
        ig.Separator()
        ig.TextWrapped("Create a simple text file named device_info.txt in the TOP folder you opened to add extra info (area, thickness, notes, etc.) for each device.")
        ig.Spacing()
        ig.TextWrapped("How to do it:")
        ig.BulletText("Create a new text file: device_info.txt")
        ig.BulletText("First line (header): device_path, area_um2, t_HZO_nm, ...")
        ig.BulletText("Add a column for each property you want to track.")
        ig.BulletText("Add one line per scope. device_path is matched as an exact slash-separated unit sequence.")
        ig.BulletText("Single units match only exact units, never substrings inside another unit.")
        ig.BulletText("Examples: D1, A9, VI/25um, RuO2test/A9/VI/D1, RuO2test_A10/VI/25um/D1")
        ig.BulletText("Longer exact unit sequences override shorter matches when they set the same field.")
        ig.Spacing()
        ig.TextDisabled("Example:")
        ig.TextDisabled("device_path,   area_um2,   t_HZO_nm,   notes,   active")
        ig.TextDisabled("D1,    12.5,   7.0,   all exact D1 units,   true")
        ig.TextDisabled("RuO2test_A10/VI/25um/D1,    12.4,   7.0,   exact device,  true")
        ig.Spacing()
        ig.TextWrapped("Save the file, then rescan or reopen the folder to load these values.")
        ig.Spacing()
        if ig.Button("Got it")
            opened = false
            ig.CloseCurrentPopup()
        end
        ig.EndPopup()
    end
    ui_state[:dev_info_modal] = opened
end
