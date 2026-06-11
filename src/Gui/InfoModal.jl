import CImGui as ig
import CImGui.CSyntax: @c
using NativeFileDialog: pick_folder

using ..Project:
    AbstractProject,
    display_label,
    kind_label
using ..MeasurementIndex: MeasurementInfo
import ..Workspace
using ..MeasurementBrowser:
    FigureScriptExistsError,
    FigureScriptIOError,
    FigureScriptValidationError,
    figure_script_path,
    write_figure_script

"""Render device and measurement details for the visible workspace selection."""
function render_info_window(state::BrowserState)::Nothing
    workspace = state.workspace
    if !(workspace isa Workspace.Workspace)
        if ig.Begin("Information Panel")
            ig.TextDisabled("Open a project folder to inspect measurements")
        end
        ig.End()
        return nothing
    end
    proj = workspace.project
    selected_devices, selected_measurements, selected_path =
        _project_visible_selection(state)
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
            sel_name = join(selected_path, "/")
            ig.Text("Location: $sel_name")
            ig.Separator()
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
                for k in sort!(collect(keys(m.parameters)); by=String)
                    v = m.parameters[k]
                    ig.BulletText("$(k) = $(v)")
                end
            else
                ig.TextDisabled("No parameters extracted")
            end

            if !isempty(m.stats)
                ig.Separator()
                ig.Text("Statistics")
                for k in sort!(collect(keys(m.stats)); by=String)
                    v = m.stats[k]
                    ig.BulletText("$(k) = $(v)")
                end
            else
                ig.TextDisabled("No stats computed")
            end

        elseif isempty(selected_measurements)
            ig.TextDisabled("Select a measurement to view details")
        else
            ig.TextDisabled("Select a single measurement to view details")
        end
        ig.EndTable()
    end
    ig.End()
    return nothing
end

"""Return the current figure-script output path when both inputs are valid."""
function _figure_script_output_path(
    state::BrowserState,
)::Union{Nothing,String}
    output_directory = _buffer_string(state.figure_scripts.output_dir_buffer)
    script_name = _buffer_string(state.figure_scripts.script_name_buffer)
    try
        return figure_script_path(output_directory, script_name)
    catch err
        err isa FigureScriptValidationError || rethrow()
        return nothing
    end
end

"""Write the configured figure script and update its browser status."""
function _write_figure_script_from_ui!(
    state::BrowserState;
    overwrite::Bool=false,
)::String
    workspace = state.workspace::Workspace.Workspace
    path = write_figure_script(
        _buffer_string(state.figure_scripts.output_dir_buffer),
        workspace.root_path,
        workspace.project,
        _buffer_string(state.figure_scripts.script_name_buffer),
        copy(state.figure_scripts.groups),
        _current_scan_measurements(state);
        overwrite=overwrite,
    )
    state.figure_scripts.overwrite_confirm = ""
    _set_figure_script_status!(state, "Wrote $(basename(path))")
    return path
end

"""Show a short preview of measurements matched by one figure-script group."""
function _render_figure_script_group_tooltip(
    project::AbstractProject,
    preview_measurements::Vector{MeasurementInfo},
)::Nothing
    ig.BeginItemTooltip() || return nothing
    for measurement in preview_measurements[1:min(6, end)]
        ig.BulletText("$(display_label(project, measurement))")
    end
    length(preview_measurements) > 6 && ig.TextDisabled("...")
    ig.EndTooltip()
    return nothing
end

"""Render the temporary figure-script editor until Workflow replaces it."""
function render_figure_script_window(state::BrowserState)::Nothing
    figure_scripts = state.figure_scripts
    figure_scripts.visible || return nothing

    workspace = state.workspace::Workspace.Workspace
    proj = workspace.project
    selected_measurements = _selected_measurements_in_panel_order(state)
    selected_count = length(selected_measurements)
    job_running = _figure_script_job_running(state)
    groups = figure_scripts.groups
    group_matches = isempty(groups) ?
        Dict{String,Vector{MeasurementInfo}}() :
        _ensure_figure_script_group_matches!(state)
    output_path = _figure_script_output_path(state)
    overwrite_path = figure_scripts.overwrite_confirm
    if output_path === nothing || output_path != overwrite_path
        figure_scripts.overwrite_confirm = ""
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
            figure_scripts.output_dir_buffer,
            length(figure_scripts.output_dir_buffer),
        )
        output_dir_changed && _persist_current_project_preferences!(state)
        ig.SameLine()
        if ig.Button("Choose...")
            selected_dir = pick_folder()
            if !isnothing(selected_dir) && !isempty(selected_dir)
                _set_buffer_string!(
                    figure_scripts.output_dir_buffer,
                    _normalize_project_path(selected_dir),
                )
                _persist_current_project_preferences!(state)
            end
        end
        ig.SameLine()
        _helpmarker("Scripts are written to this directory.")

        ig.Spacing()
        ig.Text("Script Name")
        ig.InputText(
            "##figure_script_name",
            figure_scripts.script_name_buffer,
            length(figure_scripts.script_name_buffer),
        )
        ig.TextDisabled(output_path === nothing ? "<choose output directory and enter script name>" : output_path)

        ig.Spacing()
        ig.Text("Group Name")
        ig.InputText(
            "##figure_group_name",
            figure_scripts.group_name_buffer,
            length(figure_scripts.group_name_buffer),
        )

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
                    if ig.Selectable(label, figure_scripts.selected_group == index)
                        _set_selected_figure_script_group!(state, index)
                    end
                    ig.IsItemHovered() && _render_figure_script_group_tooltip(proj, preview_measurements)
                end
            end
        end
        ig.EndChild()

        group_selected = _selected_figure_script_group(state) !== nothing
        can_apply_selection = selected_count > 0

        (!can_apply_selection || job_running) && ig.BeginDisabled()
        if ig.Button("New Group From Selection")
            _clear_figure_script_messages!(state)
            try
                _create_figure_script_group_from_selection!(state)
            catch err
                if err isa FigureScriptValidationError
                    _set_figure_script_error!(state, err)
                else
                    rethrow()
                end
            end
        end
        (!can_apply_selection || job_running) && ig.EndDisabled()

        (!group_selected || job_running) && ig.BeginDisabled()
        ig.SameLine()
        if ig.Button("Rename Group")
            _clear_figure_script_messages!(state)
            try
                _rename_selected_figure_script_group!(state)
            catch err
                if err isa FigureScriptValidationError
                    _set_figure_script_error!(state, err)
                else
                    rethrow()
                end
            end
        end
        ig.SameLine()
        if ig.Button("Delete Group")
            _clear_figure_script_messages!(state)
            try
                _delete_selected_figure_script_group!(state)
            catch err
                if err isa FigureScriptValidationError
                    _set_figure_script_error!(state, err)
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
            _clear_figure_script_messages!(state)
            try
                _add_selection_to_figure_script_group!(state)
            catch err
                if err isa FigureScriptValidationError
                    _set_figure_script_error!(state, err)
                else
                    rethrow()
                end
            end
        end
        ig.SameLine()
        if ig.Button("Remove Selection")
            _clear_figure_script_messages!(state)
            try
                _remove_selection_from_figure_script_group!(state)
            catch err
                if err isa FigureScriptValidationError
                    _set_figure_script_error!(state, err)
                else
                    rethrow()
                end
            end
        end
        (!group_selected || !can_apply_selection || job_running) && ig.EndDisabled()

        ig.Separator()
        job_running && ig.BeginDisabled()
        if ig.Button("Generate Script")
            _clear_figure_script_messages!(state)
            try
                _write_figure_script_from_ui!(state)
            catch err
                if err isa FigureScriptExistsError
                    figure_scripts.overwrite_confirm = err.path
                    _set_figure_script_error!(state, err)
                elseif err isa FigureScriptValidationError || err isa FigureScriptIOError
                    _set_figure_script_error!(state, err)
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
                _clear_figure_script_messages!(state)
                try
                    _write_figure_script_from_ui!(state; overwrite=true)
                catch err
                    if err isa FigureScriptValidationError || err isa FigureScriptIOError
                        _set_figure_script_error!(state, err)
                    else
                        rethrow()
                    end
                end
            end
            ig.SameLine()
            if ig.Button("Cancel Overwrite")
                figure_scripts.overwrite_confirm = ""
            end
            job_running && ig.EndDisabled()
        end

        error_message = figure_scripts.error
        !isempty(error_message) && ig.TextColored((1.0, 0.45, 0.45, 1.0), error_message)
        status_message = figure_scripts.status
        !isempty(status_message) && ig.TextColored((0.45, 0.85, 0.55, 1.0), status_message)
        job_running && ig.TextDisabled("Figure-script worker is running...")
    end
    open_ref[] || (figure_scripts.visible = false)
    ig.End()
    return nothing
end

"""Explain how to add device metadata when the current source root has none."""
function render_device_info_modal(state::BrowserState)::Nothing
    workspace = state.workspace
    workspace isa Workspace.Workspace || return nothing
    # Reset dismissal when root path changes
    current_root = workspace.root_path
    if state.modal_root_path != current_root
        state.modal_root_path = current_root
        state.device_info_modal = true
    end
    # always center
    center = ig.ImVec2(0.5, 0.5)
    @c ig.ImGuiViewport_GetCenter(&center, ig.GetMainViewport())
    ig.SetNextWindowPos(center, ig.ImGuiCond_Always, (0.5, 0.5))

    # Show modal if: missing metadata and user hasn't dismissed it this scan
    if state.device_info_modal &&
       !workspace.index.hierarchy.has_device_metadata
        ig.OpenPopup("Device Metadata Missing")
    end

    opened = state.device_info_modal

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
    state.device_info_modal = opened
    return nothing
end
