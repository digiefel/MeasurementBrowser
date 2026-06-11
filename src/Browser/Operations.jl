using ..Project:
    AbstractProject,
    DEFAULT_PROJECT,
    PROJECTS,
    project_name
using ..MeasurementIndex: device_path_key
import ..Workspace
using ..Workspace:
    close_workspace!,
    open_workspace

"""Return the project selected by the saved project preference."""
function _project_for_preference(pref::AbstractString)::AbstractProject
    pref == "auto" && return something(DEFAULT_PROJECT[])
    for project in PROJECTS
        project_name(project) == pref && return project
    end
    error("Unknown project preference '$pref'")
end

"""Replace the current workspace with the project selected for one source root."""
function _open_project_path!(
    state::BrowserState,
    path::String;
    project::Union{Nothing,AbstractProject}=nothing,
    persist::Bool=true,
)::Nothing
    norm_path = _normalize_project_path(path)
    if project === nothing && state.project_locked
        project = state.workspace.project
    end
    if project === nothing
        state.project_preference = _project_preference_for_path(state, norm_path)
        project = _project_for_preference(state.project_preference)
    else
        state.project_preference = project_name(project)
    end
    previous_workspace = state.workspace
    previous_root = previous_workspace isa Workspace.Workspace ?
        previous_workspace.root_path :
        ""
    _cancel_figure_script_job!(state)
    previous_workspace isa Workspace.Workspace && close_workspace!(previous_workspace)
    state.plots = PlotState(debug=state.plots.debug)
    previous_root == norm_path ||
        _reset_figure_script_state!(state, norm_path)
    view = _load_project_view(norm_path)
    !isempty(view.project) && view.project != project_name(project) &&
        (view = PersistedProjectView(project=project_name(project)))
    _load_tag_state_for_root!(state, norm_path)
    state.workspace = open_workspace(project, norm_path)
    _apply_project_view!(state, view)
    state.saved_project_view = view
    _invalidate_figure_script_scan_cache!(state)
    persist && _persist_preferences!(state; path=norm_path)
    return nothing
end

"""Select and reveal every measurement produced by one source file."""
function select_source_file!(
    state::BrowserState,
    filepath::AbstractString,
)::Bool
    workspace = state.workspace
    workspace isa Workspace.Workspace || return false
    path = String(filepath)
    measurements = [
        measurement
        for measurement in workspace.index.hierarchy.all_measurements
        if measurement.filepath == path
    ]
    isempty(measurements) && return false

    device_paths =
        unique([device_path_key(measurement.device_info) for measurement in measurements])
    expanded_paths = copy(state.expanded_device_paths)
    for device_path in device_paths
        parts = split(device_path, '/')
        for depth in 1:(length(parts) - 1)
            parent_path = join(parts[1:depth], '/')
            parent_path in expanded_paths || push!(expanded_paths, parent_path)
        end
    end

    state.expanded_device_paths = expanded_paths
    workspace.selection.device_paths = device_paths
    workspace.selection.measurement_ids =
        [measurement.unique_id for measurement in measurements]
    state.scroll_to_device_path = first(device_paths)
    state.scroll_to_measurement_id = first(measurements).unique_id
    return true
end

"""Stop browser and workspace work before the render loop exits."""
function _shutdown_background_jobs!(state::BrowserState)::Nothing
    state.shutdown_complete && return nothing
    _cancel_figure_script_job!(state)
    workspace = state.workspace
    workspace isa Workspace.Workspace && close_workspace!(workspace)
    state.shutdown_complete = true
    return nothing
end
