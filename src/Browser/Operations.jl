using ..Projects:
    Project,
    DEFAULT_PROJECT,
    PROJECTS,
    project_name
using ..ItemIndex: collection_path_key
import ..Workspace
using ..Workspace:
    close_workspace!,
    open_workspace

"""Return the project selected by the saved project preference."""
function _project_for_preference(pref::AbstractString)::Project
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
    project::Union{Nothing,Project}=nothing,
    persist::Bool=true,
)::Nothing
    norm_path = _normalize_project_path(path)
    if project === nothing && state.project_locked
        workspace = state.workspace
        if workspace isa Workspace.Workspace
            project = workspace.project
        else
            error("Cannot open a new folder before a project workspace exists")
        end
    end
    if project === nothing
        state.project_preference = _project_preference_for_path(state, norm_path)
        project = _project_for_preference(state.project_preference)
    else
        state.project_preference = project_name(project)
    end
    _attach_workspace!(state, open_workspace(project, norm_path); persist)
    return nothing
end

"""
Make an already-opened workspace the browser's current one, loading its saved view, tag state, and
figure-script context. Shared by `_open_project_path!` (which opens the workspace itself) and
`open_browser` (which is handed a caller-owned workspace).
"""
function _attach_workspace!(
    state::BrowserState,
    workspace::Workspace.Workspace;
    persist::Bool=true,
)::Nothing
    source = workspace.source
    source_root = hasproperty(source, :root_path) ? source.root_path : ""
    previous_workspace = state.workspace
    previous_workspace isa Workspace.Workspace && previous_workspace !== workspace &&
        close_workspace!(previous_workspace)
    state.plots = PlotState(debug=state.plots.debug)
    state.workspace = workspace
    view = isempty(source_root) ? PersistedProjectView() : _load_project_view(source_root)
    project = project_name(workspace.project)
    !isempty(view.project) && view.project != project &&
        (view = PersistedProjectView(project=project))
    _load_tag_state_for_root!(state, _annotation_root(workspace))
    _apply_project_view!(state, view)
    state.saved_project_view = view
    persist && !isempty(source_root) && _persist_preferences!(state; path=source_root)
    return nothing
end

"""Select and reveal every item produced by one source item."""
function select_source_item!(
    state::BrowserState,
    source_item_id::AbstractString,
)::Bool
    workspace = state.workspace
    workspace isa Workspace.Workspace || return false
    id = String(source_item_id)
    items = [
        item
        for item in workspace.index.hierarchy.all_items
        if item.source_item_id == id
    ]
    isempty(items) && return false

    collection_paths =
        unique([collection_path_key(item.collection) for item in items])
    expanded_paths = copy(state.expanded_collection_paths)
    for collection_path in collection_paths
        parts = split(collection_path, '/')
        for depth in 1:(length(parts) - 1)
            parent_path = join(parts[1:depth], '/')
            parent_path in expanded_paths || push!(expanded_paths, parent_path)
        end
    end

    state.expanded_collection_paths = expanded_paths
    workspace.selection.collection_paths = collection_paths
    workspace.selection.item_ids = [item.id for item in items]
    state.scroll_to_collection_path = first(collection_paths)
    state.scroll_to_item_id = first(items).id
    return true
end

"""Stop browser and workspace work before the render loop exits."""
function _shutdown_background_jobs!(state::BrowserState)::Nothing
    state.shutdown_complete && return nothing
    workspace = state.workspace
    workspace isa Workspace.Workspace && close_workspace!(workspace)
    state.shutdown_complete = true
    return nothing
end
