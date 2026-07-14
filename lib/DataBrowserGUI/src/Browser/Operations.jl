using DataBrowserAPI:
    Project,
    DEFAULT_PROJECT,
    PROJECTS,
    project_name
using DataBrowserAPI.ItemIndex: collection_path_keys
using DataBrowserCache: ProjectCacheSchemaError
import DataBrowserCore.Workspace
using DataBrowserCore.Workspace:
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
    rebuild_cache::Bool=false,
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
        project = _project_for_preference(state.project_preference)
    else
        state.project_preference = project_name(project)
    end
    previous_workspace = state.workspace
    reopen_options = previous_workspace isa Workspace.Workspace ?
        previous_workspace.open_options :
        (;)
    _attach_workspace!(
        state,
        open_workspace(project, norm_path; reopen_options..., rebuild=rebuild_cache),
    )
    return nothing
end

"""
Make an already-opened workspace the browser's current one, loading its saved view, tag state, and
figure-script context. Shared by `_open_project_path!` (which opens the workspace itself) and
`open_browser` (which is handed a caller-owned workspace).
"""
function _attach_workspace!(
    state::BrowserState,
    workspace::Workspace.Workspace,
)::Nothing
    source = workspace.source
    source_root = hasproperty(source, :root_path) ? source.root_path : ""
    previous_workspace = state.workspace
    previous_workspace isa Workspace.Workspace && previous_workspace !== workspace &&
        close_workspace!(previous_workspace)
    _reset_extensions!(state)
    state.workspace = workspace
    view = isempty(source_root) ? PersistedProjectView() : _load_project_view(source_root)
    project = project_name(workspace.project)
    !isempty(view.project) && view.project != project &&
        (view = PersistedProjectView(project=project))
    _load_tag_state_for_root!(state, _annotation_root(workspace))
    _apply_project_view!(state, view)
    state.saved_project_view = view
    if workspace.cache.disk_error isa ProjectCacheSchemaError
        state.cache_rebuild_modal = true
        state.cache_rebuild_path = source_root
        state.cache_rebuild_project = workspace.project
        state.cache_rebuild_error = sprint(showerror, workspace.cache.disk_error)
    end
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
        for item in values(workspace.index.items)
        if item.source_item_id == id
    ]
    isempty(items) && return false

    collections = workspace.index.collections
    collection_ids = unique(String[
        item.collection_key === nothing ?
            ROOT_COLLECTION_SELECTION_ID :
            collections.records[item.collection_key].id
        for item in items
    ])
    expanded_ids = copy(state.expanded_collection_ids)
    for item in items
        item.collection_key === nothing && continue
        path = collection_path_keys(collections, item.collection_key)
        for key in path[1:end-1]
            parent_id = collections.records[key].id
            parent_id in expanded_ids || push!(expanded_ids, parent_id)
        end
    end

    state.expanded_collection_ids = expanded_ids
    workspace.selection.collection_ids = collection_ids
    workspace.selection.item_ids = [item.id for item in items]
    state.scroll_to_collection_id = isempty(collection_ids) ? nothing : first(collection_ids)
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
