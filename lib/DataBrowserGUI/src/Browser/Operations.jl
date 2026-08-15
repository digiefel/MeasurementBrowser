using DataBrowserAPI:
    AbstractProject,
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
function _project_for_preference(pref::AbstractString)::AbstractProject
    pref == "auto" && return something(DEFAULT_PROJECT[])
    for project in PROJECTS
        project_name(project) == pref && return project
    end
    error("Unknown project preference '$pref'")
end

"""Project used when the browser reopens a workspace."""
function _open_project(state::BrowserState)::AbstractProject
    if state.project_locked
        workspace = state.workspace
        workspace isa Workspace.Workspace ||
            error("Cannot reopen before a workspace exists")
        return workspace.project
    end
    return _project_for_preference(state.project_preference)
end

"""Open a new workspace for the current source. Used for project change and cache rebuild."""
function _reopen_workspace!(
    state::BrowserState;
    rebuild_cache::Bool=false,
)::Nothing
    previous = state.workspace
    previous isa Workspace.Workspace || error("Cannot reopen before a workspace exists")
    _attach_workspace!(
        state,
        open_workspace(
            _open_project(state),
            copy(previous.source);
            rebuild=rebuild_cache,
            cache=previous.disk_cache,
            background_processing=previous.background_processing,
        ),
    )
    return nothing
end

"""
Make an already-opened workspace the browser's current one, loading its saved view, tag state, and
figure-script context. Shared by `_reopen_workspace!` and `open_browser`.
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
    key = Workspace.source_item_key(workspace, String(source_item_id))
    key === nothing && return false
    items = [
        item
        for item in values(workspace.index.items)
        if item.source_item_key == key
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
