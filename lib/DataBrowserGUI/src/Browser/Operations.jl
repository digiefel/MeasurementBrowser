using DataBrowserAPI:
    project_name,
    source_id
using DataBrowserAPI.ItemIndex: collection_path_keys
using DataBrowserCache: ProjectCacheSchemaError
import DataBrowserCore.Workspace
using DataBrowserCore.Workspace:
    close_workspace!

"""Reload tags and the persisted view when the workspace source identity changes."""
function _follow_source_identity!(
    state::BrowserState,
    workspace::Workspace.Workspace,
)::Nothing
    current = source_id(workspace.source)
    current == state.loaded_source_id && return nothing
    isempty(state.loaded_source_id) || _reset_extensions!(state)
    source_root = hasproperty(workspace.source, :root_path) ? workspace.source.root_path : ""
    view = isempty(source_root) ? PersistedProjectView() : _load_project_view(source_root)
    project = project_name(workspace.project)
    !isempty(view.project) && view.project != project &&
        (view = PersistedProjectView(project=project))
    _load_tag_state_for_root!(state, _annotation_root(workspace))
    _apply_project_view!(state, view)
    state.saved_project_view = view
    state.loaded_source_id = current
    return nothing
end

"""Open the cache-rebuild modal when a schema error appears; do not reopen it after dismiss."""
function _follow_disk_error!(
    state::BrowserState,
    workspace::Workspace.Workspace,
)::Nothing
    err = workspace.cache.disk_error
    if err isa ProjectCacheSchemaError
        if !state.cache_schema_prompted
            state.cache_rebuild_modal = true
            state.cache_rebuild_error = sprint(showerror, err)
            state.cache_schema_prompted = true
        end
    else
        state.cache_schema_prompted = false
        state.cache_rebuild_modal = false
        state.cache_rebuild_error = ""
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
