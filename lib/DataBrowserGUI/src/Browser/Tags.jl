using DataBrowserAnnotations

import DataBrowserCore.Workspace
using DataBrowserAPI.ItemIndex:
    CollectionRecord,
    ItemRecord,
    collection_item_ids,
    collection_location,
    collection_path_keys,
    has_collection_children,
    item_timestamp_key

const BAD_TAG_NAME = "bad"
const BAD_TAG_COLOR = (UInt8(0xff), UInt8(0x30), UInt8(0x30))
const BAD_TAG_PRIORITY = 100
const ROOT_COLLECTION_SELECTION_ID = ""

"""
Load tags for one source root.

Parse or I/O failures leave tags unavailable but keep bad items visible.
"""
function _annotation_root(workspace::Workspace.Workspace)::String
    root = dirname(workspace.cache.identity.cache_path)
    mkpath(root)
    return root
end

function _load_tag_state_for_root!(
    state::BrowserState,
    root_path::String,
)::Nothing
    if isempty(root_path)
        state.tag_state = nothing
        state.tag_state_error = ""
        return nothing
    end

    try
        state.tag_state = DataBrowserAnnotations.Tags.load(root_path)
        state.tag_state_error = ""
    catch err
        if err isa DataBrowserAnnotations.Tags.TagsParseError || err isa IOError
            state.tag_state = nothing
            state.tag_state_error = sprint(showerror, err)
            return nothing
        end
        rethrow()
    end
    return nothing
end

function _tag_state_ready(state::BrowserState)::Bool
    return state.tag_state !== nothing && isempty(state.tag_state_error)
end

"""
Return whether bad-tagged items should be visible in the current UI state.
Bad items stay visible when tags cannot be loaded, so missing tag metadata never hides data.
"""
function _show_bad_effective(state::BrowserState)::Bool
    return state.show_bad || !_tag_state_ready(state)
end

"""Return loaded tags or fail with the error that prevented loading them."""
function _tag_state_or_error(state::BrowserState)::DataBrowserAnnotations.Tags.TagState
    state.tag_state === nothing &&
        error("Tag state unavailable: $(state.tag_state_error)")
    return state.tag_state
end

"""Return the persisted annotation key for one collection record."""
function _collection_annotation_key(
    collections,
    collection_record::CollectionRecord,
)::String
    path = CollectionRecord[
        collections.records[key]
        for key in collection_path_keys(collections, collection_record.key)
    ]
    if all(segment -> segment.registration_name !== nothing, path)
        return join((segment.registration_name::String for segment in path), '/')
    end
    return "@collection/$(collection_record.id)"
end

"""Return every parent annotation key used when resolving inherited tags."""
function _ancestor_annotation_keys(
    collections,
    collection_record::CollectionRecord,
)::Vector{String}
    path = collection_path_keys(collections, collection_record.key)
    return String[
        _collection_annotation_key(collections, collections.records[key])
        for key in path[1:end-1]
    ]
end

"""Return whether an item or one of its supplied parents has the bad tag."""
function _has_bad_tag(
    tag_state::DataBrowserAnnotations.Tags.TagState,
    key::AbstractString,
    ancestor_keys::Vector{String},
)::Bool
    return "bad" in DataBrowserAnnotations.Tags.effective(tag_state, key, ancestor_keys)
end

"""Return whether a collection remains visible after applying the bad-tag filter."""
function _collection_is_visible(
    state::BrowserState,
    collection_record::CollectionRecord,
)::Bool
    _show_bad_effective(state) && return true
    tag_state = _tag_state_or_error(state)
    collections = (state.workspace::Workspace.Workspace).index.collections
    collection_key = _collection_annotation_key(collections, collection_record)
    return !_has_bad_tag(
        tag_state,
        collection_key,
        _ancestor_annotation_keys(collections, collection_record),
    )
end

"""Return whether an item remains visible after applying inherited bad tags."""
function _item_is_visible(
    state::BrowserState,
    item::ItemRecord,
)::Bool
    _show_bad_effective(state) && return true
    tag_state = _tag_state_or_error(state)
    collections = (state.workspace::Workspace.Workspace).index.collections
    item.collection_key === nothing &&
        return !_has_bad_tag(tag_state, item.id, String[])
    collection_record = collections.records[item.collection_key]
    collection_key = _collection_annotation_key(collections, collection_record)
    ancestors = [
        collection_key;
        _ancestor_annotation_keys(collections, collection_record)
    ]
    return !_has_bad_tag(tag_state, item.id, ancestors)
end

"""
Resolve the current selected collections and items after applying tag visibility.
The persisted selection ids are left untouched; this returns only the currently visible subset.
"""
function _project_visible_selection(
    state::BrowserState,
)::Tuple{Vector{CollectionRecord},Vector{ItemRecord},Vector{String}}
    workspace = state.workspace
    if !(workspace isa Workspace.Workspace)
        return CollectionRecord[], ItemRecord[], String[]
    end
    collections = workspace.index.collections

    selected_collections = CollectionRecord[]
    for collection_id in workspace.selection.collection_ids
        key = get(collections.key_by_id, collection_id, nothing)
        key === nothing && continue
        collection_record = collections.records[key]
        !has_collection_children(collections, key) ||
            error("Selected collection ID '$collection_id' does not point to a leaf collection")
        !_collection_is_visible(state, collection_record) && continue
        push!(selected_collections, collection_record)
    end

    visible_collection_ids = Set(collection_record.id for collection_record in selected_collections)
    root_selected = ROOT_COLLECTION_SELECTION_ID in workspace.selection.collection_ids
    item_index = workspace.index.items
    selected_items = ItemRecord[]
    for id in workspace.selection.item_ids
        item = get(item_index, id, nothing)
        item === nothing && continue
        if item.collection_key === nothing
            root_selected || continue
        else
            collections.records[item.collection_key].id in visible_collection_ids || continue
        end
        !_item_is_visible(state, item) && continue
        push!(selected_items, item)
    end
    selected_path = length(selected_collections) == 1 ?
        collection_location(collections, selected_collections[1].key) : String[]
    return selected_collections, selected_items, selected_path
end

"""
Return every item belonging to the selected visible collections, ordered by time.
"""
function _items_of_selected_collections(
    state::BrowserState,
)::Vector{ItemRecord}
    selected_collections, _, _ = _project_visible_selection(state)
    all_items = ItemRecord[]
    collections = (state.workspace::Workspace.Workspace).index.collections
    items = (state.workspace::Workspace.Workspace).index.items
    sizehint!(all_items, sum(
        length(collection_item_ids(collections, collection_record.key))
        for collection_record in selected_collections; init=0) +
        (ROOT_COLLECTION_SELECTION_ID in
         (state.workspace::Workspace.Workspace).selection.collection_ids ?
            length(collection_item_ids(collections, nothing)) : 0))
    if ROOT_COLLECTION_SELECTION_ID in
       (state.workspace::Workspace.Workspace).selection.collection_ids
        for id in collection_item_ids(collections, nothing)
            haskey(items, id) && push!(all_items, items[id])
        end
    end
    for collection_record in selected_collections
        for id in collection_item_ids(collections, collection_record.key)
            haskey(items, id) && push!(all_items, items[id])
        end
    end
    sort!(all_items, by=item_timestamp_key)
    return all_items
end

"""
Ensure the bad tag exists before writing bad-collection or bad-item assignments.
"""
function _ensure_bad_catalog_entry!(
    tag_state::DataBrowserAnnotations.Tags.TagState,
)::Nothing
    any(t -> t.name == BAD_TAG_NAME, tag_state.catalog) && return
    pushfirst!(tag_state.catalog,
        DataBrowserAnnotations.Tags.TagDef(BAD_TAG_NAME, BAD_TAG_COLOR, BAD_TAG_PRIORITY))
    return nothing
end

"""
Set or clear the bad tag on collections and persist the updated tag file.
Returns `false` when tag state is unavailable or there is nothing to change.
"""
function _set_collections_bad!(
    state::BrowserState,
    annotation_keys::Vector{String},
    bad::Bool,
)::Bool
    unique_annotation_keys = unique(copy(annotation_keys))
    isempty(unique_annotation_keys) && return false
    _tag_state_ready(state) || return false

    workspace = state.workspace::Workspace.Workspace

    tag_state = _tag_state_or_error(state)
    _ensure_bad_catalog_entry!(tag_state)

    for annotation_key in unique_annotation_keys
        if bad
            set = get!(() -> Set{String}(), tag_state.assignments, annotation_key)
            push!(set, BAD_TAG_NAME)
        else
            tags = get(tag_state.assignments, annotation_key, nothing)
            if tags !== nothing
                delete!(tags, BAD_TAG_NAME)
                isempty(tags) && delete!(tag_state.assignments, annotation_key)
            end
        end
    end

    DataBrowserAnnotations.Tags.save(_annotation_root(workspace), tag_state)
    return true
end

"""
Set or clear the bad tag on items and persist the updated tag file.
Returns `false` when tag state is unavailable or there is nothing to change.
"""
function _set_items_bad!(
    state::BrowserState,
    item_ids::Vector{String},
    bad::Bool,
)::Bool
    unique_ids = unique(copy(item_ids))
    isempty(unique_ids) && return false
    _tag_state_ready(state) || return false

    workspace = state.workspace::Workspace.Workspace

    tag_state = _tag_state_or_error(state)
    _ensure_bad_catalog_entry!(tag_state)

    for id in unique_ids
        if bad
            set = get!(() -> Set{String}(), tag_state.assignments, id)
            push!(set, BAD_TAG_NAME)
        else
            tags = get(tag_state.assignments, id, nothing)
            if tags !== nothing
                delete!(tags, BAD_TAG_NAME)
                isempty(tags) && delete!(tag_state.assignments, id)
            end
        end
    end

    DataBrowserAnnotations.Tags.save(_annotation_root(workspace), tag_state)
    return true
end

"""Apply a context-menu action to the selection when the clicked item belongs to it."""
function _selection_targets(
    selected_items::Vector{T},
    clicked_item::T,
)::Vector{T} where {T}
    if clicked_item in selected_items
        return selected_items
    end
    return T[clicked_item]
end
