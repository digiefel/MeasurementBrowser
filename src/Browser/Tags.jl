using Annotations

import ..Workspace
using ..ItemIndex:
    HierarchyNode,
    ItemRecord,
    collection_path_key,
    collection_path_tuple,
    item_record_key,
    item_timestamp_key

const BAD_TAG_NAME = "bad"
const BAD_TAG_COLOR = (UInt8(0xff), UInt8(0x30), UInt8(0x30))
const BAD_TAG_PRIORITY = 100

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
        state.tag_state = Annotations.Tags.load(root_path)
        state.tag_state_error = ""
    catch err
        if err isa Annotations.Tags.TagsParseError || err isa IOError
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
function _tag_state_or_error(state::BrowserState)::Annotations.Tags.TagState
    state.tag_state === nothing &&
        error("Tag state unavailable: $(state.tag_state_error)")
    return state.tag_state
end

"""Return the stable collection path represented by a leaf hierarchy node."""
function _collection_path_key(node::HierarchyNode)::String
    isempty(node.items) && error("Leaf node '$(node.name)' has no items")
    return collection_path_key(first(node.items).collection)
end

"""Return the collection path segments represented by a leaf hierarchy node."""
function _collection_location(node::HierarchyNode)::Vector{String}
    isempty(node.items) && error("Leaf node '$(node.name)' has no items")
    return copy(first(node.items).collection)
end

"""Return every parent path used when resolving inherited tags."""
function _ancestor_keys_for_path(path::AbstractString)::Vector{String}
    parts = split(path, '/')
    length(parts) <= 1 && return String[]
    out = String[]
    for i in 1:(length(parts) - 1)
        push!(out, join(parts[1:i], '/'))
    end
    return out
end

"""Return whether an item or one of its supplied parents has the bad tag."""
function _has_bad_tag(
    tag_state::Annotations.Tags.TagState,
    key::AbstractString,
    ancestor_keys::Vector{String},
)::Bool
    return "bad" in Annotations.Tags.effective(tag_state, key, ancestor_keys)
end

"""Return whether a collection remains visible after applying the bad-tag filter."""
function _collection_is_visible(state::BrowserState, collection_key::String)::Bool
    _show_bad_effective(state) && return true
    tag_state = _tag_state_or_error(state)
    return !_has_bad_tag(
        tag_state,
        collection_key,
        _ancestor_keys_for_path(collection_key),
    )
end

"""Return whether an item remains visible after applying inherited bad tags."""
function _item_is_visible(
    state::BrowserState,
    item::ItemRecord,
)::Bool
    _show_bad_effective(state) && return true
    tag_state = _tag_state_or_error(state)
    collection_key = collection_path_key(item.collection)
    ancestors = [collection_key; _ancestor_keys_for_path(collection_key)]
    return !_has_bad_tag(tag_state, item_record_key(item), ancestors)
end

"""
Resolve the current selected collections and items after applying tag visibility.
The persisted selection ids are left untouched; this returns only the currently visible subset.
"""
function _project_visible_selection(
    state::BrowserState,
)::Tuple{Vector{HierarchyNode},Vector{ItemRecord},Vector{String}}
    workspace = state.workspace
    if !(workspace isa Workspace.Workspace)
        return HierarchyNode[], ItemRecord[], String[]
    end
    hierarchy = workspace.index.hierarchy

    selected_collections = HierarchyNode[]
    for path_key in workspace.selection.collection_paths
        node = get(hierarchy.index, collection_path_tuple(path_key), nothing)
        node === nothing && continue
        node.kind == :leaf ||
            error("Selected collection path '$path_key' does not point to a leaf collection")
        !_collection_is_visible(state, path_key) && continue
        push!(selected_collections, node)
    end

    visible_collection_keys = Set(_collection_path_key(node) for node in selected_collections)
    item_index = workspace.index.items
    selected_items = ItemRecord[]
    for item_key in workspace.selection.item_keys
        item = get(item_index, item_key, nothing)
        item === nothing && continue
        collection_path_key(item.collection) in visible_collection_keys || continue
        !_item_is_visible(state, item) && continue
        push!(selected_items, item)
    end
    selected_path = length(selected_collections) == 1 ? _collection_location(selected_collections[1]) : String[]
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
    sizehint!(all_items, sum(length(collection.items) for collection in selected_collections; init=0))
    for collection in selected_collections
        for item in collection.items
            push!(all_items, item)
        end
    end
    sort!(all_items, by=item_timestamp_key)
    return all_items
end

"""
Ensure the bad tag exists before writing bad-collection or bad-item assignments.
"""
function _ensure_bad_catalog_entry!(
    tag_state::Annotations.Tags.TagState,
)::Nothing
    any(t -> t.name == BAD_TAG_NAME, tag_state.catalog) && return
    pushfirst!(tag_state.catalog,
        Annotations.Tags.TagDef(BAD_TAG_NAME, BAD_TAG_COLOR, BAD_TAG_PRIORITY))
    return nothing
end

"""
Set or clear the bad tag on collections and persist the updated tag file.
Returns `false` when tag state is unavailable or there is nothing to change.
"""
function _set_collections_bad!(
    state::BrowserState,
    collection_keys::Vector{String},
    bad::Bool,
)::Bool
    unique_keys = unique(copy(collection_keys))
    isempty(unique_keys) && return false
    _tag_state_ready(state) || return false

    workspace = state.workspace::Workspace.Workspace

    tag_state = _tag_state_or_error(state)
    _ensure_bad_catalog_entry!(tag_state)

    for collection_key in unique_keys
        if bad
            set = get!(() -> Set{String}(), tag_state.assignments, collection_key)
            push!(set, BAD_TAG_NAME)
        else
            tags = get(tag_state.assignments, collection_key, nothing)
            if tags !== nothing
                delete!(tags, BAD_TAG_NAME)
                isempty(tags) && delete!(tag_state.assignments, collection_key)
            end
        end
    end

    Annotations.Tags.save(_annotation_root(workspace), tag_state)
    return true
end

"""
Set or clear the bad tag on items and persist the updated tag file.
Returns `false` when tag state is unavailable or there is nothing to change.
"""
function _set_items_bad!(
    state::BrowserState,
    item_keys::Vector{String},
    bad::Bool,
)::Bool
    unique_keys = unique(copy(item_keys))
    isempty(unique_keys) && return false
    _tag_state_ready(state) || return false

    workspace = state.workspace::Workspace.Workspace

    tag_state = _tag_state_or_error(state)
    _ensure_bad_catalog_entry!(tag_state)

    for item_key in unique_keys
        if bad
            set = get!(() -> Set{String}(), tag_state.assignments, item_key)
            push!(set, BAD_TAG_NAME)
        else
            tags = get(tag_state.assignments, item_key, nothing)
            if tags !== nothing
                delete!(tags, BAD_TAG_NAME)
                isempty(tags) && delete!(tag_state.assignments, item_key)
            end
        end
    end

    Annotations.Tags.save(_annotation_root(workspace), tag_state)
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
