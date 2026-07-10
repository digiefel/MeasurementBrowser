import CImGui as ig

using DataBrowserAPI:
    Project,
    display_label,
    kind_label
using DataBrowserCore.ItemIndex:
    HierarchyNode,
    ItemRecord,
    children,
    collection_path_key
import DataBrowserCore.Workspace

"""Destroy text-filter widgets before their saved text is replaced."""
function _reset_project_filter_widgets!(state::BrowserState)::Nothing
    state.tree_filter_widget === nothing || ig.Destroy(state.tree_filter_widget)
    state.item_filter_widget === nothing ||
        ig.Destroy(state.item_filter_widget)
    state.tree_filter_widget = nothing
    state.item_filter_widget = nothing
    return nothing
end

"""Read the current text from a CImGui filter widget."""
function _imgui_text_filter_text(
    filter::Union{Nothing,Ptr{ig.lib.ImGuiTextFilter}},
)::String
    filter === nothing && return ""
    bytes = UInt8[]
    for char in unsafe_load(filter).InputBuf
        char == 0 && break
        push!(bytes, UInt8(mod(Int(char), 256)))
    end
    return String(bytes)
end

"""Return whether an item's visible labels pass the item filter."""
function _item_matches_filter(
    workspace::Workspace.Workspace,
    item::ItemRecord,
    filter_obj::Ptr{ig.lib.ImGuiTextFilter},
)::Bool
    if !ig.ImGuiTextFilter_IsActive(filter_obj)
        return true
    end
    text = string(
        # Label callbacks are project callbacks: they receive the effective parameters, not the
        # item-local subset the record stores.
        display_label(workspace.project,
            Workspace.effective_record(workspace.index.hierarchy, item)),
        "\n",
        item.item_label,
        "\n",
        kind_label(workspace.project, item.kind),
    )
    return ig.ImGuiTextFilter_PassFilter(filter_obj, text, C_NULL)
end

"""Filter items by tags and the text entered in the item panel."""
function _visible_items(
    state::BrowserState,
    workspace::Workspace.Workspace,
    items::Vector{ItemRecord},
    filter_item::Ptr{ig.lib.ImGuiTextFilter},
)::Vector{ItemRecord}
    if _show_bad_effective(state) && !ig.ImGuiTextFilter_IsActive(filter_item)
        return items
    end

    visible = ItemRecord[]
    sizehint!(visible, length(items))
    for item in items
        _item_is_visible(state, item) || continue
        _item_matches_filter(workspace, item, filter_item) || continue
        push!(visible, item)
    end
    return visible
end

"""Render the collection hierarchy and write collection selection into the workspace."""
function _render_hierarchy_tree_panel(
    state::BrowserState,
    filter_tree::Ptr{ig.lib.ImGuiTextFilter},
)::Nothing
    ig.BeginChild("Tree", (0, 0), true)
    ig.SeparatorText("Collection Selection")
    _render_tag_state_error!(state)

    workspace = state.workspace
    root = workspace isa Workspace.Workspace ? workspace.index.hierarchy.root : nothing
    parameter_keys = workspace isa Workspace.Workspace ?
        workspace.index.collection_metadata_keys :
        Symbol[]
    selected_collections, _, _ = _project_visible_selection(state)

    visible_collections = HierarchyNode[]
    visible_collection_keys_ref = Ref{Union{Nothing,Vector{String}}}(nothing)
    expanded_collection_path_set = Set(state.expanded_collection_paths)
    selected_collection_path_set = workspace isa Workspace.Workspace ?
        Set(workspace.selection.collection_paths) :
        Set{String}()
    all_collection_count = Ref(0)
    filter_active = ig.ImGuiTextFilter_IsActive(filter_tree)

    if root !== nothing
        has_visible_leaf_cache = IdDict{HierarchyNode,Bool}()
        subtree_matches_cache = IdDict{HierarchyNode,Bool}()
        collection_key_cache = IdDict{HierarchyNode,String}()

        collection_key(node::HierarchyNode) = get!(collection_key_cache, node) do
            _collection_path_key(node)
        end

        """Return whether this subtree contains a visible collection."""
        function has_visible_leaf(node::HierarchyNode)::Bool
            return get!(has_visible_leaf_cache, node) do
                if isempty(children(node))
                    return _collection_is_visible(state, collection_key(node))
                end
                for child in children(node)
                    has_visible_leaf(child) && return true
                end
                return false
            end
        end

        """Return whether this node or one of its descendants matches the filter."""
        function subtree_matches(node::HierarchyNode)::Bool
            return get!(subtree_matches_cache, node) do
                ig.ImGuiTextFilter_PassFilter(filter_tree, node.name, C_NULL) && return true
                for child in children(node)
                    subtree_matches(child) && return true
                end
                return false
            end
        end

        """Collect visible leaf collections while respecting inherited filter matches."""
        function collect_visible_collections!(
            node::HierarchyNode,
            force_show::Bool=false,
        )::Nothing
            has_visible_leaf(node) || return
            force_show || subtree_matches(node) || return

            direct_match = force_show ||
                ig.ImGuiTextFilter_PassFilter(filter_tree, node.name, C_NULL)
            if isempty(children(node))
                push!(visible_collections, node)
                return nothing
            end
            for child in children(node)
                collect_visible_collections!(child, direct_match)
            end
            return nothing
        end

        """Count leaf collections below one hierarchy node."""
        function count_leaf_nodes(node::HierarchyNode)::Int
            if isempty(children(node))
                return 1
            end
            total = 0
            for child in children(node)
                total += count_leaf_nodes(child)
            end
            return total
        end

        _time!(state, :hierarchy_prep) do
            all_collection_count[] = count_leaf_nodes(root)
            for child in children(root)
                collect_visible_collections!(child, false)
            end
            return nothing
        end

        visible_collection_keys = function ()
            keys = visible_collection_keys_ref[]
            if keys === nothing
                keys = [collection_key(node) for node in visible_collections]
                visible_collection_keys_ref[] = keys
            end
            return keys
        end

        _render_selection_toolbar!(
            length(selected_collections),
            length(visible_collections),
            all_collection_count[],
            filter_tree,
            () -> begin
                workspace.selection.collection_paths = copy(visible_collection_keys())
            end;
            item_label="collections", filter_id="##tree_filter"
        )

        node_seed = UInt64(0x9e3779b97f4a7c15)
        next_node_id(parent_id::UInt64, node::HierarchyNode) = hash(node.name, parent_id)
        to_imgui_id(node_id::UInt64) = Int32(node_id % UInt64(typemax(Int32)))

        """Render one hierarchy row and recursively render its open descendants."""
        function render_node(
            node::HierarchyNode,
            node_id::UInt64,
            path_segments::Vector{String},
            force_show::Bool=false,
        )::Nothing
            has_visible_leaf(node) || return
            force_show || subtree_matches(node) || return

            state.performance.node_count += 1
            ig.TableNextRow()

            direct_match = force_show ||
                ig.ImGuiTextFilter_PassFilter(filter_tree, node.name, C_NULL)
            ig.PushID(to_imgui_id(node_id))

            is_leaf = isempty(children(node))
            path_key = join(path_segments, '/')
            leaf_collection_key = is_leaf ? collection_key(node) : nothing
            selected = is_leaf && (leaf_collection_key in selected_collection_path_set)

            flags = (
                ig.ImGuiTreeNodeFlags_OpenOnArrow |
                ig.ImGuiTreeNodeFlags_OpenOnDoubleClick |
                ig.ImGuiTreeNodeFlags_NavLeftJumpsToParent |
                ig.ImGuiTreeNodeFlags_SpanFullWidth |
                ig.ImGuiTreeNodeFlags_DrawLinesToNodes |
                ig.ImGuiTreeNodeFlags_SpanAllColumns
            )
            if is_leaf
                flags |= (
                    ig.ImGuiTreeNodeFlags_Leaf |
                    ig.ImGuiTreeNodeFlags_Bullet |
                    ig.ImGuiTreeNodeFlags_NoTreePushOnOpen
                )
            end
            if selected
                flags |= ig.ImGuiTreeNodeFlags_Selected
            end
            if filter_active && direct_match && !is_leaf
                flags |= ig.ImGuiTreeNodeFlags_DefaultOpen
            end

            ig.TableSetColumnIndex(0)
            bad_text_pushed = false
            if is_leaf && leaf_collection_key !== nothing && _tag_state_ready(state)
                tag_state = _tag_state_or_error(state)
                ancestor_keys = _ancestor_keys_for_path(leaf_collection_key)
                bad_text_pushed = _push_tag_text_style!(tag_state, leaf_collection_key, ancestor_keys)
            end
            if !is_leaf && !filter_active
                ig.SetNextItemOpen(path_key in expanded_collection_path_set, ig.ImGuiCond_Always)
            end
            opened = ig.TreeNodeEx(is_leaf ? "" : node.name, flags, node.name)
            if is_leaf && state.scroll_to_collection_path == leaf_collection_key
                ig.SetScrollHereY(0.5)
                state.scroll_to_collection_path = nothing
            end
            if !is_leaf && !filter_active && ig.IsItemToggledOpen()
                expanded_collection_paths = copy(state.expanded_collection_paths)
                if opened
                    path_key in expanded_collection_paths || push!(expanded_collection_paths, path_key)
                else
                    filter!(key -> key != path_key, expanded_collection_paths)
                end
                state.expanded_collection_paths = expanded_collection_paths
                expanded_collection_path_set = Set(expanded_collection_paths)
            end
            bad_text_pushed && ig.PopStyleColor()

            if ig.IsItemClicked()
                if is_leaf
                    io = ig.GetIO()
                    shift_held = unsafe_load(io.KeyShift)
                    ctrl_held = unsafe_load(io.KeyCtrl)
                    selected_collection_paths = copy(workspace.selection.collection_paths)
                    _update_multi_selection!(selected_collection_paths, leaf_collection_key, visible_collection_keys(), shift_held, ctrl_held)
                    workspace.selection.collection_paths = selected_collection_paths
                elseif !isempty(workspace.selection.collection_paths)
                    empty!(workspace.selection.collection_paths)
                end
            end

            if is_leaf && ig.BeginPopupContextItem()
                target_nodes = _selection_targets(selected_collections, node)
                target_keys = [collection_key(target) for target in target_nodes]
                selected_count = length(target_keys)
                selected_count > 1 && ig.TextDisabled("Apply to $selected_count collections")

                editable = _tag_state_ready(state)
                if !editable
                    ig.TextDisabled("Fix tags.txt and rescan to edit")
                    ig.Separator()
                end

                !editable && ig.BeginDisabled()
                if ig.MenuItem("Mark Bad")
                    _set_collections_bad!(state, target_keys, true)
                end
                if ig.MenuItem("Unmark Bad")
                    _set_collections_bad!(state, target_keys, false)
                end
                !editable && ig.EndDisabled()
                ig.EndPopup()
            end

            for (i, k) in enumerate(parameter_keys)
                ig.TableSetColumnIndex(i)
                if haskey(node.metadata, k)
                    ig.Text(string(node.metadata[k]))
                elseif is_leaf
                    ig.TextDisabled("--")
                end
            end

            if opened && !is_leaf
                for child in children(node)
                    render_node(child, next_node_id(node_id, child), [path_segments; child.name], direct_match)
                end
                ig.TreePop()
            end
            ig.PopID()
            return nothing
        end

        local table_flags = ig.ImGuiTableFlags_BordersV | ig.ImGuiTableFlags_BordersOuterH |
                            ig.ImGuiTableFlags_Resizable | ig.ImGuiTableFlags_RowBg |
                            ig.ImGuiTableFlags_Reorderable | ig.ImGuiTableFlags_Hideable
        ncols = 1 + length(parameter_keys) + 1
        if ig.BeginTable("tree_table", ncols, table_flags)
            local index_flags = ig.ImGuiTableColumnFlags_NoHide | ig.ImGuiTableColumnFlags_NoReorder |
                                ig.ImGuiTableColumnFlags_NoSort | ig.ImGuiTableColumnFlags_WidthStretch
            ig.TableSetupColumn("Collection", index_flags, 5.0)
            for k in parameter_keys
                ig.TableSetupColumn(String(k), ig.ImGuiTableColumnFlags_AngledHeader | ig.ImGuiTableFlags_SizingFixedFit)
            end
            ig.TableSetupColumn("")
            ig.TableAngledHeadersRow()
            ig.TableHeadersRow()
            for child in children(root)
                render_node(child, next_node_id(node_seed, child), [child.name], false)
            end
            ig.EndTable()
        end
    else
        _render_selection_toolbar!(
            0,
            0,
            0,
            filter_tree,
            () -> nothing;
            item_label="collections", filter_id="##tree_filter"
        )
        ig.Text("No data loaded")
    end

    isempty(visible_collections) && all_collection_count[] > 0 && ig.TextDisabled("No collections match filter")
    ig.EndChild()
end

"""Render selection counts shared by the collection and item panels."""
function _render_selection_status!(
    selected_count::Int,
    filtered_count::Int,
    total_count::Int,
    item_type::String,
)::Nothing
    if selected_count > 1
        ig.TextColored((0.2, 0.8, 0.2, 1.0), "$selected_count/$filtered_count $item_type selected")
    elseif selected_count == 1
        ig.TextColored((0.6, 0.8, 1.0, 1.0), "1/$filtered_count $item_type selected")
    else
        ig.TextDisabled("0/$filtered_count $item_type selected")
    end

    if filtered_count != total_count
        ig.SameLine()
        ig.TextDisabled("($total_count total)")
    end
    ig.SameLine()
    _helpmarker("Multi-select: Shift+click=range, Ctrl+click=toggle, Ctrl+A=all")
    return nothing
end

"""Render filtering and select-all controls shared by both selection panels."""
function _render_selection_toolbar!(
    selected_count::Int,
    visible_count::Int,
    total_count::Int,
    filter_obj::Ptr{ig.lib.ImGuiTextFilter},
    on_select_all!::Function;
    item_label::String,
    filter_id::String,
)::Nothing
    _render_selection_status!(selected_count, visible_count, total_count, item_label)

    ig.Text("Filter")
    ig.SameLine()
    _helpmarker("incl,-excl")
    ig.SameLine()
    ig.SetNextItemShortcut(
        ig.ImGuiMod_Ctrl | ig.ImGuiKey_F,
        ig.ImGuiInputFlags_Tooltip
    )
    ig.ImGuiTextFilter_Draw(filter_obj, filter_id, -1)
    if ig.IsKeyPressed(ig.ImGuiKey_A) && ig.IsWindowFocused()
        io = ig.GetIO()
        unsafe_load(io.KeyCtrl) && on_select_all!()
    end
    return nothing
end

# Right panel (items list) rendering
"""Render items for the selected collections and update workspace selection."""
function _render_items_panel(
    state::BrowserState,
    filter_item::Ptr{ig.lib.ImGuiTextFilter},
)::Nothing
    ig.BeginChild("Items", (0, 0), true)
    ig.SeparatorText("Item Selection")
    workspace = state.workspace
    if !(workspace isa Workspace.Workspace)
        ig.TextDisabled("Open a project folder to browse items")
        ig.EndChild()
        return nothing
    end

    selected_collections, selected_items, selected_path =
        _project_visible_selection(state)
    all_items_ref = Ref{Vector{ItemRecord}}()
    _time!(state, :items_panel) do
        all_items_ref[] = _items_of_selected_collections(state)
        return nothing
    end
    all_items = all_items_ref[]
    visible_items_ref = Ref{Vector{ItemRecord}}()
    _time!(state, :visible_items) do
        visible_items_ref[] = _visible_items(state, workspace, all_items, filter_item)
        return nothing
    end
    visible_items = visible_items_ref[]
    selected_id_set = Set(workspace.selection.item_ids)
    registry_ready = _tag_state_ready(state)

    _render_selection_toolbar!(
        length(selected_items),
        length(visible_items),
        length(all_items),
        filter_item,
        () -> begin
            workspace.selection.item_ids = [item.id for item in visible_items]
        end;
        item_label="items", filter_id="##items_filter"
    )

    if !isempty(selected_collections)
        if length(selected_collections) == 1
            sel_name = join(selected_path, "/")
            ig.Text("Items for $sel_name")
            ig.Separator()
        else
            shown = min(3, length(selected_collections))
            first_names = join((selected_collections[i].name for i in 1:shown), ", ")
            ig.Text("Items from $(length(selected_collections)) collections: $first_names$(length(selected_collections) > 3 ? "..." : "")")
            ig.Separator()
        end
    end

    if isempty(selected_collections)
        ig.Text("Select one or more collections to view their items")
        ig.EndChild()
        return nothing
    end

    state.performance.item_rows_visible = length(visible_items)
    state.performance.item_rows_rendered = 0

    if !isempty(visible_items)
        table_flags = ig.ImGuiTableFlags_BordersOuterH | ig.ImGuiTableFlags_BordersOuterV |
                      ig.ImGuiTableFlags_RowBg | ig.ImGuiTableFlags_SizingStretchProp
        if ig.BeginTable("items_table", 1, table_flags)
            ig.TableSetupColumn("Item", ig.ImGuiTableColumnFlags_WidthStretch)

            rows_rendered = 0
            clipper = ig.lib.ImGuiListClipper()
            try
                ig.Begin(clipper, length(visible_items), ig.GetTextLineHeightWithSpacing())
                while ig.Step(clipper)
                    start_idx = Int(unsafe_load(clipper.DisplayStart)) + 1
                    end_idx = Int(unsafe_load(clipper.DisplayEnd))
                    for idx in start_idx:end_idx
                        item = visible_items[idx]
                        id = item.id
                        ig.PushID(id)
                        ig.TableNextRow()
                        ig.TableSetColumnIndex(0)
                        rows_rendered += 1

                        is_selected = id in selected_id_set
                        bad_text_pushed = false
                        if registry_ready
                            tag_state = _tag_state_or_error(state)
                            dev_key = collection_path_key(item.collection)
                            ancestor_keys = _ancestor_keys_for_path(dev_key)
                            bad_text_pushed = _push_tag_text_style!(
                                tag_state, id, [dev_key; ancestor_keys])
                        end
                        if ig.Selectable(
                            display_label(workspace.project,
                                Workspace.effective_record(workspace.index.hierarchy, item)),
                            is_selected,
                            ig.ImGuiSelectableFlags_SpanAllColumns,
                        )
                            io = ig.GetIO()
                            shift_held = unsafe_load(io.KeyShift)
                            ctrl_held = unsafe_load(io.KeyCtrl)
                            selected_ids = copy(workspace.selection.item_ids)
                            visible_ids = [visible.id for visible in visible_items]
                            _update_multi_selection!(
                                selected_ids,
                                id,
                                visible_ids,
                                shift_held,
                                ctrl_held,
                            )
                            workspace.selection.item_ids = selected_ids
                        end
                        if state.scroll_to_item_id == id
                            ig.SetScrollHereY(0.5)
                            state.scroll_to_item_id = nothing
                        end
                        bad_text_pushed && ig.PopStyleColor()

                        if ig.BeginPopupContextItem()
                            target_items = _selection_targets(selected_items, item)
                            target_ids = [target.id for target in target_items]
                            selected_count = length(target_ids)
                            selected_count > 1 && ig.TextDisabled("Apply to $selected_count items")

                            if ig.MenuItem("Open Plot in New Window")
                                for ext in state.extensions
                                    Browser.open_detached_plot!(
                                        ext,
                                        state;
                                        item_ids=[id],
                                        title=item.item_label,
                                        kind=item.kind,
                                    ) && break
                                end
                            end

                            editable = _tag_state_ready(state)
                            if !editable
                                ig.TextDisabled("Fix tags.txt and rescan to edit")
                                ig.Separator()
                            end

                            !editable && ig.BeginDisabled()
                            if ig.MenuItem("Mark Bad")
                                _set_items_bad!(state, target_ids, true)
                            end
                            if ig.MenuItem("Unmark Bad")
                                _set_items_bad!(state, target_ids, false)
                            end
                            !editable && ig.EndDisabled()
                            ig.EndPopup()
                        end
                        ig.PopID()
                    end
                end
            finally
                ig.Destroy(clipper)
            end

            state.performance.item_rows_rendered = rows_rendered
            state.performance.item_rows_visible = length(visible_items)
            ig.EndTable()
        end
    else
        ig.TextDisabled("No items match filter")
    end
    ig.EndChild()
    return nothing
end

"""Render both selection panels and keep their text filters in browser state."""
function render_selection_window(state::BrowserState)::Nothing
    state.performance.node_count = 0
    if state.reset_project_filters
        _reset_project_filter_widgets!(state)
        state.reset_project_filters = false
    end
    if state.tree_filter_widget === nothing
        state.tree_filter_widget =
            ig.ImGuiTextFilter_ImGuiTextFilter(state.tree_filter)
    end
    if state.item_filter_widget === nothing
        state.item_filter_widget =
            ig.ImGuiTextFilter_ImGuiTextFilter(state.item_filter)
    end
    filter_tree = state.tree_filter_widget
    filter_item = state.item_filter_widget

    if ig.Begin("Hierarchy", C_NULL, ig.ImGuiWindowFlags_MenuBar)
        render_menu_bar(state)
        ig.Columns(2, "main_layout")
        _time!(state, :collection_tree) do
            _render_hierarchy_tree_panel(state, filter_tree)
        end
        ig.NextColumn()
        _time!(state, :item_panel) do
            _render_items_panel(state, filter_item)
        end
    end
    ig.End()
    state.tree_filter = _imgui_text_filter_text(filter_tree)
    state.item_filter = _imgui_text_filter_text(filter_item)
    return nothing
end
