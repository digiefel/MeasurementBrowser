import CImGui as ig

using ..Projects:
    Project,
    display_label,
    kind_label
using ..MeasurementIndex:
    HierarchyNode,
    MeasurementInfo,
    children,
    device_path_key
import ..Workspace

"""Destroy text-filter widgets before their saved text is replaced."""
function _reset_project_filter_widgets!(state::BrowserState)::Nothing
    state.tree_filter_widget === nothing || ig.Destroy(state.tree_filter_widget)
    state.measurement_filter_widget === nothing ||
        ig.Destroy(state.measurement_filter_widget)
    state.tree_filter_widget = nothing
    state.measurement_filter_widget = nothing
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

# Left panel (hierarchy tree) rendering
"""Return whether a hierarchy label passes the tree filter."""
function _tree_node_matches_filter(
    filter_tree::Ptr{ig.lib.ImGuiTextFilter},
    node::HierarchyNode,
)::Bool
    return ig.ImGuiTextFilter_PassFilter(filter_tree, node.name, C_NULL)
end

"""Return whether a measurement's visible labels pass the measurement filter."""
function _measurement_matches_filter(
    project::Project,
    measurement::MeasurementInfo,
    filter_obj::Ptr{ig.lib.ImGuiTextFilter},
)::Bool
    if !ig.ImGuiTextFilter_IsActive(filter_obj)
        return true
    end
    text = string(
        display_label(project, measurement),
        "\n",
        measurement.clean_title,
        "\n",
        kind_label(project, measurement.measurement_kind),
    )
    return ig.ImGuiTextFilter_PassFilter(filter_obj, text, C_NULL)
end

"""Filter measurements by tags and the text entered in the measurement panel."""
function _visible_measurements(
    state::BrowserState,
    project::Project,
    measurements::Vector{MeasurementInfo},
    filter_meas::Ptr{ig.lib.ImGuiTextFilter},
)::Vector{MeasurementInfo}
    if _show_bad_effective(state) && !ig.ImGuiTextFilter_IsActive(filter_meas)
        return measurements
    end

    visible = MeasurementInfo[]
    sizehint!(visible, length(measurements))
    for measurement in measurements
        _measurement_is_visible(state, measurement) || continue
        _measurement_matches_filter(project, measurement, filter_meas) || continue
        push!(visible, measurement)
    end
    return visible
end

"""Render the device hierarchy and write device selection into the workspace."""
function _render_hierarchy_tree_panel(
    state::BrowserState,
    filter_tree::Ptr{ig.lib.ImGuiTextFilter},
)::Nothing
    ig.BeginChild("Tree", (0, 0), true)
    ig.SeparatorText("Device Selection")
    _render_tag_state_error!(state)

    workspace = state.workspace
    root = workspace isa Workspace.Workspace ? workspace.index.hierarchy.root : nothing
    meta_keys = workspace isa Workspace.Workspace ?
        workspace.index.device_metadata_keys :
        Symbol[]
    selected_devices, _, _ = _project_visible_selection(state)

    visible_devices = HierarchyNode[]
    visible_device_keys_ref = Ref{Union{Nothing,Vector{String}}}(nothing)
    expanded_device_path_set = Set(state.expanded_device_paths)
    selected_device_path_set = workspace isa Workspace.Workspace ?
        Set(workspace.selection.device_paths) :
        Set{String}()
    all_device_count = 0
    filter_active = ig.ImGuiTextFilter_IsActive(filter_tree)

    if root !== nothing
        has_visible_leaf_cache = IdDict{HierarchyNode,Bool}()
        subtree_matches_cache = IdDict{HierarchyNode,Bool}()
        device_key_cache = IdDict{HierarchyNode,String}()

        device_key(node::HierarchyNode) = get!(device_key_cache, node) do
            _device_path_key(node)
        end

        """Return whether this subtree contains a visible device."""
        function has_visible_leaf(node::HierarchyNode)::Bool
            return get!(has_visible_leaf_cache, node) do
                if isempty(children(node))
                    return _device_is_visible(state, device_key(node))
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
                _tree_node_matches_filter(filter_tree, node) && return true
                for child in children(node)
                    subtree_matches(child) && return true
                end
                return false
            end
        end

        """Collect visible leaf devices while respecting inherited filter matches."""
        function collect_visible_devices!(
            node::HierarchyNode,
            force_show::Bool=false,
        )::Nothing
            has_visible_leaf(node) || return
            force_show || subtree_matches(node) || return

            direct_match = force_show || _tree_node_matches_filter(filter_tree, node)
            if isempty(children(node))
                push!(visible_devices, node)
                return nothing
            end
            for child in children(node)
                collect_visible_devices!(child, direct_match)
            end
            return nothing
        end

        """Count device leaves below one hierarchy node."""
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

        all_device_count = count_leaf_nodes(root)
        for child in children(root)
            collect_visible_devices!(child, false)
        end

        visible_device_keys = function ()
            keys = visible_device_keys_ref[]
            if keys === nothing
                keys = [device_key(node) for node in visible_devices]
                visible_device_keys_ref[] = keys
            end
            return keys
        end

        _render_selection_toolbar!(
            length(selected_devices),
            length(visible_devices),
            all_device_count,
            filter_tree,
            () -> begin
                workspace.selection.device_paths = copy(visible_device_keys())
            end;
            item_label="devices", filter_id="##tree_filter"
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

            direct_match = force_show || _tree_node_matches_filter(filter_tree, node)
            ig.PushID(to_imgui_id(node_id))

            is_leaf = isempty(children(node))
            path_key = join(path_segments, '/')
            leaf_device_key = is_leaf ? device_key(node) : nothing
            selected = is_leaf && (leaf_device_key in selected_device_path_set)

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
            if is_leaf && leaf_device_key !== nothing && _tag_state_ready(state)
                tag_state = _tag_state_or_error(state)
                ancestor_keys = _ancestor_keys_for_path(leaf_device_key)
                bad_text_pushed = _push_tag_text_style!(tag_state, leaf_device_key, ancestor_keys)
            end
            if !is_leaf && !filter_active
                ig.SetNextItemOpen(path_key in expanded_device_path_set, ig.ImGuiCond_Always)
            end
            opened = ig.TreeNodeEx(is_leaf ? "" : node.name, flags, node.name)
            if is_leaf && state.scroll_to_device_path == leaf_device_key
                ig.SetScrollHereY(0.5)
                state.scroll_to_device_path = nothing
            end
            if !is_leaf && !filter_active && ig.IsItemToggledOpen()
                expanded_device_paths = copy(state.expanded_device_paths)
                if opened
                    path_key in expanded_device_paths || push!(expanded_device_paths, path_key)
                else
                    filter!(key -> key != path_key, expanded_device_paths)
                end
                state.expanded_device_paths = expanded_device_paths
                expanded_device_path_set = Set(expanded_device_paths)
            end
            bad_text_pushed && ig.PopStyleColor()

            if ig.IsItemClicked()
                if is_leaf
                    io = ig.GetIO()
                    shift_held = unsafe_load(io.KeyShift)
                    ctrl_held = unsafe_load(io.KeyCtrl)
                    selected_device_paths = copy(workspace.selection.device_paths)
                    _update_multi_selection!(selected_device_paths, leaf_device_key, visible_device_keys(), shift_held, ctrl_held)
                    workspace.selection.device_paths = selected_device_paths
                elseif !isempty(workspace.selection.device_paths)
                    empty!(workspace.selection.device_paths)
                end
            end

            if is_leaf && ig.BeginPopupContextItem()
                target_nodes = _selection_targets(selected_devices, node)
                target_keys = [device_key(target) for target in target_nodes]
                selected_count = length(target_keys)
                selected_count > 1 && ig.TextDisabled("Apply to $selected_count devices")

                editable = _tag_state_ready(state)
                if !editable
                    ig.TextDisabled("Fix tags.txt and rescan to edit")
                    ig.Separator()
                end

                !editable && ig.BeginDisabled()
                if ig.MenuItem("Mark Bad")
                    _set_devices_bad!(state, target_keys, true)
                end
                if ig.MenuItem("Unmark Bad")
                    _set_devices_bad!(state, target_keys, false)
                end
                !editable && ig.EndDisabled()
                ig.EndPopup()
            end

            dev_meta = nothing
            if is_leaf && !isempty(node.measurements)
                dev_meta = first(node.measurements).device_info.parameters
            end
            for (i, k) in enumerate(meta_keys)
                ig.TableSetColumnIndex(i)
                if dev_meta !== nothing && haskey(dev_meta, k)
                    ig.Text(string(dev_meta[k]))
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
        ncols = 1 + length(meta_keys) + 1
        if ig.BeginTable("tree_table", ncols, table_flags)
            local index_flags = ig.ImGuiTableColumnFlags_NoHide | ig.ImGuiTableColumnFlags_NoReorder |
                                ig.ImGuiTableColumnFlags_NoSort | ig.ImGuiTableColumnFlags_WidthStretch
            ig.TableSetupColumn("Device", index_flags, 5.0)
            for k in meta_keys
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
            item_label="devices", filter_id="##tree_filter"
        )
        ig.Text("No data loaded")
    end

    isempty(visible_devices) && all_device_count > 0 && ig.TextDisabled("No devices match filter")
    ig.EndChild()
end

# Shared multi-select utility functions

"""Apply range, toggle, or replacement selection to a list of stable item ids."""
function _update_multi_selection!(
    selected_items::Vector{T},
    item::T,
    all_items::Vector{T},
    shift_held::Bool,
    ctrl_held::Bool,
)::Nothing where {T}
    if shift_held && !isempty(selected_items)
        # Range selection: select from last selected to current
        last_item = selected_items[end]
        start_idx = findfirst(x -> x == last_item, all_items)
        end_idx = findfirst(x -> x == item, all_items)

        if start_idx !== nothing && end_idx !== nothing
            if start_idx > end_idx
                start_idx, end_idx = end_idx, start_idx
            end
            selected_range = all_items[start_idx:end_idx]
            # Merge with existing selection
            append!(selected_items, selected_range)
            unique!(selected_items)
        end
    elseif ctrl_held
        # Toggle selection
        if item in selected_items
            filter!(x -> x != item, selected_items)
        else
            push!(selected_items, item)
        end
    else
        # Single selection (replace existing)
        empty!(selected_items)
        push!(selected_items, item)
    end
    return nothing
end

"""Render selection counts shared by the device and measurement panels."""
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

# Right panel (measurements list) rendering
"""Render measurements for the selected devices and update workspace selection."""
function _render_measurements_panel(
    state::BrowserState,
    filter_meas::Ptr{ig.lib.ImGuiTextFilter},
)::Nothing
    ig.BeginChild("Measurements", (0, 0), true)
    ig.SeparatorText("Measurement Selection")
    workspace = state.workspace
    if !(workspace isa Workspace.Workspace)
        ig.TextDisabled("Open a project folder to browse measurements")
        ig.EndChild()
        return nothing
    end
    proj = workspace.project

    selected_devices, selected_measurements, selected_path =
        _project_visible_selection(state)
    all_measurements = _measurements_of_selected_devices(state)
    visible_measurements =
        _visible_measurements(state, proj, all_measurements, filter_meas)
    selected_measurement_id_set = Set(workspace.selection.measurement_ids)
    registry_ready = _tag_state_ready(state)

    _render_selection_toolbar!(
        length(selected_measurements),
        length(visible_measurements),
        length(all_measurements),
        filter_meas,
        () -> begin
            workspace.selection.measurement_ids =
                [measurement.unique_id for measurement in visible_measurements]
        end;
        item_label="measurements", filter_id="##measurements_filter"
    )

    if !isempty(selected_devices)
        if length(selected_devices) == 1
            sel_name = join(selected_path, "/")
            ig.Text("Measurements for $sel_name")
            ig.Separator()
        else
            shown = min(3, length(selected_devices))
            first_names = join((selected_devices[i].name for i in 1:shown), ", ")
            ig.Text("Measurements from $(length(selected_devices)) devices: $first_names$(length(selected_devices) > 3 ? "..." : "")")
            ig.Separator()
        end
    end

    if isempty(selected_devices)
        ig.Text("Select one or more devices to view their measurements")
        ig.EndChild()
        return nothing
    end

    state.performance.measurement_rows_visible = length(visible_measurements)
    state.performance.measurement_rows_rendered = 0

    if !isempty(visible_measurements)
        table_flags = ig.ImGuiTableFlags_BordersOuterH | ig.ImGuiTableFlags_BordersOuterV |
                      ig.ImGuiTableFlags_RowBg | ig.ImGuiTableFlags_SizingStretchProp
        if ig.BeginTable("measurements_table", 1, table_flags)
            ig.TableSetupColumn("Measurement", ig.ImGuiTableColumnFlags_WidthStretch)

            rows_rendered = 0
            clipper = ig.lib.ImGuiListClipper()
            try
                ig.Begin(clipper, length(visible_measurements), ig.GetTextLineHeightWithSpacing())
                while ig.Step(clipper)
                    start_idx = Int(unsafe_load(clipper.DisplayStart)) + 1
                    end_idx = Int(unsafe_load(clipper.DisplayEnd))
                    for idx in start_idx:end_idx
                        measurement = visible_measurements[idx]
                        ig.PushID(measurement.unique_id)
                        ig.TableNextRow()
                        ig.TableSetColumnIndex(0)
                        rows_rendered += 1

                        is_selected = measurement.unique_id in selected_measurement_id_set
                        bad_text_pushed = false
                        if registry_ready
                            tag_state = _tag_state_or_error(state)
                            dev_key = device_path_key(measurement.device_info)
                            ancestor_keys = _ancestor_keys_for_path(dev_key)
                            bad_text_pushed = _push_tag_text_style!(
                                tag_state, measurement.unique_id, [dev_key; ancestor_keys])
                        end
                        if ig.Selectable(display_label(proj, measurement), is_selected, ig.ImGuiSelectableFlags_SpanAllColumns)
                            io = ig.GetIO()
                            shift_held = unsafe_load(io.KeyShift)
                            ctrl_held = unsafe_load(io.KeyCtrl)
                            selected_measurement_ids = copy(workspace.selection.measurement_ids)
                            visible_measurement_ids = [visible.unique_id for visible in visible_measurements]
                            _update_multi_selection!(
                                selected_measurement_ids,
                                measurement.unique_id,
                                visible_measurement_ids,
                                shift_held,
                                ctrl_held,
                            )
                            workspace.selection.measurement_ids = selected_measurement_ids
                        end
                        if state.scroll_to_measurement_id == measurement.unique_id
                            ig.SetScrollHereY(0.5)
                            state.scroll_to_measurement_id = nothing
                        end
                        bad_text_pushed && ig.PopStyleColor()

                        if ig.BeginPopupContextItem()
                            target_measurements =
                                _selection_targets(selected_measurements, measurement)
                            target_ids = [target.unique_id for target in target_measurements]
                            selected_count = length(target_ids)
                            selected_count > 1 && ig.TextDisabled("Apply to $selected_count measurements")

                            if ig.MenuItem("Open Plot in New Window")
                                plots = state.plots
                                plots.next_window_id += 1
                                push!(
                                    plots.windows,
                                    PlotViewState(
                                        id="plot_$(plots.next_window_id)",
                                        title=measurement.clean_title,
                                        live=false,
                                        measurement_ids=[measurement.unique_id],
                                        plot_kind=get(
                                            plots.kind_by_measurement,
                                            measurement.measurement_kind,
                                            nothing,
                                        ),
                                    ),
                                )
                            end

                            editable = _tag_state_ready(state)
                            if !editable
                                ig.TextDisabled("Fix tags.txt and rescan to edit")
                                ig.Separator()
                            end

                            !editable && ig.BeginDisabled()
                            if ig.MenuItem("Mark Bad")
                                _set_measurements_bad!(state, target_ids, true)
                            end
                            if ig.MenuItem("Unmark Bad")
                                _set_measurements_bad!(state, target_ids, false)
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

            state.performance.measurement_rows_rendered = rows_rendered
            state.performance.measurement_rows_visible = length(visible_measurements)
            ig.EndTable()
        end
    else
        ig.TextDisabled("No measurements match filter")
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
    if state.measurement_filter_widget === nothing
        state.measurement_filter_widget =
            ig.ImGuiTextFilter_ImGuiTextFilter(state.measurement_filter)
    end
    filter_tree = state.tree_filter_widget
    filter_meas = state.measurement_filter_widget

    if ig.Begin("Hierarchy", C_NULL, ig.ImGuiWindowFlags_MenuBar)
        render_menu_bar(state)
        ig.Columns(2, "main_layout")
        _time!(state, :device_tree) do
            _render_hierarchy_tree_panel(state, filter_tree)
        end
        ig.NextColumn()
        _time!(state, :measurement_panel) do
            _render_measurements_panel(state, filter_meas)
        end
    end
    ig.End()
    state.tree_filter = _imgui_text_filter_text(filter_tree)
    state.measurement_filter = _imgui_text_filter_text(filter_meas)
    return nothing
end
