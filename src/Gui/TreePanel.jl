# Left panel (hierarchy tree) rendering
function _tree_node_matches_filter(filter_tree, node::HierarchyNode)
    return ig.ImGuiTextFilter_PassFilter(filter_tree, node.name, C_NULL)
end

function _measurement_filter_text(proj, measurement::MeasurementInfo)
    return string(
        display_label(proj, measurement), "\n",
        measurement.clean_title, "\n",
        kind_label(proj, measurement.measurement_kind),
    )
end

function _measurement_matches_filter(proj, measurement::MeasurementInfo, filter_obj)
    if !ig.ImGuiTextFilter_IsActive(filter_obj)
        return true
    end
    return ig.ImGuiTextFilter_PassFilter(filter_obj, _measurement_filter_text(proj, measurement), C_NULL)
end

function _visible_measurements(ui_state, proj, measurements, filter_meas)
    if get(ui_state, :show_bad, true) && !ig.ImGuiTextFilter_IsActive(filter_meas)
        return measurements
    end

    visible = MeasurementInfo[]
    sizehint!(visible, length(measurements))
    for measurement in measurements
        _measurement_is_visible(ui_state, measurement) || continue
        _measurement_matches_filter(proj, measurement, filter_meas) || continue
        push!(visible, measurement)
    end
    return visible
end

function _render_hierarchy_tree_panel(ui_state, filter_tree)
    ig.BeginChild("Tree", (0, 0), true)
    ig.SeparatorText("Device Selection")
    _render_bad_registry_error!(ui_state)

    root = get(ui_state, :hierarchy_root, nothing)
    meta_keys = get(ui_state, :device_metadata_keys, Symbol[])

    visible_devices = HierarchyNode[]
    visible_device_keys_ref = Ref{Union{Nothing,Vector{String}}}(nothing)
    selected_device_path_set = Set(get(ui_state, :selected_device_paths, String[]))
    all_device_count = 0

    if root !== nothing
        has_visible_leaf_cache = IdDict{HierarchyNode,Bool}()
        subtree_matches_cache = IdDict{HierarchyNode,Bool}()
        device_key_cache = IdDict{HierarchyNode,String}()

        device_key(node::HierarchyNode) = get!(device_key_cache, node) do
            _device_path_key(node)
        end

        function has_visible_leaf(node::HierarchyNode)
            return get!(has_visible_leaf_cache, node) do
                if isempty(children(node))
                    return _device_is_visible(ui_state, device_key(node))
                end
                for child in children(node)
                    has_visible_leaf(child) && return true
                end
                return false
            end
        end

        function subtree_matches(node::HierarchyNode)
            return get!(subtree_matches_cache, node) do
                _tree_node_matches_filter(filter_tree, node) && return true
                for child in children(node)
                    subtree_matches(child) && return true
                end
                return false
            end
        end

        function collect_visible_devices!(node::HierarchyNode, force_show::Bool=false)
            has_visible_leaf(node) || return
            force_show || subtree_matches(node) || return

            direct_match = force_show || _tree_node_matches_filter(filter_tree, node)
            if isempty(children(node))
                push!(visible_devices, node)
                return
            end
            for child in children(node)
                collect_visible_devices!(child, direct_match)
            end
        end

        function count_leaf_nodes(node::HierarchyNode)
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
            length(get(ui_state, :selected_devices, HierarchyNode[])),
            length(visible_devices),
            all_device_count,
            filter_tree,
            () -> begin
                ui_state[:selected_device_paths] = copy(visible_device_keys())
                _apply_visible_selection!(ui_state)
            end;
            item_label="devices", filter_id="##tree_filter"
        )

        node_seed = UInt64(0x9e3779b97f4a7c15)
        next_node_id(parent_id::UInt64, node::HierarchyNode) = hash(node.name, parent_id)
        to_imgui_id(node_id::UInt64) = Int32(node_id % UInt64(typemax(Int32)))

        function render_node(node::HierarchyNode, node_id::UInt64, force_show::Bool=false)
            has_visible_leaf(node) || return
            force_show || subtree_matches(node) || return

            ui_state[:_node_count] += 1
            ig.TableNextRow()

            direct_match = force_show || _tree_node_matches_filter(filter_tree, node)
            ig.PushID(to_imgui_id(node_id))

            is_leaf = isempty(children(node))
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
            if ig.ImGuiTextFilter_IsActive(filter_tree) && direct_match && !is_leaf
                flags |= ig.ImGuiTreeNodeFlags_DefaultOpen
            end

            ig.TableSetColumnIndex(0)
            bad_text_pushed = is_leaf && leaf_device_key !== nothing &&
                              _bad_registry_ready(ui_state) &&
                              _push_bad_text_style!(_device_is_explicitly_bad(ui_state, leaf_device_key))
            opened = ig.TreeNodeEx(is_leaf ? "" : node.name, flags, node.name)
            bad_text_pushed && ig.PopStyleColor()

            if ig.IsItemClicked()
                if is_leaf
                    io = ig.GetIO()
                    shift_held = unsafe_load(io.KeyShift)
                    ctrl_held = unsafe_load(io.KeyCtrl)
                    selected_device_paths = copy(get(ui_state, :selected_device_paths, String[]))
                    _update_multi_selection!(selected_device_paths, leaf_device_key, visible_device_keys(), shift_held, ctrl_held)
                    ui_state[:selected_device_paths] = selected_device_paths
                    _apply_visible_selection!(ui_state)
                elseif !isempty(get(ui_state, :selected_device_paths, String[]))
                    ui_state[:selected_device_paths] = String[]
                    _apply_visible_selection!(ui_state)
                end
            end

            if is_leaf && ig.BeginPopupContextItem()
                target_nodes = _selection_targets(get(ui_state, :selected_devices, HierarchyNode[]), node)
                target_keys = [device_key(target) for target in target_nodes]
                selected_count = length(target_keys)
                selected_count > 1 && ig.TextDisabled("Apply to $selected_count devices")

                editable = _bad_registry_ready(ui_state)
                if !editable
                    ig.TextDisabled("Fix bad_measurements and rescan to edit")
                    ig.Separator()
                end

                !editable && ig.BeginDisabled()
                if ig.MenuItem("Mark Bad")
                    _set_devices_bad!(ui_state, target_keys, true)
                end
                if ig.MenuItem("Unmark Bad")
                    _set_devices_bad!(ui_state, target_keys, false)
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
                    render_node(child, next_node_id(node_id, child), direct_match)
                end
                ig.TreePop()
            end
            ig.PopID()
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
                render_node(child, next_node_id(node_seed, child), false)
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

"""
    _update_multi_selection!(selected_items, item, all_items, shift_held, ctrl_held)

Common multi-select logic for both device and measurement panels.
Handles range selection (Shift), toggle selection (Ctrl), and single selection.
Modifies selected_items in place.
"""
function _update_multi_selection!(selected_items::Vector{T}, item::T, all_items::Vector{T}, shift_held::Bool, ctrl_held::Bool) where {T}
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
end

"""
    _render_selection_status!(selected_count, filtered_count, total_count, item_type)

Renders consistent selection status display for both device and measurement panels.
Shows "X/Y items selected" with total count if filtering is active.
Uses appropriate colors: green for multi-select, blue for single, disabled for none.
"""
function _render_selection_status!(selected_count::Int, filtered_count::Int, total_count::Int, item_type::String)
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
end

function _render_selection_toolbar!(
    selected_count::Int,
    visible_count::Int,
    total_count::Int,
    filter_obj,
    on_select_all!::Function;
    item_label::String,
    filter_id::String,
)
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
end

function _collect_leaf_nodes!(devices::Vector{HierarchyNode}, node::HierarchyNode)
    if isempty(children(node))
        push!(devices, node)
        return
    end
    for child in children(node)
        _collect_leaf_nodes!(devices, child)
    end
end

function _all_devices(ui_state)
    root = get(ui_state, :hierarchy_root, nothing)
    root === nothing && return HierarchyNode[]

    devices = HierarchyNode[]
    _collect_leaf_nodes!(devices, root)
    return devices
end

function _selected_measurements(ui_state)
    haskey(ui_state, :selected_all_measurements) ||
        error("UI selection state is missing selected_all_measurements")
    return ui_state[:selected_all_measurements]
end

# Right panel (measurements list) rendering
function _render_measurements_panel(ui_state, filter_meas)
    proj = ui_state[:project]
    ig.BeginChild("Measurements", (0, 0), true)
    ig.SeparatorText("Measurement Selection")

    selected_devices = get(ui_state, :selected_devices, HierarchyNode[])
    all_measurements = _selected_measurements(ui_state)
    visible_measurements = _visible_measurements(ui_state, proj, all_measurements, filter_meas)
    haskey(ui_state, :selected_measurement_id_set) ||
        error("UI selection state is missing selected_measurement_id_set")
    selected_measurement_id_set = ui_state[:selected_measurement_id_set]
    registry_ready = _bad_registry_ready(ui_state)

    _render_selection_toolbar!(
        length(get(ui_state, :selected_measurements, MeasurementInfo[])),
        length(visible_measurements),
        length(all_measurements),
        filter_meas,
        () -> begin
            ui_state[:selected_measurement_ids] = [measurement.id for measurement in visible_measurements]
            _apply_visible_selection!(ui_state)
        end;
        item_label="measurements", filter_id="##measurements_filter"
    )

    if !isempty(selected_devices)
        if length(selected_devices) == 1
            sel_name = join(get(ui_state, :selected_path, [""]), "/")
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
        return
    end

    ui_state[:_measurement_rows_visible] = length(visible_measurements)
    ui_state[:_measurement_rows_rendered] = 0

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
                        ig.PushID(measurement.id)
                        ig.TableNextRow()
                        ig.TableSetColumnIndex(0)
                        rows_rendered += 1

                        is_selected = measurement.id in selected_measurement_id_set
                        bad_text_pushed = registry_ready &&
                                          _push_bad_text_style!(_measurement_is_bad(ui_state, measurement))
                        if ig.Selectable(display_label(proj, measurement), is_selected, ig.ImGuiSelectableFlags_SpanAllColumns)
                            io = ig.GetIO()
                            shift_held = unsafe_load(io.KeyShift)
                            ctrl_held = unsafe_load(io.KeyCtrl)
                            selected_measurement_ids = copy(get(ui_state, :selected_measurement_ids, String[]))
                            visible_measurement_ids = [visible.id for visible in visible_measurements]
                            _update_multi_selection!(
                                selected_measurement_ids,
                                measurement.id,
                                visible_measurement_ids,
                                shift_held,
                                ctrl_held,
                            )
                            ui_state[:selected_measurement_ids] = selected_measurement_ids
                            _apply_visible_selection!(ui_state)
                        end
                        bad_text_pushed && ig.PopStyleColor()

                        if ig.BeginPopupContextItem()
                            target_measurements = _selection_targets(get(ui_state, :selected_measurements, MeasurementInfo[]), measurement)
                            target_ids = [target.id for target in target_measurements]
                            selected_count = length(target_ids)
                            selected_count > 1 && ig.TextDisabled("Apply to $selected_count measurements")

                            if ig.MenuItem("Open Plot in New Window")
                                open_plots = get!(ui_state, :open_plot_windows) do
                                    Vector{Dict{Symbol,Any}}()
                                end
                                push!(open_plots, _measurement_plot_window_entry(measurement))
                            end

                            editable = _bad_registry_ready(ui_state)
                            if !editable
                                ig.TextDisabled("Fix bad_measurements and rescan to edit")
                                ig.Separator()
                            end

                            !editable && ig.BeginDisabled()
                            if ig.MenuItem("Mark Bad")
                                _set_measurements_bad!(ui_state, target_ids, true)
                            end
                            if ig.MenuItem("Unmark Bad")
                                _set_measurements_bad!(ui_state, target_ids, false)
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

            ui_state[:_measurement_rows_rendered] = rows_rendered
            ui_state[:_measurement_rows_visible] = length(visible_measurements)
            ig.EndTable()
        end
    else
        ig.TextDisabled("No measurements match filter")
    end
    ig.EndChild()
end

function render_selection_window(ui_state)
    ui_state[:_node_count] = 0
    filter_tree = get!(ui_state, :_imgui_text_filter_tree) do
        ig.ImGuiTextFilter_ImGuiTextFilter(C_NULL)
    end
    filter_meas = get!(ui_state, :_imgui_text_filter_meas) do
        ig.ImGuiTextFilter_ImGuiTextFilter(C_NULL)
    end

    if ig.Begin("Hierarchy", C_NULL, ig.ImGuiWindowFlags_MenuBar)
        render_menu_bar(ui_state)
        ig.Columns(2, "main_layout")
        _time!(ui_state, :device_tree) do
            _render_hierarchy_tree_panel(ui_state, filter_tree)
        end
        ig.NextColumn()
        _time!(ui_state, :measurement_panel) do
            _render_measurements_panel(ui_state, filter_meas)
        end
    end
    ig.End()
end
