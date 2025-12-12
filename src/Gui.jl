import GLFW
using GLMakie
import GLMakie.Makie as Makie
import CImGui as ig
import CImGui.CSyntax: @c
import ModernGL as gl
using Printf

using Statistics: mean

include("MakieIntegration.jl")
using .MakieImguiIntegration

using DataPlotter: figure_for_file, figure_for_files, get_combined_plot_types

# Timing & allocation utilities
# usage: _time!(ui_state, :key) do ... end
function _time!(f::Function, ui_state, key::Symbol)
    timings = get!(() -> Dict{Symbol,Vector{Float64}}(), ui_state, :_timings)
    allocs = get!(() -> Dict{Symbol,Vector{Int}}(), ui_state, :_allocs)
    t0 = time_ns()
    bytes = @allocated f()
    dt_ms = (time_ns() - t0) / 1e6
    vec = get!(() -> Float64[], timings, key)
    push!(vec, dt_ms)
    length(vec) > 400 && popfirst!(vec)
    avec = get!(() -> Int[], allocs, key)
    push!(avec, bytes)
    length(avec) > 400 && popfirst!(avec)
    nothing
end

function _collect_gl_info!()
    try
        Dict(
            :vendor => unsafe_string(gl.glGetString(gl.GL_VENDOR)),
            :renderer => unsafe_string(gl.glGetString(gl.GL_RENDERER)),
            :version => unsafe_string(gl.glGetString(gl.GL_VERSION)),
            :sl => unsafe_string(gl.glGetString(gl.GL_SHADING_LANGUAGE_VERSION)),
        )
    catch err
        @warn "GL info query failed" error = err
        Dict{Symbol,String}()
    end
end

function _print_perf_summary(ui_state)
    @debug begin
        gi = ui_state[:_gl_info]
        timings = get(ui_state, :_timings, Dict{Symbol,Vector{Float64}}())
        allocs = get(ui_state, :_allocs, Dict{Symbol,Vector{Int}}())
        msg = """\n
        ==== Performance Summary ====
        GL Vendor:   $(get(gi, :vendor, "?"))
        GL Renderer: $(get(gi, :renderer, "?"))
        GL Version:  $(get(gi, :version, "?"))
        """
        if !isempty(timings)
            msg = msg * @sprintf "%-12s %5s %9s %9s %9s %12s %12s\n" "Key" "n" "Mean(ms)" "Max(ms)" "Last(ms)" "AllocMean(KB)" "AllocLast(KB)"
        end
        for k in sort(collect(keys(timings)))
            v = timings[k]
            isempty(v) && continue
            a = get(allocs, k, Int[])
            n = length(v)
            mean_ms = mean(v)
            max_ms = maximum(v)
            last_ms = v[end]
            mean_alloc = isempty(a) ? 0.0 : mean(a) / 1024
            last_alloc = isempty(a) ? 0.0 : a[end] / 1024
            msg = msg * @sprintf "%-12s %5d %9.2f %9.2f %9.2f %12.1f %12.1f\n" String(k) n mean_ms max_ms last_ms mean_alloc last_alloc
        end
        msg * "=============================="
    end
end

function _helpmarker(desc::String)
    ig.TextDisabled("(?)")
    if ig.BeginItemTooltip()
        ig.PushTextWrapPos(ig.GetFontSize() * 35.0)
        ig.TextUnformatted(desc)
        ig.PopTextWrapPos()
        ig.EndTooltip()
    end
end

function render_perf_window(ui_state)
    if !get(ui_state, :show_performance_window, false)
        return
    end

    if ig.Begin("Performance")
        raw_io = ig.GetIO()
        fps = unsafe_load(raw_io.Framerate)
        if fps > 0
            ig.Text(
                "FPS: $(round(fps; digits=1))  Frame: " *
                "$(round(1000 / fps; digits=2)) ms"
            )
        end

        if haskey(ui_state, :_gl_info)
            gi = ui_state[:_gl_info]
            for k in (:vendor, :renderer, :version)
                haskey(gi, k) && ig.Text("GL $(k): $(gi[k])")
            end
        end

        if haskey(ui_state, :_node_count)
            ig.Text("Tree nodes rendered: $(ui_state[:_node_count])")
        end

        timings = get(ui_state, :_timings,
            Dict{Symbol,Vector{Float64}}())
        allocs = get(ui_state, :_allocs,
            Dict{Symbol,Vector{Int}}())

        for (k, v) in timings
            isempty(v) && continue
            a = get(allocs, k, Int[])
            last_ms = round(v[end]; digits=2)
            mean_ms = round(mean(v); digits=2)
            last_alloc = isempty(a) ? 0.0 : round(a[end] / 1024; digits=1)
            mean_alloc = isempty(a) ? 0.0 : round(mean(a) / 1024; digits=1)
            msg = @sprintf "%s: last=%.2f ms  mean=%.2f ms  alloc_last=%.1f KB  alloc_mean=%.1f KB" String(k) last_ms mean_ms last_alloc mean_alloc

            ig.BulletText(msg)
        end

        if ig.Button("Clear timings")
            empty!(get!(() -> Dict{Symbol,Vector{Float64}}(), ui_state, :_timings))
            empty!(get!(() -> Dict{Symbol,Vector{Int}}(), ui_state, :_allocs))
        end
    end
    ig.End()
end

function render_main_window(ui_state)
    if ig.Begin("Measurement Browser", C_NULL,
        ig.ImGuiWindowFlags_MenuBar)
    end
    ig.End()
end

function render_menu_bar(ui_state)
    if ig.BeginMenuBar()
        if ig.BeginMenu("File")
            if ig.MenuItem("Open Folder...")
                path = try
                    String(readchomp(`kdialog --getexistingdirectory`))
                catch
                end
                if !isnothing(path) && !isempty(path)
                    @info "Selected path: $path"
                    hierarchy = scan_directory(path)
                    ui_state[:hierarchy_root] = hierarchy.root
                    ui_state[:all_measurements] = hierarchy.all_measurements
                    ui_state[:root_path] = path
                    ui_state[:has_device_metadata] = hierarchy.has_device_metadata
                    # Collect device-level metadata keys (union)
                    all_params = Set{Symbol}()
                    for m in hierarchy.all_measurements
                        for k in keys(m.device_info.parameters)
                            push!(all_params, k)
                        end
                    end
                    ui_state[:device_metadata_keys] = sort!(collect(all_params); by=String)
                end
            end
            if ig.MenuItem("Reload")
                if haskey(ui_state, :root_path) && !isempty(ui_state[:root_path])
                    @info "Reloading path: $(ui_state[:root_path])"
                    hierarchy = scan_directory(ui_state[:root_path])
                    ui_state[:hierarchy_root] = hierarchy.root
                    ui_state[:all_measurements] = hierarchy.all_measurements
                    ui_state[:has_device_metadata] = hierarchy.has_device_metadata
                    # Collect device-level metadata keys (union)
                    all_params = Set{Symbol}()
                    for m in hierarchy.all_measurements
                        for k in keys(m.device_info.parameters)
                            push!(all_params, k)
                        end
                    end
                    ui_state[:device_metadata_keys] = sort!(collect(all_params); by=String)
                end
            end
            ig.EndMenu()
        end
        if ig.BeginMenu("Debug")
            if ig.MenuItem("Performance Window", C_NULL, get(ui_state, :show_performance_window, false))
                ui_state[:show_performance_window] = !get(ui_state, :show_performance_window, false)
            end
            if ig.MenuItem("Debug Plot Mode", C_NULL, get(ui_state, :debug_plot_mode, false))
                ui_state[:debug_plot_mode] = !get(ui_state, :debug_plot_mode, false)
                # Invalidate cached figures when toggled
                delete!(ui_state, :plot_figure)
                delete!(ui_state, :_last_plotted_path)
                delete!(ui_state, :_last_plotted_mtime)
                # Invalidate additional plot windows
                # (disabled for now since we want to print from the active figure only, usually)
                # if haskey(ui_state, :open_plot_windows)
                #     for entry in ui_state[:open_plot_windows]
                #         if haskey(entry, :figure)
                #             delete!(entry, :figure)
                #         end
                #     end
                # end
            end
            ig.EndMenu()
        end
        ig.EndMenuBar()
    end
end

# Left panel (hierarchy tree) rendering
function _render_hierarchy_tree_panel(ui_state, filter_tree)
    ig.BeginChild("Tree", (0, 0), true)
    ig.SeparatorText("Device Selection")

    # Show selection count
    selected_devices = get!(ui_state, :selected_devices, HierarchyNode[])
    all_devices = HierarchyNode[]
    if haskey(ui_state, :hierarchy_root)
        collect_leaf_nodes!(all_devices, ui_state[:hierarchy_root])
    end
    device_filter_func = (device, filter_obj) -> ig.ImGuiTextFilter_PassFilter(filter_obj, device.name, C_NULL)
    filtered_devices = count_filtered_items(all_devices, filter_tree, device_filter_func)

    render_selection_status(length(selected_devices), filtered_devices, length(all_devices), "devices")

    ig.Text("Filter")
    ig.SameLine()
    _helpmarker("incl,-excl")
    ig.SameLine()
    ig.SetNextItemShortcut(
        ig.ImGuiMod_Ctrl | ig.ImGuiKey_F,
        ig.ImGuiInputFlags_Tooltip
    )
    ig.ImGuiTextFilter_Draw(filter_tree, "##tree_filter", -1)

    # Handle Ctrl+A to select all devices
    get_all_devices_func = (ui_state) -> begin
        devices = HierarchyNode[]
        if haskey(ui_state, :hierarchy_root)
            collect_leaf_nodes!(devices, ui_state[:hierarchy_root])
        end
        devices
    end
    _handle_ctrl_a_selection!(ui_state, :selected_devices, get_all_devices_func, device_filter_func, filter_tree)

    meta_keys = get(ui_state, :device_metadata_keys, Symbol[])

    # Local helpers tied to filter object
    node_matches(node::HierarchyNode) = ig.ImGuiTextFilter_PassFilter(filter_tree, node.name, C_NULL)
    subtree_match(node::HierarchyNode) = node_matches(node) || any(subtree_match(c) for c in children(node))

    function render_node(node::HierarchyNode, path::Vector{String}=String[], force_show::Bool=false)
        # return if neither the node nor any of its descendants match the filter
        force_show || subtree_match(node) || return

        ui_state[:_node_count] += 1
        ig.TableNextRow()
        ig.TableSetColumnIndex(0)

        full_path = vcat(path, node.name)
        direct_match = force_show || node_matches(node)
        unique_id = join(full_path, "/")
        ig.PushID(unique_id)

        # appearance flags
        is_leaf = isempty(children(node))
        selected_devices = get(ui_state, :selected_devices, HierarchyNode[])
        selected = is_leaf && node in selected_devices
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

        opened = ig.TreeNodeEx(is_leaf ? "" : node.name, flags, node.name)
        # Handle device selection with multi-select support
        if ig.IsItemClicked()
            ui_state[:selected_path] = full_path

            if is_leaf
                io = ig.GetIO()
                shift_held = unsafe_load(io.KeyShift)
                ctrl_held = unsafe_load(io.KeyCtrl)
                selected_devices = get!(ui_state, :selected_devices, HierarchyNode[])

                all_devices = get_all_devices(ui_state)
                _handle_multi_select!(selected_devices, node, all_devices, shift_held, ctrl_held)
                ui_state[:selected_devices] = selected_devices

                # Maintain backwards compatibility
                if length(selected_devices) == 1
                    ui_state[:selected_device] = selected_devices[1]
                elseif haskey(ui_state, :selected_device)
                    delete!(ui_state, :selected_device)
                end

                # Note: measurements are selected independently in the measurements panel
            elseif haskey(ui_state, :selected_device)
                delete!(ui_state, :selected_device)
                ui_state[:selected_devices] = HierarchyNode[]
            end
        end

        # Fill metadata columns
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
            else
                # non-leaf left blank
            end
        end

        # render children
        if opened && !is_leaf
            for c in children(node)
                render_node(c, full_path, direct_match)
            end
            ig.TreePop()
        end
        ig.PopID()
    end

    if haskey(ui_state, :hierarchy_root)
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
            for child in children(ui_state[:hierarchy_root])
                render_node(child, String[], false)
            end
            ig.EndTable()
        end
    else
        ig.Text("No data loaded")
    end
    ig.EndChild()
end

# Shared multi-select utility functions

"""
    _handle_multi_select!(selected_items, item, all_items, shift_held, ctrl_held)

Common multi-select logic for both device and measurement panels.
Handles range selection (Shift), toggle selection (Ctrl), and single selection.
Modifies selected_items in place.
"""
function _handle_multi_select!(selected_items::Vector{T}, item::T, all_items::Vector{T}, shift_held::Bool, ctrl_held::Bool) where {T}
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
    _handle_ctrl_a_selection!(ui_state, state_key, all_items_func, filter_func, filter_obj)

Handles Ctrl+A "select all visible" functionality for both device and measurement panels.
Only selects items that pass the current filter. Updates ui_state[state_key] with filtered selection.
"""
function _handle_ctrl_a_selection!(ui_state, state_key::Symbol, all_items_func::Function, filter_func::Function, filter_obj)
    if ig.IsKeyPressed(ig.ImGuiKey_A) && ig.IsWindowFocused()
        io = ig.GetIO()
        if unsafe_load(io.KeyCtrl)
            selected_items = get!(ui_state, state_key, [])
            empty!(selected_items)

            all_items = all_items_func(ui_state)
            for item in all_items
                if filter_func(item, filter_obj)
                    push!(selected_items, item)
                end
            end

            ui_state[state_key] = selected_items
        end
    end
end

"""
    render_selection_status(selected_count, filtered_count, total_count, item_type)

Renders consistent selection status display for both device and measurement panels.
Shows "X/Y items selected" with total count if filtering is active.
Uses appropriate colors: green for multi-select, blue for single, disabled for none.
"""
function render_selection_status(selected_count::Int, filtered_count::Int, total_count::Int, item_type::String)
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

"""
    count_filtered_items(items, filter_obj, filter_func)

Counts how many items pass the current filter. Used for consistent status display
across device and measurement panels. filter_func should return true for items that pass.
"""
function count_filtered_items(items::Vector{T}, filter_obj, filter_func::Function) where {T}
    return count(item -> filter_func(item, filter_obj), items)
end

# Helper function to collect all leaf nodes (devices) from hierarchy
function collect_leaf_nodes!(devices::Vector{HierarchyNode}, node::HierarchyNode)
    if isempty(children(node))
        # This is a leaf node (device)
        push!(devices, node)
    else
        # Recursively collect from children
        for child in children(node)
            collect_leaf_nodes!(devices, child)
        end
    end
end

# Helper function to get device path from hierarchy traversal
function get_device_path(target_device::HierarchyNode, current_node::HierarchyNode, path::Vector{String}=String[])
    if current_node == target_device
        return vcat(path, current_node.name)
    end

    for child in children(current_node)
        result = get_device_path(target_device, child, vcat(path, current_node.name))
        if !isnothing(result)
            return result
        end
    end
    return nothing
end


# Helper function to get all devices for range selection
function get_all_devices(ui_state)
    if !haskey(ui_state, :hierarchy_root)
        return HierarchyNode[]
    end

    all_devices = HierarchyNode[]
    collect_leaf_nodes!(all_devices, ui_state[:hierarchy_root])
    return all_devices
end

# Helper function to update measurements from selected devices - removed as measurements are selected independently

# Right panel (measurements list) rendering
function _render_measurements_panel(ui_state, filter_meas)
    ig.BeginChild("Measurements", (0, 0), true)
    ig.SeparatorText("Measurement Selection")

    # Show selection count and calculate measurement statistics
    selected_measurements = get(ui_state, :selected_measurements, MeasurementInfo[])

    # Get all measurements from currently selected devices
    selected_devices = get(ui_state, :selected_devices, HierarchyNode[])
    all_measurements = MeasurementInfo[]
    if !isempty(selected_devices)
        for device in selected_devices
            append!(all_measurements, device.measurements)
        end
    elseif haskey(ui_state, :selected_device)
        append!(all_measurements, ui_state[:selected_device].measurements)
    end

    measurement_filter_func = (m, filter_obj) -> begin
        if !ig.ImGuiTextFilter_IsActive(filter_obj)
            return true
        end
        # Evaluate the filter against a single combined string so negative tokens work reliably
        filter_text = string(
            display_label(m), "\n",
            m.clean_title, "\n",
            measurement_label(m.measurement_kind)
        )
        return ig.ImGuiTextFilter_PassFilter(filter_obj, filter_text, C_NULL)
    end
    filtered_count = count_filtered_items(all_measurements, filter_meas, measurement_filter_func)

    render_selection_status(length(selected_measurements), filtered_count, length(all_measurements), "measurements")

    ig.Text("Filter")
    ig.SameLine()
    _helpmarker("incl,-excl")
    ig.SameLine()
    ig.SetNextItemShortcut(
        ig.ImGuiMod_Ctrl | ig.ImGuiKey_F,
        ig.ImGuiInputFlags_Tooltip
    )
    ig.ImGuiTextFilter_Draw(filter_meas, "##measurements_filter", -1)

    # Handle Ctrl+A for measurements
    get_all_measurements_func = (ui_state) -> begin
        selected_devices = get(ui_state, :selected_devices, HierarchyNode[])
        candidate_measurements = MeasurementInfo[]
        if !isempty(selected_devices)
            for device in selected_devices
                append!(candidate_measurements, device.measurements)
            end
        elseif haskey(ui_state, :selected_device)
            append!(candidate_measurements, ui_state[:selected_device].measurements)
        end
        candidate_measurements
    end
    _handle_ctrl_a_selection!(ui_state, :selected_measurements, get_all_measurements_func, measurement_filter_func, filter_meas)


    # Determine which measurements to show based on selected devices
    meas_vec = MeasurementInfo[]

    if !isempty(selected_devices)
        if length(selected_devices) == 1
            # Single device selected
            device = selected_devices[1]
            sel_name = join(get(ui_state, :selected_path, [""]), "/")
            ig.Text("Measurements for $sel_name")
            ig.Separator()
            meas_vec = device.measurements
        else
            # Multiple devices selected - show combined measurements
            device_names = [d.name for d in selected_devices]
            ig.Text("Measurements from $(length(selected_devices)) devices: $(join(device_names[1:min(3, end)], ", "))$(length(device_names) > 3 ? "..." : "")")
            ig.Separator()
            for device in selected_devices
                append!(meas_vec, device.measurements)
            end
        end
    elseif haskey(ui_state, :selected_device)
        # Backwards compatibility - single device selection
        device = ui_state[:selected_device]
        sel_name = join(get(ui_state, :selected_path, [""]), "/")
        ig.Text("Measurements for $sel_name")
        ig.Separator()
        meas_vec = device.measurements
    end

    if !isempty(meas_vec)
        any_shown = false
        for m in meas_vec
            passes = measurement_filter_func(m, filter_meas)
            passes || continue
            any_shown = true
            selected_measurements = get!(ui_state, :selected_measurements, MeasurementInfo[])
            is_selected = m in selected_measurements

            if ig.Selectable(display_label(m), is_selected)
                io = ig.GetIO()
                shift_held = unsafe_load(io.KeyShift)
                ctrl_held = unsafe_load(io.KeyCtrl)

                # Use common multi-select logic
                _handle_multi_select!(selected_measurements, m, meas_vec, shift_held, ctrl_held)
                ui_state[:selected_measurements] = selected_measurements

                # Maintain backwards compatibility
                if length(selected_measurements) == 1
                    ui_state[:selected_measurement] = selected_measurements[1]
                elseif haskey(ui_state, :selected_measurement)
                    delete!(ui_state, :selected_measurement)
                end
            end
            # Right-click context menu per measurement entry
            if ig.BeginPopupContextItem()
                if ig.MenuItem("Open Plot in New Window")
                    open_plots = get!(ui_state, :open_plot_windows) do
                        Vector{Dict{Symbol,Any}}()
                    end
                    push!(open_plots, Dict(
                        :filepath => m.filepath,
                        :title => m.clean_title,
                        :params => deepcopy(m.device_info.parameters),
                    ))
                end
                ig.EndPopup()
            end
        end
        !any_shown && ig.TextDisabled("No measurements match filter")
    else
        ig.Text("Select one or more devices to view their measurements")
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
        _render_hierarchy_tree_panel(ui_state, filter_tree)
        ig.NextColumn()
        _render_measurements_panel(ui_state, filter_meas)
    end
    ig.End()
end

# Unified helper: ensure (and cache) a Figure for a filepath with params
function _ensure_plot_figure(ui_state, filepath; kind=nothing, params...)
    # Produce a fresh Figure every time this is called (caller controls call frequency).
    # Avoid caching and reusing the same Figure object across multiple ImGui/Makie
    # windows because sharing a single GLMakie.Figure/Screen texture in multiple
    # ImGui contexts can trigger crashes.
    isfile(filepath) || return nothing
    try
        if get(ui_state, :debug_plot_mode, false)
            @info "Debug mode is ON: plots pass DEBUG flag."
        end
        return figure_for_file(filepath, kind; DEBUG=get(ui_state, :debug_plot_mode, false), params...)
    catch err
        @warn "figure_for_file failed" filepath error = err
        return nothing
    end
end

function render_plot_window(ui_state)
    selected_measurements = get(ui_state, :selected_measurements, MeasurementInfo[])
    combined_type = get(ui_state, :combined_plot_type, nothing)

    # Determine what to plot
    should_update_plot = false
    current_plot_key = nothing

    # Check if user requested combined plot generation
    if get(ui_state, :generate_combined_plot, false) &&
       length(selected_measurements) > 1 && combined_type !== nothing

        ui_state[:generate_combined_plot] = false  # Reset flag

        # Filter measurements for the selected plot type
        compatible_measurements = filter_measurements_for_plot_type(selected_measurements, combined_type)

        if length(compatible_measurements) >= 2
            # Generate combined plot
            fig = nothing
            try
                filepaths = [m.filepath for m in compatible_measurements]
                device_params_list = [m.device_info.parameters for m in compatible_measurements]
                fig = figure_for_files(filepaths, combined_type; device_params_list=device_params_list)
            catch err
                @warn "Combined plot generation failed" error = err
            end

            if fig !== nothing
                ui_state[:plot_figure] = fig
                # Store info for reference
                ui_state[:_last_plot_key] = (combined_type, sort([m.filepath for m in compatible_measurements]))
            else
                delete!(ui_state, :plot_figure)
                delete!(ui_state, :_last_plot_key)
            end
        end
    elseif length(selected_measurements) == 1
        # Single file plot mode (backwards compatibility)
        m = selected_measurements[1]
        filepath = m.filepath
        mtime = Dates.unix2datetime(stat(filepath).mtime)
        current_plot_key = (filepath, mtime)
        last_plot_key = get(ui_state, :_last_plot_key, nothing)

        if current_plot_key != last_plot_key
            # Merge device parameters and measurement parameters
            all_params = merge(m.device_info.parameters, m.parameters)
            fig = _ensure_plot_figure(ui_state, filepath; kind=detect_measurement_kind(m.filename), device_params=all_params)
            if fig !== nothing
                ui_state[:plot_figure] = fig
                ui_state[:_last_plot_key] = current_plot_key
            else
                delete!(ui_state, :plot_figure)
                delete!(ui_state, :_last_plot_key)
            end
        end

        # Maintain backwards compatibility
        ui_state[:selected_measurement] = m
    end

    if ig.Begin("Plot Area")
        if get(ui_state, :debug_plot_mode, false)
            ig.TextColored((0.2, 0.8, 0.2, 1.0), "Debug Plot Mode")
            ig.SameLine()
            _helpmarker("Debug mode is ON: plots have DEBUG flag.")
        end

        if haskey(ui_state, :plot_figure)
            f = ui_state[:plot_figure]
            _time!(ui_state, :makie_fig) do
                MakieFigure("measurement_plot", f; auto_resize_x=true, auto_resize_y=true)
            end
        else
            # Show appropriate message based on state
            if length(selected_measurements) > 1 && combined_type !== nothing
                compatible_files = filter_measurements_for_plot_type(selected_measurements, combined_type)
                if length(compatible_files) < 2
                    ig.TextColored((1.0, 0.6, 0.2, 1.0), "Not enough compatible measurements for $(combined_type)")
                    ig.Text("Selected: $(length(selected_measurements)), Compatible: $(length(compatible_files))")
                    if combined_type === :tlm_analysis
                        ig.TextDisabled("TLM Analysis requires ≥2 TLM 4-point measurements")
                    elseif combined_type === :pund_fatigue
                        ig.TextDisabled("PUND Fatigue requires ≥2 PUND measurements")
                    end
                else
                    ig.TextColored((1.0, 0.4, 0.4, 1.0), "Combined plot generation failed")
                    ig.Text("Check file formats and data quality")
                end
            elseif length(selected_measurements) > 1 && combined_type === nothing
                ig.TextColored((0.6, 0.8, 1.0, 1.0), "Multiple measurements selected ($(length(selected_measurements)))")
                ig.Text("Choose a combined plot type and click Generate in the Combined Plots window")
            elseif length(selected_measurements) == 0
                ig.TextDisabled("Select measurements from the Hierarchy panel")
                ig.TextDisabled("• Single measurement: regular plot")
                ig.TextDisabled("• Multiple measurements + plot type: combined plot")
            else
                ig.TextColored((1.0, 0.4, 0.4, 1.0), "Plot generation failed")
                ig.Text("Check file format and measurement type")
            end
        end

        ig.Separator()
    end
    ig.End()
end



function render_combined_plots_window(ui_state)
    if ig.Begin("Combined Plots")
        ig.Text("Select Combined Plot Type:")

        # Instructions
        ig.TextDisabled("1. Select measurements from the Hierarchy panel")
        ig.TextDisabled("2. Choose a plot type below")
        ig.TextDisabled("3. Click Generate to create combined plot")
        ig.Separator()

        # Plot type selector
        current_type = get(ui_state, :combined_plot_type, nothing)
        type_label = current_type === nothing ? "None" : string(current_type)

        if ig.BeginCombo("Plot Type", type_label)
            # Use extensible plot type system
            for (type_key, type_name, type_description) in get_combined_plot_types()
                if ig.Selectable(type_name, current_type === type_key)
                    ui_state[:combined_plot_type] = type_key
                end
                if type_key !== nothing  # Don't show help for "None"
                    ig.SameLine()
                    _helpmarker(type_description)
                end
            end
            ig.EndCombo()
        end



        # Show info about current selection
        selected_measurements = get(ui_state, :selected_measurements, MeasurementInfo[])
        ig.Separator()

        # Selection status
        if !isempty(selected_measurements)
            ig.TextColored((0.6, 0.8, 1.0, 1.0), "Selected: $(length(selected_measurements)) measurements")
        else
            ig.TextDisabled("No measurements selected")
        end

        # Compatibility check and generate button
        if current_type !== nothing && length(selected_measurements) > 1
            compatible_count = count_compatible_measurements(selected_measurements, current_type)
            if compatible_count >= 2
                ig.TextColored((0.2, 0.8, 0.2, 1.0), "Ready: $compatible_count compatible measurements")
                ig.Separator()
                if ig.Button("Generate Combined Plot", (-1, 0))
                    ui_state[:generate_combined_plot] = true
                end
                if ig.BeginPopupContextItem()
                    if ig.MenuItem("Open in New Window")
                        # Prepare compatible selections
                        local sel_meas = get(ui_state, :selected_measurements, MeasurementInfo[])
                        local comp_meas = filter_measurements_for_plot_type(sel_meas, current_type)
                        if length(comp_meas) >= 2
                            try
                                local paths = [m.filepath for m in comp_meas]
                                local dev_params = [m.device_info.parameters for m in comp_meas]
                                local fig = figure_for_files(paths, current_type; device_params_list=dev_params)
                                if fig !== nothing
                                    local open_plots = get!(ui_state, :open_plot_windows) do
                                        Vector{Dict{Symbol,Any}}()
                                    end
                                    local counter = get!(ui_state, :_combined_plot_counter) do
                                        0
                                    end
                                    ui_state[:_combined_plot_counter] = counter + 1
                                    push!(open_plots, Dict(
                                        :figure => fig,
                                        :title => "Combined: $(string(current_type))",
                                        :id => "combined_$(counter + 1)",
                                    ))
                                end
                            catch err
                                @warn "Failed to generate combined plot in new window" error = err
                            end
                        end
                    end
                    ig.EndPopup()
                end
            else
                ig.TextColored((1.0, 0.6, 0.2, 1.0), "Warning: Only $compatible_count compatible (need 2+)")
                ig.Separator()
                ig.BeginDisabled()
                ig.Button("Generate Combined Plot", (-1, 0))
                ig.EndDisabled()
            end
        elseif current_type !== nothing && length(selected_measurements) <= 1
            ig.TextDisabled("Select multiple measurements for combined plots")
            ig.Separator()
            ig.BeginDisabled()
            ig.Button("Generate Combined Plot", (-1, 0))
            ig.EndDisabled()
        elseif length(selected_measurements) > 1
            ig.TextColored((0.8, 0.8, 0.2, 1.0), "Select a plot type above to continue")
            ig.Separator()
            ig.BeginDisabled()
            ig.Button("Generate Combined Plot", (-1, 0))
            ig.EndDisabled()
        end
    end
    ig.End()
end

function count_compatible_measurements(measurements::Vector{MeasurementInfo}, plot_type::Symbol)
    if plot_type === :tlm_analysis || plot_type === :tlm_temperature
        return count(m -> m.measurement_kind === :tlm4p, measurements)
    elseif plot_type === :pund_fatigue
        return count(m -> (m.measurement_kind === :pund || m.measurement_kind === :wakeup), measurements)
    end
    return 0
end

function filter_measurements_for_plot_type(measurements::Vector{MeasurementInfo}, plot_type::Symbol)
    if plot_type === :tlm_analysis || plot_type === :tlm_temperature
        return filter(m -> m.measurement_kind === :tlm4p, measurements)
    elseif plot_type === :pund_fatigue
        return filter(m -> (m.measurement_kind === :pund || m.measurement_kind === :wakeup), measurements)
    end
    return MeasurementInfo[]
end



function render_info_window(ui_state)
    if ig.Begin("Information Panel")
        flags = ig.ImGuiTableFlags_Borders | ig.ImGuiTableFlags_RowBg | ig.ImGuiTableFlags_ScrollY
        ig.BeginTable("info_cols", 2, flags)
        ig.TableSetupColumn("Device")
        ig.TableSetupColumn("Measurement")
        ig.TableHeadersRow()
        ig.TableNextRow()
        ig.TableNextColumn()

        if haskey(ui_state, :selected_device)
            meas_vec = ui_state[:selected_device].measurements
            sel_name = join(get(ui_state, :selected_path, [""]), "/")
            ig.Text("Location: $sel_name")
            ig.Separator()
            stats = begin
                try
                    get_measurements_stats(meas_vec)
                catch err
                    @warn "Failed to compute stats" error = err
                    Dict{Symbol,Any}()
                end
            end
            if !isempty(stats)
                ig.Text("Stats")
                ig.BulletText("Total: $(stats[:total_measurements])")
                ig.BulletText("Types: $(join(stats[:measurement_types], ", "))")
                if haskey(stats, :first_measurement)
                    ig.BulletText("First: $(stats[:first_measurement]) ")
                    ig.BulletText("Last:  $(stats[:last_measurement]) ")
                end
            else
                ig.TextDisabled("No stats available")
            end
            ig.Separator()
            # Device-level metadata
            if !isempty(meas_vec)
                dev_meta = first(meas_vec).device_info.parameters
                if !isempty(dev_meta)
                    ig.Text("Device metadata")
                    for (k, v) in dev_meta
                        ig.BulletText("$(k): $(v)")
                    end
                else
                    ig.TextDisabled("No metadata parameters found")
                end
            end
        else
            ig.TextDisabled("Select a device to see details")
        end

        ig.TableNextColumn()
        if haskey(ui_state, :selected_measurement)
            m = ui_state[:selected_measurement]
            ig.Text("Title: $(m.clean_title)")
            ig.Separator()
            ig.BulletText("Type: $(measurement_label(m.measurement_kind))")
            ig.BulletText("Timestamp: $(m.timestamp)")
            ig.BulletText("Filename:")
            ig.SameLine()
            ig.TextLinkOpenURL(m.filename, m.filepath)
            ig.Separator()
            if !isempty(m.parameters)
                ig.Text("Parameters")
                for (k, v) in m.parameters
                    ig.BulletText("$(k) = $(v)")
                end
            else
                ig.TextDisabled("No parameters extracted")
            end

        else
            ig.TextDisabled("Select a measurement to view details")
        end
        ig.EndTable()
    end
    ig.End()
end

# ------------------------------------------------------------------
# Modal for missing device metadata (shown each scan when missing)
# ------------------------------------------------------------------
function render_device_info_modal(ui_state)
    # Reset dismissal when root path changes
    current_root = get(ui_state, :root_path, "")
    if get(ui_state, :_modal_last_root_path, "") != current_root
        ui_state[:_modal_last_root_path] = current_root
        ui_state[:dev_info_modal] = true
    end
    # always center
    center = ig.ImVec2(0.5, 0.5)
    @c ig.ImGuiViewport_GetCenter(&center, ig.GetMainViewport())
    ig.SetNextWindowPos(center, ig.ImGuiCond_Always, (0.5, 0.5))

    # Show modal if: missing metadata and user hasn't dismissed it this scan
    if get(ui_state, :dev_info_modal, true) && !get(ui_state, :has_device_metadata, true)
        ig.OpenPopup("Device Metadata Missing")
    end

    opened = get(ui_state, :dev_info_modal, true)

    if @c ig.BeginPopupModal("Device Metadata Missing", &opened, ig.ImGuiWindowFlags_AlwaysAutoResize)
        ig.Text("No device metadata file (device_info.txt) was found.")
        ig.Separator()
        ig.TextWrapped("Create a simple text file named device_info.txt in the TOP folder you opened to add extra info (area, thickness, notes, etc.) for each device.")
        ig.Spacing()
        ig.TextWrapped("How to do it:")
        ig.BulletText("Create a new text file: device_info.txt")
        ig.BulletText("First line (header): device_path, area_um2, t_HZO_nm, ...")
        ig.BulletText("Add a column for each property you want to track.")
        ig.BulletText("Add one line per device. device_path can be just a name (A1) or a full path (CHIP1/SITE1/A2)")
        ig.BulletText("A full path entry overrides a simple name entry for the same leaf.")
        ig.Spacing()
        ig.TextDisabled("Example:")
        ig.TextDisabled("device,   area_um2,   t_HZO_nm,   notes,   active")
        ig.TextDisabled("A1,    12.5,   7.0,   baseline,   true")
        ig.TextDisabled("CHIP1/SITE1/A2,    12.4,   7.0,   override,  true")
        ig.Spacing()
        ig.TextWrapped("Save the file, then rescan or reopen the folder to load these values.")
        ig.Spacing()
        if ig.Button("Got it")
            opened = false
            ig.CloseCurrentPopup()
        end
        ig.EndPopup()
    end
    ui_state[:dev_info_modal] = opened
end

# Render any additional plot windows opened via right-click context menu.
function render_additional_plot_windows(ui_state)
    open_plots = get(ui_state, :open_plot_windows, nothing)
    open_plots === nothing && return
    isempty(open_plots) && return
    to_keep = Vector{Dict{Symbol,Any}}()
    for entry in open_plots
        filepath = get(entry, :filepath, "")
        # Support figure-only entries (no filepath)
        if isempty(filepath) && haskey(entry, :figure)
            title = get(entry, :title, "Combined Plot")
            id = get(entry, :id, "combined_plot")
            open_ref = Ref(true)
            if ig.Begin("Plot: $title###plot_window_$id", open_ref)
                f = entry[:figure]
                _time!(ui_state, :makie_fig) do
                    MakieFigure("measurement_plot_$id", f; auto_resize_x=true, auto_resize_y=true)
                end
                ig.Separator()
                ig.TextDisabled(title)
            end
            ig.End()
            open_ref[] && push!(to_keep, entry)
            continue
        end
        isempty(filepath) && continue
        if !isfile(filepath)
            continue
        end
        title = get(entry, :title, basename(filepath))
        # Refresh / create figure (per-window; no global shared Figure)
        mtime = Dates.unix2datetime(stat(filepath).mtime)
        existing_mtime = get(entry, :mtime, nothing)
        refresh = !haskey(entry, :figure) || existing_mtime != mtime
        if refresh
            k = detect_measurement_kind(basename(filepath))
            fig = haskey(entry, :params) ?
                  _ensure_plot_figure(ui_state, filepath; kind=k, entry[:params]...) :
                  _ensure_plot_figure(ui_state, filepath; kind=k)
            fig === nothing && continue
            entry[:figure] = fig
            entry[:mtime] = mtime
        end
        # Window (allow user to close)
        open_ref = Ref(true)
        if ig.Begin("Plot: $title###plot_window_$filepath", open_ref)
            if haskey(entry, :figure)
                f = entry[:figure]
                # Sanitize id for ImGui (avoid slashes)
                id_str = replace(filepath, '/' => '_')
                _time!(ui_state, :makie_fig) do
                    MakieFigure("measurement_plot_$id_str", f; auto_resize_x=true, auto_resize_y=true)
                end
            else
                ig.Text("No plot available")
            end
            ig.Separator()
            ig.TextDisabled(basename(filepath))
        end
        ig.End()
        open_ref[] && push!(to_keep, entry)
    end
    ui_state[:open_plot_windows] = to_keep
end

function create_window_and_run_loop(root_path::Union{Nothing,String}=nothing; engine=nothing, spawn=1)
    ig.set_backend(:GlfwOpenGL3)
    ui_state = Dict{Symbol,Any}()
    ui_state[:_frame] = 0
    ctx = ig.CreateContext()
    io = ig.GetIO()
    io.ConfigFlags = unsafe_load(io.ConfigFlags) | ig.ImGuiConfigFlags_DockingEnable
    io.ConfigFlags = unsafe_load(io.ConfigFlags) | ig.ImGuiConfigFlags_ViewportsEnable
    io.ConfigFlags = unsafe_load(io.ConfigFlags) | ig.ImGuiConfigFlags_NavEnableKeyboard
    ig.StyleColorsDark()
    if root_path !== nothing && root_path != ""
        hierarchy = scan_directory(root_path)
        ui_state[:hierarchy_root] = hierarchy.root
        ui_state[:all_measurements] = hierarchy.all_measurements
        ui_state[:root_path] = root_path
        ui_state[:has_device_metadata] = hierarchy.has_device_metadata
        all_params = Set{Symbol}()
        for m in hierarchy.all_measurements
            for k in keys(m.device_info.parameters)
                push!(all_params, k)
            end
        end
        ui_state[:device_metadata_keys] = sort!(collect(all_params); by=String)
    end
    first_frame = Ref(true)
    ig.render(
        ctx;
        engine,
        window_size=(1920, 1080),
        window_title="Measurement Browser",
        opengl_version=v"3.3",
        spawn,
        wait_events=false,
        on_exit=() -> _print_perf_summary(ui_state),
    ) do
        ui_state[:_frame] += 1
        if first_frame[] && !haskey(ui_state, :_gl_info)
            ui_state[:_gl_info] = _collect_gl_info!()
            first_frame[] = false
            # try
            #     GLFW.SwapInterval(0)  # disable vsync
            # catch err
            #     @warn "Failed to disable vsync" error = err
            # end
        end
        ig.DockSpaceOverViewport(0, ig.GetMainViewport())
        _time!(ui_state, :device_tree) do
            render_selection_window(ui_state)
        end
        _time!(ui_state, :info) do
            render_info_window(ui_state)
        end
        _time!(ui_state, :plot) do
            render_plot_window(ui_state)
        end
        _time!(ui_state, :extra_plots) do
            render_additional_plot_windows(ui_state)
        end
        _time!(ui_state, :perf_window) do
            render_perf_window(ui_state)
        end
        _time!(ui_state, :combined_plots) do
            render_combined_plots_window(ui_state)
        end
        # Show metadata guidance modal if needed
        render_device_info_modal(ui_state)
    end
end

function start_browser(root_path::Union{Nothing,String}=nothing)
    create_window_and_run_loop(root_path)
end
