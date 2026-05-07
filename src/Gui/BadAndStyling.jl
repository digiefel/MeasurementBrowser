function _copy_bad_registry(registry::BadRegistry)
    return BadRegistry(copy(registry.devices), copy(registry.measurements))
end

function _load_bad_registry_for_root!(ui_state, root_path::String)
    if isempty(root_path)
        ui_state[:bad_registry] = nothing
        ui_state[:bad_registry_error] = ""
        return
    end

    try
        ui_state[:bad_registry] = load_bad_registry(root_path)
        ui_state[:bad_registry_error] = ""
        _apply_visible_selection!(ui_state)
    catch err
        if err isa BadRegistryParseError || err isa BadRegistryIOError
            ui_state[:bad_registry] = nothing
            ui_state[:bad_registry_error] = sprint(showerror, err)
            ui_state[:show_bad] = true
            _apply_visible_selection!(ui_state)
            return
        end
        rethrow()
    end
end

function _bad_registry_ready(ui_state)
    return get(ui_state, :bad_registry, nothing) isa BadRegistry && isempty(get(ui_state, :bad_registry_error, ""))
end

function _bad_registry_or_error(ui_state)
    registry = get(ui_state, :bad_registry, nothing)
    registry isa BadRegistry || error("Bad registry unavailable: $(get(ui_state, :bad_registry_error, ""))")
    return registry
end

function _device_path_key(node::HierarchyNode)
    isempty(node.measurements) && error("Leaf node '$(node.name)' has no measurements")
    return device_path_key(first(node.measurements).device_info)
end

function _device_location(node::HierarchyNode)
    isempty(node.measurements) && error("Leaf node '$(node.name)' has no measurements")
    return copy(first(node.measurements).device_info.location)
end

function _device_is_explicitly_bad(ui_state, device_key::String)
    return device_key in _bad_registry_or_error(ui_state).devices
end

function _measurement_is_explicitly_bad(ui_state, measurement_id::String)
    return measurement_id in _bad_registry_or_error(ui_state).measurements
end

function _measurement_is_bad(ui_state, measurement::MeasurementInfo)
    return _measurement_is_explicitly_bad(ui_state, measurement.id) ||
           _device_is_explicitly_bad(ui_state, device_path_key(measurement.device_info))
end

function _assert_bad_registry_visibility_available(ui_state)
    _bad_registry_ready(ui_state) && return
    error("Cannot hide bad items while bad registry is unavailable: $(get(ui_state, :bad_registry_error, ""))")
end

function _device_is_visible(ui_state, device_key::String)
    get(ui_state, :show_bad, true) && return true
    _assert_bad_registry_visibility_available(ui_state)
    return !_device_is_explicitly_bad(ui_state, device_key)
end

function _measurement_is_visible(ui_state, measurement::MeasurementInfo)
    get(ui_state, :show_bad, true) && return true
    _assert_bad_registry_visibility_available(ui_state)
    return !_measurement_is_bad(ui_state, measurement)
end

function _project_visible_selection(ui_state)
    hierarchy = get(ui_state, :scan_hierarchy, nothing)
    if hierarchy === nothing
        return HierarchyNode[], MeasurementInfo[], String[]
    end

    selected_devices = HierarchyNode[]
    for path_key in get(ui_state, :selected_device_paths, String[])
        node = get(hierarchy.index, device_path_tuple(path_key), nothing)
        node === nothing && continue
        node.kind == :leaf || error("Selected device path '$path_key' does not point to a leaf device")
        !_device_is_visible(ui_state, path_key) && continue
        push!(selected_devices, node)
    end

    visible_device_keys = Set(_device_path_key(node) for node in selected_devices)
    measurement_index = get(ui_state, :measurement_index, Dict{String,MeasurementInfo}())
    selected_measurements = MeasurementInfo[]
    for measurement_id in get(ui_state, :selected_measurement_ids, String[])
        measurement = get(measurement_index, measurement_id, nothing)
        measurement === nothing && continue
        device_path_key(measurement.device_info) in visible_device_keys || continue
        !_measurement_is_visible(ui_state, measurement) && continue
        push!(selected_measurements, measurement)
    end
    selected_path = length(selected_devices) == 1 ? _device_location(selected_devices[1]) : String[]
    return selected_devices, selected_measurements, selected_path
end

function _apply_visible_selection!(ui_state)
    selected_devices, selected_measurements, selected_path = _project_visible_selection(ui_state)

    all_measurements = MeasurementInfo[]
    sizehint!(all_measurements, sum(length(device.measurements) for device in selected_devices; init=0))
    for device in selected_devices
        for measurement in device.measurements
            push!(all_measurements, measurement)
        end
    end

    ui_state[:selected_devices] = selected_devices
    ui_state[:selected_measurements] = selected_measurements
    ui_state[:selected_all_measurements] = all_measurements
    ui_state[:selected_measurement_id_set] = Set(get(ui_state, :selected_measurement_ids, String[]))
    ui_state[:selected_path] = selected_path
end

function _set_devices_bad!(ui_state, device_keys::Vector{String}, bad::Bool)
    unique_keys = unique(copy(device_keys))
    isempty(unique_keys) && return false
    _bad_registry_ready(ui_state) || return false

    root_path = get(ui_state, :root_path, "")
    isempty(root_path) && error("Cannot update bad registry without an active project root")

    registry = get(ui_state, :bad_registry, nothing)
    registry isa BadRegistry || error("Bad registry is unavailable for editing")

    updated = _copy_bad_registry(registry)
    for device_key in unique_keys
        bad ? push!(updated.devices, device_key) : delete!(updated.devices, device_key)
    end

    try
        save_bad_registry(root_path, updated)
    catch err
        if err isa BadRegistryIOError
            ui_state[:bad_registry_error] = sprint(showerror, err)
            return false
        end
        rethrow()
    end

    ui_state[:bad_registry] = updated
    ui_state[:bad_registry_error] = ""
    _apply_visible_selection!(ui_state)
    return true
end

function _set_measurements_bad!(ui_state, measurement_ids::Vector{String}, bad::Bool)
    unique_ids = unique(copy(measurement_ids))
    isempty(unique_ids) && return false
    _bad_registry_ready(ui_state) || return false

    root_path = get(ui_state, :root_path, "")
    isempty(root_path) && error("Cannot update bad registry without an active project root")

    registry = get(ui_state, :bad_registry, nothing)
    registry isa BadRegistry || error("Bad registry is unavailable for editing")

    updated = _copy_bad_registry(registry)
    for measurement_id in unique_ids
        bad ? push!(updated.measurements, measurement_id) : delete!(updated.measurements, measurement_id)
    end

    try
        save_bad_registry(root_path, updated)
    catch err
        if err isa BadRegistryIOError
            ui_state[:bad_registry_error] = sprint(showerror, err)
            return false
        end
        rethrow()
    end

    ui_state[:bad_registry] = updated
    ui_state[:bad_registry_error] = ""
    _apply_visible_selection!(ui_state)
    return true
end

function _selection_targets(selected_items::Vector{T}, clicked_item::T) where {T}
    if clicked_item in selected_items
        return selected_items
    end
    return T[clicked_item]
end

function _render_bad_registry_error!(ui_state)
    message = get(ui_state, :bad_registry_error, "")
    isempty(message) && return
    ig.TextColored((1.0, 0.5, 0.5, 1.0), "bad_measurements error")
    if ig.BeginItemTooltip()
        ig.PushTextWrapPos(ig.GetFontSize() * 35.0)
        ig.TextUnformatted(message)
        ig.PopTextWrapPos()
        ig.EndTooltip()
    end
end

function _push_bad_text_style!(bad::Bool)
    bad || return false
    ig.PushStyleColor(ig.ImGuiCol_Text, (0.82, 0.35, 0.35, 1.0))
    return true
end
