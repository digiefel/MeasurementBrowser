function _load_tag_state_for_root!(ui_state, root_path::String)
    if isempty(root_path)
        ui_state[:tag_state] = nothing
        ui_state[:tag_state_error] = ""
        return
    end

    try
        ui_state[:tag_state] = Annotations.Tags.load(root_path)
        ui_state[:tag_state_error] = ""
        _apply_visible_selection!(ui_state)
    catch err
        if err isa Annotations.Tags.TagsParseError || err isa IOError
            ui_state[:tag_state] = nothing
            ui_state[:tag_state_error] = sprint(showerror, err)
            _apply_visible_selection!(ui_state)
            return
        end
        rethrow()
    end
end

function _tag_state_ready(ui_state)
    return get(ui_state, :tag_state, nothing) isa Annotations.Tags.TagState &&
           isempty(get(ui_state, :tag_state_error, ""))
end

function _show_bad_effective(ui_state)
    return get(ui_state, :show_bad, true) || !_tag_state_ready(ui_state)
end

function _tag_state_or_error(ui_state)
    state = get(ui_state, :tag_state, nothing)
    state isa Annotations.Tags.TagState ||
        error("Tag state unavailable: $(get(ui_state, :tag_state_error, ""))")
    return state
end

function _device_path_key(node::HierarchyNode)
    isempty(node.measurements) && error("Leaf node '$(node.name)' has no measurements")
    return device_path_key(first(node.measurements).device_info)
end

function _device_location(node::HierarchyNode)
    isempty(node.measurements) && error("Leaf node '$(node.name)' has no measurements")
    return copy(first(node.measurements).device_info.location)
end

function _ancestor_keys_for_path(path::AbstractString)
    parts = split(path, '/')
    length(parts) <= 1 && return String[]
    out = String[]
    for i in 1:(length(parts) - 1)
        push!(out, join(parts[1:i], '/'))
    end
    return out
end

function _has_bad_tag(tag_state::Annotations.Tags.TagState, key::AbstractString,
                      ancestor_keys)
    return "bad" in Annotations.Tags.effective(tag_state, key, ancestor_keys)
end

function _device_is_explicitly_bad(ui_state, device_key::String)
    tag_state = _tag_state_or_error(ui_state)
    return _has_bad_tag(tag_state, device_key, _ancestor_keys_for_path(device_key))
end

function _measurement_is_explicitly_bad(ui_state, measurement_id::String)
    tag_state = _tag_state_or_error(ui_state)
    return _has_bad_tag(tag_state, measurement_id, String[])
end

function _measurement_is_bad(ui_state, measurement::MeasurementInfo)
    tag_state = _tag_state_or_error(ui_state)
    dev_key = device_path_key(measurement.device_info)
    ancestor_keys = _ancestor_keys_for_path(dev_key)
    return _has_bad_tag(tag_state, measurement.id, [dev_key; ancestor_keys])
end

function _assert_tag_state_visibility_available(ui_state)
    _tag_state_ready(ui_state) && return
    error("Cannot hide bad items while tag state is unavailable: $(get(ui_state, :tag_state_error, ""))")
end

function _device_is_visible(ui_state, device_key::String)
    _show_bad_effective(ui_state) && return true
    return !_device_is_explicitly_bad(ui_state, device_key)
end

function _measurement_is_visible(ui_state, measurement::MeasurementInfo)
    _show_bad_effective(ui_state) && return true
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

function _ensure_bad_catalog_entry!(tag_state::Annotations.Tags.TagState)
    any(t -> t.name == Annotations.Tags.LEGACY_BAD_TAG, tag_state.catalog) && return
    pushfirst!(tag_state.catalog,
        Annotations.Tags.TagDef(Annotations.Tags.LEGACY_BAD_TAG,
                                Annotations.Tags.LEGACY_BAD_COLOR,
                                Annotations.Tags.LEGACY_BAD_PRIORITY))
end

function _set_devices_bad!(ui_state, device_keys::Vector{String}, bad::Bool)
    unique_keys = unique(copy(device_keys))
    isempty(unique_keys) && return false
    _tag_state_ready(ui_state) || return false

    root_path = get(ui_state, :root_path, "")
    isempty(root_path) && error("Cannot update tag state without an active project root")

    tag_state = _tag_state_or_error(ui_state)
    _ensure_bad_catalog_entry!(tag_state)

    for device_key in unique_keys
        if bad
            set = get!(() -> Set{String}(), tag_state.assignments, device_key)
            push!(set, Annotations.Tags.LEGACY_BAD_TAG)
        else
            tags = get(tag_state.assignments, device_key, nothing)
            if tags !== nothing
                delete!(tags, Annotations.Tags.LEGACY_BAD_TAG)
                isempty(tags) && delete!(tag_state.assignments, device_key)
            end
        end
    end

    Annotations.Tags.save(root_path, tag_state)
    _apply_visible_selection!(ui_state)
    return true
end

function _set_measurements_bad!(ui_state, measurement_ids::Vector{String}, bad::Bool)
    unique_ids = unique(copy(measurement_ids))
    isempty(unique_ids) && return false
    _tag_state_ready(ui_state) || return false

    root_path = get(ui_state, :root_path, "")
    isempty(root_path) && error("Cannot update tag state without an active project root")

    tag_state = _tag_state_or_error(ui_state)
    _ensure_bad_catalog_entry!(tag_state)

    for measurement_id in unique_ids
        if bad
            set = get!(() -> Set{String}(), tag_state.assignments, measurement_id)
            push!(set, Annotations.Tags.LEGACY_BAD_TAG)
        else
            tags = get(tag_state.assignments, measurement_id, nothing)
            if tags !== nothing
                delete!(tags, Annotations.Tags.LEGACY_BAD_TAG)
                isempty(tags) && delete!(tag_state.assignments, measurement_id)
            end
        end
    end

    Annotations.Tags.save(root_path, tag_state)
    _apply_visible_selection!(ui_state)
    return true
end

function _selection_targets(selected_items::Vector{T}, clicked_item::T) where {T}
    if clicked_item in selected_items
        return selected_items
    end
    return T[clicked_item]
end

function _render_tag_state_error!(ui_state)
    message = get(ui_state, :tag_state_error, "")
    isempty(message) && return
    ig.TextColored((1.0, 0.5, 0.5, 1.0), "tags error")
    if ig.BeginItemTooltip()
        ig.PushTextWrapPos(ig.GetFontSize() * 35.0)
        ig.TextUnformatted(message)
        ig.PopTextWrapPos()
        ig.EndTooltip()
    end
end

function _push_tag_text_style!(tag_state::Annotations.Tags.TagState, key::AbstractString,
                                ancestor_keys)
    eff = Annotations.Tags.effective(tag_state, key, ancestor_keys)
    color = Annotations.Tags.dominant_color(tag_state, eff)
    color === nothing && return false
    r = Float32(color[1]) / 255f0
    g = Float32(color[2]) / 255f0
    b = Float32(color[3]) / 255f0
    ig.PushStyleColor(ig.ImGuiCol_Text, (r, g, b, 1.0f0))
    return true
end
