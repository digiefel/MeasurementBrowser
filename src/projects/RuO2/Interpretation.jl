using DataLoader: read_wakeup_summary

function _ruo2_cvsweep_schema(path::AbstractString)
    header = open(path, "r") do io
        eof(io) && return nothing
        return chomp(readline(io))
    end
    header === nothing && return nothing

    columns = split(header, ',')
    required = ("Frequency_Hz", "Bias_V", "Cp (F)", "Time_sec")
    all(col -> col in columns, required) || return nothing

    if "G (S)" in columns
        return :conductance
    elseif "Rp (Ohm)" in columns
        return :parallel_resistance
    end
    return nothing
end

function device_path_label(::RuO2Project, device_info::DeviceInfo)
    length(device_info.location) <= 1 && return join(device_info.location, "_")
    return join(device_info.location[2:end], "_")
end

function _ruo2_prefixed_filename_device(filename::AbstractString, prefix_units::Vector{String})
    stem = replace(String(filename), r"\.csv$"i => "")
    prefix = join(prefix_units, "_") * "_"
    startswith(stem, prefix) || return nothing
    rest = stem[length(prefix) + 1:end]
    isempty(rest) && return nothing
    parts = split(rest, '_')
    isempty(parts) && return nothing
    return String(parts[1])
end

function _ruo2_nested_location(indexed::IndexedCsvFile)
    parent = dirname(indexed.filepath)
    tail = String[]
    while true
        name = basename(parent)
        if occursin(REGEX_RUO2_CHIP_DIR, name) && length(tail) == 3
            site, subsite, device_dir = reverse(tail)
            location = [name, site, subsite, device_dir]
            filename_device = _ruo2_prefixed_filename_device(indexed.filename, location[1:3])
            if filename_device === nothing
                return location
            end
            if startswith(device_dir, "RuO2")
                location[end] = filename_device
                return location
            end
            device_dir == filename_device || error(
                "RuO2 directory device '$device_dir' does not match filename device " *
                "'$filename_device' for $(indexed.filepath)",
            )
            return location
        end
        push!(tail, name)
        length(tail) > 3 && return nothing
        next_parent = dirname(parent)
        next_parent == parent && return nothing
        parent = next_parent
    end
end

function _ruo2_location_from_identifier(identifier::AbstractString)
    parts = split(String(identifier), '_')
    length(parts) >= 4 || error("Invalid RuO2 device identifier: $identifier")
    chip = String(join(parts[1:end-3], "_"))
    site = String(parts[end-2])
    subsite = String(parts[end-1])
    device = String(parts[end])
    isempty(chip) && error("Invalid RuO2 device identifier: $identifier")
    return [chip, site, subsite, device]
end

function _ruo2_header_identifier(indexed::IndexedCsvFile)
    raw = get(indexed.header_summary, :device_id, nothing)
    raw isa AbstractString || return nothing
    stripped = strip(String(raw))
    isempty(stripped) && return nothing
    match_obj = match(REGEX_RUO2_IDENTIFIER, stripped)
    match_obj === nothing && error("Invalid RuO2 Device ID header '$stripped' in $(indexed.filepath)")
    return String(match_obj.captures[1])
end

function _ruo2_filename_identifier(filename::AbstractString)
    match_obj = match(REGEX_RUO2_BRACKET_IDENTIFIER, String(filename))
    match_obj === nothing && return nothing
    return String(match_obj.captures[1])
end

function _ruo2_resolve_location(indexed::IndexedCsvFile)
    nested = _ruo2_nested_location(indexed)
    header_identifier = _ruo2_header_identifier(indexed)
    if nested !== nothing
        if header_identifier !== nothing
            header_location = _ruo2_location_from_identifier(header_identifier)
            header_location == nested || error(
                "RuO2 Device ID header '$header_identifier' does not match directory hierarchy " *
                "'$(join(nested, '/'))' for $(indexed.filepath)",
            )
        end
        return nested
    end

    header_identifier !== nothing && return _ruo2_location_from_identifier(header_identifier)

    filename_identifier = _ruo2_filename_identifier(indexed.filename)
    filename_identifier !== nothing && return _ruo2_location_from_identifier(filename_identifier)

    return nothing
end

function parse_device_info(::RuO2Project, indexed::IndexedCsvFile)
    location = _ruo2_resolve_location(indexed)
    location === nothing && error("Could not resolve RuO2 device path from $(indexed.filepath)")
    return DeviceInfo(location)
end

function detect_kind(::RuO2Project, filename::String)::Symbol
    lower = lowercase(filename)
    if occursin("pund_fatigue", lower) || occursin("pund fatigue", lower)
        return :pund_fatigue
    elseif occursin("fe pund", lower) || occursin("fepund", lower) || occursin("_pund", lower)
        return :pund
    elseif occursin("cvsweep", lower)
        return :cvsweep
    elseif occursin("tlm_4p", lower) || occursin("tlm", lower)
        return :tlm4p
    elseif occursin("i_v sweep", lower) || occursin("iv sweep", lower) || occursin("fourterminaliv", lower)
        return :iv
    elseif occursin("break", lower) || occursin("breakdown", lower)
        return :breakdown
    elseif occursin("wakeup", lower)
        return :wakeup
    else
        return :unknown
    end
end

function interpret_file(::RuO2Project, indexed::IndexedCsvFile; should_cancel::Union{Nothing,Function}=nothing)::Vector{MeasurementItem}
    kind = detect_kind(RUO2_PROJECT, indexed.filename)
    kind == :unknown && return MeasurementItem[]
    if kind == :cvsweep && _ruo2_cvsweep_schema(indexed.filepath) === nothing
        @warn "Ignoring unsupported RuO2 CVSweep file" path = indexed.filepath
        return MeasurementItem[]
    end
    device_info = parse_device_info(RUO2_PROJECT, indexed)

    params = parse_parameters(indexed.filename)
    if kind == :wakeup
        summary = read_wakeup_summary(indexed.filename, dirname(indexed.filepath))
        params[:wakeup_pulse_count] = summary.pulse_count
        params[:amplitude_V] = summary.amplitude
    end

    title = build_clean_title(RUO2_PROJECT, indexed.filename, kind, device_info, indexed.header_summary)
    base = MeasurementItem(
        item_id(indexed.id),
        indexed.id,
        indexed.filepath,
        kind,
        copy(device_info.location),
        indexed.timestamp,
        Dict{Symbol,Any}(),
        params,
        title,
    )

    expanded = _ruo2_expand_multi_device_item(base)
    if kind == :pund_fatigue
        return vcat([_ruo2_expand_pund_fatigue_item(item; should_cancel=should_cancel) for item in expanded]...)
    end
    return expanded
end

function _ruo2_expand_multi_device_item(item::MeasurementItem)::Vector{MeasurementItem}
    item.kind == :breakdown || return [item]
    dev = last(item.device_path)
    if (m = match(r"^([A-Z][0-9]+)([A-Z][0-9]+)$", dev)) === nothing
        return [item]
    end
    parts = m.captures
    loc = copy(item.device_path)
    return [MeasurementItem(
        item_id(item.source_file_id; split=p),
        item.source_file_id,
        item.filepath,
        item.kind,
        vcat(loc[1:end-1], [p]),
        item.timestamp,
        deepcopy(item.device_parameters),
        deepcopy(item.parameters),
        replace(item.title, dev => p),
    ) for p in parts]
end

function _ruo2_expand_pund_fatigue_item(item::MeasurementItem; should_cancel::Union{Nothing,Function}=nothing)::Vector{MeasurementItem}
    item.kind == :pund_fatigue || return [item]
    cycles, voltage_V = try
        _ruo2_scan_fatigue_file(item.filepath; should_cancel=should_cancel)
    catch e
        e isa ScanCancelled && rethrow()
        @warn "Could not read fatigue cycles from $(item.filepath)" error = e
        return MeasurementItem[]
    end
    isempty(cycles) && return MeasurementItem[]
    return [MeasurementItem(
        item_id(item.source_file_id; cycle=c),
        item.source_file_id,
        item.filepath,
        :pund,
        copy(item.device_path),
        item.timestamp,
        deepcopy(item.device_parameters),
        merge(deepcopy(item.parameters), Dict{Symbol,Any}(:fatigue_cycle => c, :voltage_V => voltage_V)),
        item.title * " cycle $c",
    ) for c in cycles]
end
