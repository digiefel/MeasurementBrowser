using DataLoader: read_pund_wakeup_reps, cv_sweep_has_schema

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

function _ruo2_nested_location(indexed::SourceFile)
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

function _ruo2_header_identifier(indexed::SourceFile)
    raw = get(indexed.header_summary, :device_id, nothing)
    raw isa AbstractString || return nothing
    stripped = strip(String(raw))
    isempty(stripped) && return nothing
    match_obj = match(REGEX_RUO2_IDENTIFIER, stripped)
    match_obj === nothing && error("Invalid RuO2 Device ID header '$stripped' in $(indexed.filepath)")
    return String(match_obj.captures[1])
end

function _ruo2_filename_identifier(filename::AbstractString)
    s = String(filename)
    m = match(REGEX_RUO2_BRACKET_IDENTIFIER, s)
    m !== nothing && return String(m.captures[1])
    m = match(REGEX_RUO2_TIMESTAMP_FILENAME, s)
    m !== nothing && return String(m.captures[1])
    return nothing
end

function _ruo2_resolve_location(indexed::SourceFile)
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

function parse_device_info(::RuO2Project, indexed::SourceFile)
    location = _ruo2_resolve_location(indexed)
    location === nothing && error("Could not resolve RuO2 device path from $(indexed.filepath)")
    return DeviceInfo(location)
end

function detect_kind(::RuO2Project, filename::String)::Symbol
    lower = lowercase(filename)
    if occursin("pund_fatigue", lower) || occursin("pund fatigue", lower)
        return :pund_fatigue
    elseif occursin("pund", lower) && occursin("wakeup", lower)
        return :pund_wakeup
    elseif occursin("fe pund", lower) || occursin("fepund", lower) || occursin("_pund", lower)
        return :pund
    elseif occursin("cvsweep", lower)
        return :cvsweep
    elseif occursin("tlm_4p", lower) || occursin("tlm", lower)
        return :tlm4p
    elseif occursin("i_v sweep", lower) || occursin("iv sweep", lower) ||
           occursin("ivsweep", lower) || occursin("fourterminaliv", lower)
        return :iv
    elseif occursin("break", lower) || occursin("breakdown", lower)
        return :breakdown
    else
        return :unknown
    end
end

function parse_measurement_parameters(indexed::SourceFile, kind::Symbol)
    if kind === :pund || kind === :pund_wakeup || kind === :pund_fatigue
        return Dict{Symbol,Any}(
            :wakeup_count => 0,
            :wakeup_f => NaN,
            :wakeup_V => NaN,
            :fatigue_count => 0,
            :fatigue_f => NaN,
            :fatigue_V => NaN,
        )
    end
    return Dict{Symbol,Any}()
end

function interpret_file(project::RuO2Project, indexed::SourceFile; should_cancel::Union{Nothing,Function}=nothing)::Vector{MeasurementItem}
    kind = detect_kind(project, indexed.filename)
    kind == :unknown && return MeasurementItem[]
    if kind == :cvsweep && !cv_sweep_has_schema(indexed.filepath)
        @warn "Ignoring unsupported RuO2 CVSweep file" path = indexed.filepath
        return MeasurementItem[]
    end
    device_info = parse_device_info(project, indexed)

    # assign measurement parameters (not stats)
    # that are ~directly stored in the file
    params = parse_measurement_parameters(indexed, kind)

    title = build_clean_title(project, indexed.filename, kind, device_info, indexed.header_summary)
    base = MeasurementItem(
        filepath=indexed.filepath,
        kind=kind,
        device_path=copy(device_info.location),
        timestamp=indexed.timestamp,
        parameters=params,
        title=title,
    )

    expanded = _ruo2_expand_multi_device_item(base)
    if kind == :pund_fatigue
        return vcat([_ruo2_expand_pund_fatigue_item(item; should_cancel=should_cancel) for item in expanded]...)
    elseif kind == :pund_wakeup
        return vcat([_ruo2_expand_pund_wakeup_item(item) for item in expanded]...)
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
    return [MeasurementItem(item;
        unique_id="$(item.filepath)#device=$p",
        device_path=vcat(loc[1:end-1], [p]),
        title=replace(item.title, dev => p),
    ) for p in parts]
end

function _ruo2_expand_pund_fatigue_item(item::MeasurementItem; should_cancel::Union{Nothing,Function}=nothing)::Vector{MeasurementItem}
    item.kind == :pund_fatigue || return [item]
    cycles, _ = _ruo2_scan_fatigue_file(item.filepath; should_cancel=should_cancel)
    isempty(cycles) && return MeasurementItem[]
    header = _ruo2_read_pund_header(item.filepath)
    fatigue_f = parse(Float64, first(split(header[:fatigue_freq], ',')))
    fatigue_V = parse(Float64, first(split(header[:vmax], ',')))
    return [MeasurementItem(item;
        unique_id="$(item.filepath)#fatigue_count=$(Int(c))",
        kind=:pund,
        parameters=Dict{Symbol,Any}(
            :wakeup_count => 0,
            :wakeup_f => NaN,
            :wakeup_V => NaN,
            :fatigue_count => c,
            :fatigue_f => fatigue_f,
            :fatigue_V => fatigue_V,
        ),
        title=item.title * " cycle $c",
    ) for c in cycles]
end

function _ruo2_scan_pund_wakeup(filepath::AbstractString)
    amplitudes = Float64[]
    pulse_type = :both
    open(filepath, "r") do io
        for line in eachline(io)
            if startswith(line, '#')
                m = match(r"read_pulse_type:\s*(\w+)", line)
                m !== nothing && (pulse_type = Symbol(m.captures[1]))
                continue
            end
            parts = split(line, ',')
            length(parts) >= 4 || continue
            v = tryparse(Float64, parts[1])
            v === nothing && continue
            v in amplitudes || push!(amplitudes, v)
        end
    end
    reps_per_amplitude = read_pund_wakeup_reps(basename(filepath), dirname(filepath))
    return (amplitudes=amplitudes, reps_per_amplitude=reps_per_amplitude, pulse_type=pulse_type)
end

"""Read the `#   key: value` parameter lines from a RuO2 PUND procedure file."""
function _ruo2_read_pund_header(filepath::AbstractString)
    header = Dict{Symbol,String}()
    open(filepath, "r") do io
        for raw_line in eachline(io)
            line = strip(raw_line)
            startswith(line, '#') || break
            m = match(r"^#\s+([^:]+):\s*(.*)$", line)
            m === nothing && continue
            header[Symbol(strip(m.captures[1]))] = strip(m.captures[2])
        end
    end
    return header
end

function _ruo2_expand_pund_wakeup_item(item::MeasurementItem)::Vector{MeasurementItem}
    item.kind == :pund_wakeup || return [item]
    info = _ruo2_scan_pund_wakeup(item.filepath)
    isempty(info.amplitudes) && return MeasurementItem[]
    header = _ruo2_read_pund_header(item.filepath)
    wakeup_count = parse(Float64, first(split(header[:fatigue_count], ',')))
    wakeup_f = parse(Float64, first(split(header[:fatigue_freq], ',')))

    segments = info.pulse_type === :pund ? [(:wakeup_pund, "PUND")] :
               info.pulse_type === :pn   ? [(:wakeup_pn,   "PN")]   :
                                           [(:wakeup_pn, "PN"), (:wakeup_pund, "PUND")]

    device_label = join(item.device_path[2:end], "_")
    date_str = item.timestamp === nothing ? "" : Dates.format(item.timestamp, "yyyy-mm-dd")

    result = MeasurementItem[]
    for amp in info.amplitudes
        amp_str = "$(amp)V"
        n_reps = get(info.reps_per_amplitude, amp, 1)
        for rep in 1:n_reps
            rep_suffix = n_reps > 1 ? " rep $rep" : ""
            for (kind, seg_label) in segments
                params = Dict{Symbol,Any}(
                    :wakeup_count => wakeup_count,
                    :wakeup_f => wakeup_f,
                    :wakeup_V => amp,
                    :fatigue_count => 0,
                    :fatigue_f => NaN,
                    :fatigue_V => NaN,
                )
                title = join(filter(!isempty, ["Wakeup", device_label, date_str, amp_str, seg_label * rep_suffix]), " ")
                push!(result, MeasurementItem(item;
                    unique_id="$(item.filepath)#wakeup_V=$(amp),rep=$(rep),kind=$(kind)",
                    kind=kind,
                    parameters=params,
                    title=title,
                ))
            end
        end
    end
    return result
end
