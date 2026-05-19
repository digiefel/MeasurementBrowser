using DataLoader: read_pund_wakeup_amplitude, read_pund_wakeup_reps, cv_sweep_has_schema

# function device_path_label(::RuO2Project, device_info::DeviceInfo)
#     length(device_info.location) <= 1 && return join(device_info.location, "_")
#     return join(device_info.location[2:end], "_")
# end

"""Resolve the RuO2 device path for one indexed CSV file."""
function parse_device_info(::RuO2Project, file::SourceFile)
    location = _ruo2_resolve_location(file)
    location === nothing && error("Could not resolve RuO2 device path from $(file.filepath)")
    return DeviceInfo(location)
end

"""Classify a RuO2 filename into the measurement kind used by the browser."""
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

"""Create the six PUND measurement parameters for one logical measurement."""
function pund_measurement_parameters(;
    wakeup_count=0,
    wakeup_f=NaN,
    wakeup_V=NaN,
    fatigue_count=0,
    fatigue_f=NaN,
    fatigue_V=NaN,
)
    return Dict{Symbol,Any}(
        :wakeup_count => wakeup_count,
        :wakeup_f => wakeup_f,
        :wakeup_V => wakeup_V,
        :fatigue_count => fatigue_count,
        :fatigue_f => fatigue_f,
        :fatigue_V => fatigue_V,
    )
end

"""Interpret one RuO2 CSV file into the logical measurements shown by the browser."""
function interpret_file(project::RuO2Project, file::SourceFile; should_cancel::Union{Nothing,Function}=nothing)::Vector{MeasurementInfo}
    kind = detect_kind(project, file.filename)
    kind == :unknown && return MeasurementInfo[]
    if kind == :cvsweep && !cv_sweep_has_schema(file.filepath)
        @warn "Ignoring unsupported RuO2 CVSweep file" path = file.filepath
        return MeasurementInfo[]
    end
    device_info = parse_device_info(project, file)

    params = kind === :pund ? pund_measurement_parameters() : Dict{Symbol,Any}()
    stats = kind === :pund ? compute_pund_stats(file.filepath, params) : Dict{Symbol,Any}()

    title = build_clean_title(project, file.filename, kind, device_info, file.header_summary)
    base = MeasurementInfo(
        filepath=file.filepath,
        measurement_kind=kind,
        device_info=DeviceInfo(copy(device_info.location)),
        timestamp=file.timestamp,
        parameters=params,
        stats=stats,
        clean_title=title,
    )

    expanded = _ruo2_expand_multi_device_item(base)
    if kind == :pund_fatigue
        return vcat([_ruo2_expand_pund_fatigue_item(item; should_cancel=should_cancel) for item in expanded]...)
    elseif kind == :pund_wakeup
        return vcat([_ruo2_expand_pund_wakeup_item(item) for item in expanded]...)
    end
    return expanded
end

"""Return the device token embedded after a known chip/site/subsite filename prefix."""
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

"""Resolve a RuO2 device path from a nested chip/site/subsite/device directory layout."""
function _ruo2_nested_location(file::SourceFile)
    parent = dirname(file.filepath)
    tail = String[]
    while true
        name = basename(parent)
        if occursin(REGEX_RUO2_CHIP_DIR, name) && length(tail) == 3
            site, subsite, device_dir = reverse(tail)
            location = [name, site, subsite, device_dir]
            filename_device = _ruo2_prefixed_filename_device(file.filename, location[1:3])
            if filename_device === nothing
                return location
            end
            if startswith(device_dir, "RuO2")
                location[end] = filename_device
                return location
            end
            device_dir == filename_device || error(
                "RuO2 directory device '$device_dir' does not match filename device " *
                "'$filename_device' for $(file.filepath)",
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

"""Split a RuO2 Device ID string into chip, site, subsite, and device path parts."""
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

"""Read and validate the optional `Device ID` header from a RuO2 CSV file."""
function _ruo2_header_identifier(file::SourceFile)
    raw = get(file.header_summary, :device_id, nothing)
    raw isa AbstractString || return nothing
    stripped = strip(String(raw))
    isempty(stripped) && return nothing
    match_obj = match(REGEX_RUO2_IDENTIFIER, stripped)
    match_obj === nothing && error("Invalid RuO2 Device ID header '$stripped' in $(file.filepath)")
    return String(match_obj.captures[1])
end

"""Extract a RuO2 device identifier from supported filename conventions."""
function _ruo2_filename_identifier(filename::AbstractString)
    s = String(filename)
    m = match(REGEX_RUO2_BRACKET_IDENTIFIER, s)
    m !== nothing && return String(m.captures[1])
    m = match(REGEX_RUO2_TIMESTAMP_FILENAME, s)
    m !== nothing && return String(m.captures[1])
    return nothing
end

"""Choose the RuO2 device path from directory, header, or filename information."""
function _ruo2_resolve_location(file::SourceFile)
    nested = _ruo2_nested_location(file)
    header_identifier = _ruo2_header_identifier(file)
    if nested !== nothing
        if header_identifier !== nothing
            header_location = _ruo2_location_from_identifier(header_identifier)
            header_location == nested || error(
                "RuO2 Device ID header '$header_identifier' does not match directory hierarchy " *
                "'$(join(nested, '/'))' for $(file.filepath)",
            )
        end
        return nested
    end

    header_identifier !== nothing && return _ruo2_location_from_identifier(header_identifier)

    filename_identifier = _ruo2_filename_identifier(file.filename)
    filename_identifier !== nothing && return _ruo2_location_from_identifier(filename_identifier)

    return nothing
end

"""Split paired-device breakdown measurements into one browser measurement per device."""
function _ruo2_expand_multi_device_item(measurement::MeasurementInfo)::Vector{MeasurementInfo}
    measurement.measurement_kind == :breakdown || return [measurement]
    dev = last(measurement.device_info.location)
    if (m = match(r"^([A-Z][0-9]+)([A-Z][0-9]+)$", dev)) === nothing
        return [measurement]
    end
    parts = m.captures
    loc = copy(measurement.device_info.location)
    return [MeasurementInfo(measurement;
        unique_id="$(measurement.filepath)#device=$p",
        device_info=DeviceInfo(
            vcat(loc[1:end-1], [p]),
            deepcopy(measurement.device_info.parameters),
        ),
        clean_title=replace(measurement.clean_title, dev => p),
    ) for p in parts]
end

"""Expand a PUND fatigue file into one PUND measurement per fatigue count."""
function _ruo2_expand_pund_fatigue_item(measurement::MeasurementInfo; should_cancel::Union{Nothing,Function}=nothing)::Vector{MeasurementInfo}
    measurement.measurement_kind == :pund_fatigue || return [measurement]
    cycles, _ = _ruo2_scan_fatigue_file(measurement.filepath; should_cancel=should_cancel)
    isempty(cycles) && return MeasurementInfo[]
    header = _ruo2_read_pund_header(measurement.filepath)
    fatigue_f = haskey(header, :fatigue_freq) ?
        parse(Float64, first(split(header[:fatigue_freq], ','))) :
        NaN
    fatigue_V = haskey(header, :vmax) ?
        parse(Float64, first(split(header[:vmax], ','))) :
        NaN
    fatigue_df = _load_ruo2_pund_fatigue_file(measurement.filepath; should_cancel=should_cancel)
    return [MeasurementInfo(measurement;
        unique_id="$(measurement.filepath)#fatigue_count=$(Int(c))",
        measurement_kind=:pund,
        parameters=pund_measurement_parameters(
            fatigue_count=c,
            fatigue_f=fatigue_f,
            fatigue_V=fatigue_V,
        ),
        stats=pund_stats_from_waveform(_select_pund_fatigue_cycle(fatigue_df, c)),
        clean_title=measurement.clean_title * " cycle $c",
    ) for c in cycles]
end

"""Scan a PUND wakeup file for readout amplitudes, repetitions, and readout kind."""
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

"""Expand a PUND wakeup file into PN/PUND readout measurements for each wakeup voltage."""
function _ruo2_expand_pund_wakeup_item(measurement::MeasurementInfo)::Vector{MeasurementInfo}
    measurement.measurement_kind == :pund_wakeup || return [measurement]
    info = _ruo2_scan_pund_wakeup(measurement.filepath)
    isempty(info.amplitudes) && return MeasurementInfo[]
    header = _ruo2_read_pund_header(measurement.filepath)
    wakeup_count = parse(Float64, first(split(header[:fatigue_count], ',')))
    wakeup_f = parse(Float64, first(split(header[:fatigue_freq], ',')))

    segments = info.pulse_type === :pund ? [(:wakeup_pund, "PUND")] :
               info.pulse_type === :pn   ? [(:wakeup_pn,   "PN")]   :
                                           [(:wakeup_pn, "PN"), (:wakeup_pund, "PUND")]

    device_label = join(measurement.device_info.location[2:end], "_")
    date_str = measurement.timestamp === nothing ? "" : Dates.format(measurement.timestamp, "yyyy-mm-dd")

    result = MeasurementInfo[]
    for amp in info.amplitudes
        amp_str = "$(amp)V"
        n_reps = get(info.reps_per_amplitude, amp, 1)
        for rep in 1:n_reps
            rep_suffix = n_reps > 1 ? " rep $rep" : ""
            for (kind, seg_label) in segments
                params = pund_measurement_parameters(
                    wakeup_count=wakeup_count,
                    wakeup_f=wakeup_f,
                    wakeup_V=amp,
                )
                df = read_pund_wakeup_amplitude(
                    basename(measurement.filepath),
                    dirname(measurement.filepath),
                    amp,
                    rep,
                )
                title = join(filter(!isempty, ["Wakeup", device_label, date_str, amp_str, seg_label * rep_suffix]), " ")
                push!(result, MeasurementInfo(measurement;
                    unique_id="$(measurement.filepath)#wakeup_V=$(amp),rep=$(rep),kind=$(kind)",
                    measurement_kind=kind,
                    parameters=params,
                    stats=pund_stats_from_waveform(df),
                    clean_title=title,
                ))
            end
        end
    end
    return result
end
