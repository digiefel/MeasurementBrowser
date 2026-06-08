# function device_path_label(::RuO2Project, device_info::DeviceInfo)
#     length(device_info.location) <= 1 && return join(device_info.location, "_")
#     return join(device_info.location[2:end], "_")
# end

"""Resolve the RuO2 device path for one indexed CSV file."""
function parse_device_info(::RuO2Project, file::SourceFile)
    location = _ruo2_location_from_filename(file.filename)
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

"""True when the source CSV stores many PUND readouts indexed by fatigue count."""
function is_pund_fatigue_file(filepath::AbstractString)
    return detect_kind(RUO2_PROJECT, basename(String(filepath))) === :pund_fatigue
end

"""Parse filename/header settings before waveform data or device history are analyzed."""
function parse_measurement_parameters(file::SourceFile, kind::Symbol)
    params = Dict{Symbol,Any}()
    if (m = match(r"(\d+(?:\.\d+)?)K", file.filename)) !== nothing
        params[:temperature_K] = parse(Float64, m.captures[1])
    end
    header = file.header_summary
    if kind === :pund_fatigue
        params[:fatigue_f] = parse(Float64, get(header, "fatigue_freq", "NaN"))
        params[:fatigue_Vamp] = parse(Float64, get(header, "vmax", "NaN"))
        params[:fatigue_Vbase] = parse(Float64, get(header, "base_voltage", "NaN"))
    elseif kind === :pund_wakeup
        # the keys in the wakeup headers are still called fatigue_* for historical reasons
        params[:wakeup_f] = parse(Float64, get(header, "fatigue_freq", "NaN"))
        params[:wakeup_count] = parse(Float64, get(header, "fatigue_count", "NaN"))
    end
    return params
end

"""Interpret one RuO2 CSV file into the logical measurements shown by the browser."""
function interpret_file(project::RuO2Project, file::SourceFile)::Vector{MeasurementInfo}
    kind = detect_kind(project, file.filename)
    kind == :unknown && return MeasurementInfo[]
    if kind == :cvsweep && !cv_sweep_has_schema(file.filepath)
        @warn "Ignoring unsupported RuO2 CVSweep file" path = file.filepath
        return MeasurementInfo[]
    end
    device_info = parse_device_info(project, file)

    title = build_clean_title(project, file.filename, kind, device_info, file.header_summary)
    base = MeasurementInfo(
        filepath=file.filepath,
        measurement_kind=kind,
        device_info=DeviceInfo(copy(device_info.location)),
        timestamp=file.timestamp,
        parameters=parse_measurement_parameters(file, kind),
        clean_title=title,
    )

    expanded = _ruo2_expand_multi_device_item(base)
    if kind == :pund_fatigue
        return vcat([_ruo2_expand_pund_fatigue_item(item) for item in expanded]...)
    elseif kind == :pund_wakeup
        return vcat([_ruo2_expand_pund_wakeup_item(item, file.header_summary) for item in expanded]...)
    end
    return expanded
end

"""Extract the RuO2 device path from supported filename conventions."""
function _ruo2_location_from_filename(filename::AbstractString)
    s = String(filename)
    # Match old (and probably unused) format
    m = match(r"\[((?:RuO2)[^()\[\];\s]+)", s)
    identifier = m === nothing ? nothing : m.captures[1]
    if identifier === nothing
        m = match(r"^(RuO2.+?)_\d{8}_\d{6}_", s)
        identifier = m === nothing ? nothing : m.captures[1]
    end
    identifier === nothing && return nothing

    # match things like RuO2test_{A1}_{RomanNumeralSite}_{Subsite}_{Device}
    m = match(r"^(RuO2test_?[A-Z0-9]+)_([XVI]+)_(.+)_([A-Z][0-9]+(?:[A-Z][0-9]+)*)$", identifier)
    m !== nothing && return [String(m[1]), String(m[2]), replace(String(m[3]), "_" => ""), String(m[4])]

    parts = split(identifier, '_')
    length(parts) >= 4 || error("Invalid RuO2 filename identifier: $identifier")
    return [String(join(parts[1:end-3], "_")), String(parts[end-2]), String(parts[end-1]), String(parts[end])]
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
function _ruo2_expand_pund_fatigue_item(measurement::MeasurementInfo)::Vector{MeasurementInfo}
    measurement.measurement_kind == :pund_fatigue || return [measurement]
    cycles = unique(read_pund_fatigue_file(
        measurement.filepath,
    ).cycle)
    isempty(cycles) && return MeasurementInfo[]
    return [MeasurementInfo(measurement;
        unique_id="$(measurement.filepath)#fatigue_count=$(Int(c))",
        measurement_kind=:pund,
        parameters=merge(deepcopy(measurement.parameters), Dict{Symbol,Any}(:fatigue_idx => c)),
        clean_title=measurement.clean_title * " cycle $c",
    ) for c in cycles]
end

"""Expand a PUND wakeup file into PN/PUND readout measurements for each wakeup voltage."""
function _ruo2_expand_pund_wakeup_item(
    measurement::MeasurementInfo,
    header::Dict{String,String},
)::Vector{MeasurementInfo}
    measurement.measurement_kind == :pund_wakeup || return [measurement]
    amplitudes = parse.(Float64, strip.(split(header["vmax"], ',')))
    isempty(amplitudes) && return MeasurementInfo[]

    pulse_type = Symbol(header["read_pulse_type"])
    segments = pulse_type === :pund ? [(:wakeup_pund, "PUND")] :
               pulse_type === :pn   ? [(:wakeup_pn,   "PN")]   :
                                      [(:wakeup_pn, "PN"), (:wakeup_pund, "PUND")]

    device_label = join(measurement.device_info.location[2:end], "_")
    date_str = measurement.timestamp === nothing ? "" : Dates.format(measurement.timestamp, "yyyy-mm-dd")

    result = MeasurementInfo[]
    for amp in amplitudes
        amp_str = "$(amp)V"
        for (kind, seg_label) in segments
            params = merge(deepcopy(measurement.parameters), Dict{Symbol,Any}(:wakeup_V => amp))
            title = join(
                filter(!isempty, ["Wakeup", device_label, date_str, amp_str, seg_label]),
                " ",
            )
            push!(result, MeasurementInfo(measurement;
                unique_id="$(measurement.filepath)#wakeup_V=$(amp),kind=$(kind)",
                measurement_kind=kind,
                parameters=params,
                clean_title=title,
            ))
        end
    end
    return result
end
