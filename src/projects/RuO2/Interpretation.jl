using DataLoader: read_wakeup_summary

function parse_device_info(::RuO2Project, filename::String)
    if (m = match(REGEX_DEVICE_NEW, filename)) !== nothing
        caps = filter(!isnothing, collect(m.captures))
        return DeviceInfo(String.(caps))
    elseif (m = match(REGEX_DEVICE, filename)) !== nothing
        caps = filter(!isnothing, collect(m.captures))
        return DeviceInfo(String.(caps))
    end
    error("Unrecognized RuO2 device filename format: $filename")
end

function detect_kind(::RuO2Project, filename::String)::Symbol
    lower = lowercase(filename)
    if occursin("pund_fatigue", lower) || occursin("pund fatigue", lower)
        return :pund_fatigue
    elseif occursin("fe pund", lower) || occursin("fepund", lower)
        return :pund
    elseif occursin("i_v sweep", lower) || occursin("iv sweep", lower)
        return :iv
    elseif occursin("tlm_4p", lower) || occursin("tlm", lower)
        return :tlm4p
    elseif occursin("break", lower) || occursin("breakdown", lower)
        return :breakdown
    elseif occursin("wakeup", lower)
        return :wakeup
    else
        return :unknown
    end
end

function interpret_file(::RuO2Project, indexed::IndexedCsvFile; should_cancel::Union{Nothing,Function}=nothing)::Vector{MeasurementItem}
    occursin("ruo2test_", lowercase(indexed.filename)) || return MeasurementItem[]
    device_info = try
        parse_device_info(RUO2_PROJECT, indexed.filename)
    catch
        return MeasurementItem[]
    end
    kind = detect_kind(RUO2_PROJECT, indexed.filename)
    kind == :unknown && return MeasurementItem[]

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
