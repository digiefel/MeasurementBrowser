function parse_device_info(::TASEProject, filename::String)
    m = match(REGEX_TASE, filename)
    m === nothing && error("Unrecognized TASE filename format: $filename")
    chip, facet, device_type, device_id = String.(m.captures)
    return DeviceInfo([chip, facet, device_type, device_id])
end

function parse_device_info(project::TASEProject, indexed::IndexedCsvFile)
    return parse_device_info(project, indexed.filename)
end

detect_kind(::TASEProject, filename::String)::Symbol =
    match(REGEX_TASE, filename) !== nothing ? :four_terminal_iv : :unknown

function interpret_file(::TASEProject, indexed::IndexedCsvFile; should_cancel::Union{Nothing,Function}=nothing)::Vector{MeasurementItem}
    match(REGEX_TASE, indexed.filename) === nothing && return MeasurementItem[]
    device_info = parse_device_info(TASE_PROJECT, indexed)
    kind = detect_kind(TASE_PROJECT, indexed.filename)
    kind == :unknown && return MeasurementItem[]
    params = parse_parameters(indexed.filename)
    title = build_clean_title(TASE_PROJECT, indexed.filename, kind, device_info, indexed.header_summary)
    return [MeasurementItem(
        item_id(indexed.id),
        indexed.id,
        indexed.filepath,
        kind,
        copy(device_info.location),
        indexed.timestamp,
        Dict{Symbol,Any}(),
        params,
        title,
    )]
end
