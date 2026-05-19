function parse_device_info(::TASEProject, filename::String)
    m = match(REGEX_TASE, filename)
    m === nothing && error("Unrecognized TASE filename format: $filename")
    chip, facet, device_type, device_id = String.(m.captures)
    return DeviceInfo([chip, facet, device_type, device_id])
end

function parse_device_info(project::TASEProject, indexed::SourceFile)
    return parse_device_info(project, indexed.filename)
end

detect_kind(::TASEProject, filename::String)::Symbol =
    match(REGEX_TASE, filename) !== nothing ? :four_terminal_iv : :unknown

function interpret_file(project::TASEProject, indexed::SourceFile; should_cancel::Union{Nothing,Function}=nothing)::Vector{MeasurementItem}
    match(REGEX_TASE, indexed.filename) === nothing && return MeasurementItem[]
    device_info = parse_device_info(project, indexed)
    kind = detect_kind(project, indexed.filename)
    kind == :unknown && return MeasurementItem[]
    params = Dict{Symbol,Any}()
    title = build_clean_title(project, indexed.filename, kind, device_info, indexed.header_summary)
    return [MeasurementItem(
        filepath=indexed.filepath,
        kind=kind,
        device_path=copy(device_info.location),
        timestamp=indexed.timestamp,
        parameters=params,
        title=title,
    )]
end
