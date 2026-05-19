function parse_device_info(::TASEProject, filename::String)
    m = match(REGEX_TASE, filename)
    m === nothing && error("Unrecognized TASE filename format: $filename")
    chip, facet, device_type, device_id = String.(m.captures)
    return DeviceInfo([chip, facet, device_type, device_id])
end

function parse_device_info(project::TASEProject, file::SourceFile)
    return parse_device_info(project, file.filename)
end

detect_kind(::TASEProject, filename::String)::Symbol =
    match(REGEX_TASE, filename) !== nothing ? :four_terminal_iv : :unknown

function interpret_file(project::TASEProject, file::SourceFile; should_cancel::Union{Nothing,Function}=nothing)::Vector{MeasurementInfo}
    match(REGEX_TASE, file.filename) === nothing && return MeasurementInfo[]
    device_info = parse_device_info(project, file)
    kind = detect_kind(project, file.filename)
    kind == :unknown && return MeasurementInfo[]
    params = Dict{Symbol,Any}()
    title = build_clean_title(project, file.filename, kind, device_info, file.header_summary)
    return [MeasurementInfo(
        filepath=file.filepath,
        measurement_kind=kind,
        device_info=DeviceInfo(copy(device_info.location)),
        timestamp=file.timestamp,
        parameters=params,
        clean_title=title,
    )]
end
