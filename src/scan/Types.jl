using Dates

struct FileFingerprint
    path::String
    size_bytes::Int64
    mtime_ns::Int64
end

struct SourceFile
    unique_id::String
    filepath::String
    filename::String
    timestamp::Union{DateTime,Nothing}
    header_summary::Dict{Symbol,Any}
    fingerprint::FileFingerprint
    measurements::Vector
end

function SourceFile(source::SourceFile, measurements::Vector)
    return SourceFile(
        source.unique_id,
        source.filepath,
        source.filename,
        source.timestamp,
        source.header_summary,
        source.fingerprint,
        measurements,
    )
end

struct MeasurementItem
    unique_id::String
    filepath::String
    kind::Symbol
    device_path::Vector{String}
    timestamp::Union{DateTime,Nothing}
    device_parameters::Dict{Symbol,Any}
    parameters::Dict{Symbol,Any}
    stats::Dict{Symbol,Any}
    title::String
end

function MeasurementItem(
    item::Union{Nothing,MeasurementItem}=nothing;
    filepath=item === nothing ? nothing : item.filepath,
    kind=item === nothing ? nothing : item.kind,
    device_path=item === nothing ? nothing : deepcopy(item.device_path),
    timestamp=item === nothing ? nothing : item.timestamp,
    title=item === nothing ? nothing : item.title,
    unique_id=item === nothing ? filepath : item.unique_id,
    device_parameters=item === nothing ? Dict{Symbol,Any}() : deepcopy(item.device_parameters),
    parameters=item === nothing ? Dict{Symbol,Any}() : deepcopy(item.parameters),
    stats=item === nothing ? Dict{Symbol,Any}() : deepcopy(item.stats),
)
    filepath === nothing && error("MeasurementItem requires filepath")
    kind === nothing && error("MeasurementItem requires kind")
    device_path === nothing && error("MeasurementItem requires device_path")
    title === nothing && error("MeasurementItem requires title")
    unique_id === nothing && error("MeasurementItem requires unique_id")
    return MeasurementItem(
        String(unique_id),
        String(filepath),
        kind,
        device_path,
        timestamp,
        device_parameters,
        parameters,
        stats,
        title,
    )
end
