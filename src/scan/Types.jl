using Dates

struct FileFingerprint
    path::String
    size_bytes::Int64
    mtime_ns::Int64
end

struct SourceFile
    id::String
    filepath::String
    filename::String
    timestamp::Union{DateTime,Nothing}
    header_summary::Dict{Symbol,Any}
    fingerprint::FileFingerprint
    measurements::Vector
end

struct MeasurementItem
    id::String
    source_file_id::String
    filepath::String
    kind::Symbol
    device_path::Vector{String}
    timestamp::Union{DateTime,Nothing}
    device_parameters::Dict{Symbol,Any}
    parameters::Dict{Symbol,Any}
    title::String
end

source_path(path::AbstractString) = normpath(abspath(expanduser(String(path))))
file_id(path::AbstractString) = source_path(path)

function source_file_with_measurements(source::SourceFile, measurements::Vector)
    return SourceFile(
        source.id,
        source.filepath,
        source.filename,
        source.timestamp,
        source.header_summary,
        source.fingerprint,
        measurements,
    )
end

function item_id(
    source_file_id::AbstractString;
    cycle::Union{Nothing,Integer}=nothing,
    split::Union{Nothing,AbstractString}=nothing,
)
    base = String(source_file_id)
    if cycle !== nothing
        return base * "#cycle=$(Int(cycle))"
    elseif split !== nothing
        return base * "#split=$(String(split))"
    end
    return base
end
