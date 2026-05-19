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
