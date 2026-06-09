using Dates

struct FileFingerprint
    path::String
    size_bytes::Int64
    mtime_ns::Int64
end

"""Return true when two fingerprints describe the same source-file version."""
function Base.:(==)(left::FileFingerprint, right::FileFingerprint)::Bool
    return left.path == right.path &&
        left.size_bytes == right.size_bytes &&
        left.mtime_ns == right.mtime_ns
end

struct SourceFile
    unique_id::String
    filepath::String
    filename::String
    timestamp::Union{DateTime,Nothing}
    header_summary::Dict{String,String}
    fingerprint::FileFingerprint
    measurements::Vector
end

function SourceFile(file::SourceFile, measurements::Vector)
    return SourceFile(
        file.unique_id,
        file.filepath,
        file.filename,
        file.timestamp,
        file.header_summary,
        file.fingerprint,
        measurements,
    )
end
