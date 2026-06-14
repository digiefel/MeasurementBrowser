"""
Fingerprint used to decide whether a physical source file has changed.
"""
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

"""
One physical source file and the logical measurements interpreted from it.
"""
struct SourceFile
    unique_id::String
    filepath::String
    filename::String
    timestamp::Union{DateTime,Nothing}
    fingerprint::FileFingerprint
    measurements::Vector{MeasurementInfo}
end

"""Copy indexed file metadata and attach its interpreted measurements."""
function SourceFile(
    file::SourceFile,
    measurements::Vector{MeasurementInfo},
)::SourceFile
    return SourceFile(
        file.unique_id,
        file.filepath,
        file.filename,
        file.timestamp,
        file.fingerprint,
        measurements,
    )
end

"""Read the current fingerprint of one physical source file."""
function file_fingerprint(path::AbstractString)::FileFingerprint
    normalized = normpath(abspath(expanduser(String(path))))
    stat_info = stat(normalized)
    return FileFingerprint(
        normalized,
        Int64(stat_info.size),
        Int64(stat_info.mtime * 1_000_000_000),
    )
end

"""Index one physical source file without reading its measurement data."""
function index_source_file(path::AbstractString)::SourceFile
    normalized = normpath(abspath(expanduser(String(path))))
    filename = basename(normalized)
    return SourceFile(
        normalized,
        normalized,
        filename,
        parse_timestamp(filename),
        file_fingerprint(normalized),
        MeasurementInfo[],
    )
end

"""Return whether a directory entry is a visible CSV source file."""
function is_source_filename(name::AbstractString)::Bool
    return endswith(lowercase(name), ".csv") && !startswith(name, ".")
end

"""Collect every indexed source file below a root."""
function collect_source_files(
    root_path::AbstractString;
    on_file::Union{Nothing,Function}=nothing,
)::Vector{SourceFile}
    source_files = SourceFile[]
    for (root, _, names) in walkdir(root_path)
        for name in names
            check_cancel()
            is_source_filename(name) || continue
            file = index_source_file(joinpath(root, name))
            push!(source_files, file)
            on_file !== nothing && on_file(file, length(source_files))
        end
    end
    return source_files
end
