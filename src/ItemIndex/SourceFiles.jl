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

"""One physical source file discovered inside a file-backed data source."""
struct SourceFile <: AbstractDataSourceItem
    filepath::String
    filename::String
    timestamp::Union{DateTime,Nothing}
    fingerprint::FileFingerprint
end

Projects.source_item_id(file::SourceFile)::String = file.filepath
Projects.source_item_label(file::SourceFile)::String = file.filename
Projects.source_item_fingerprint(file::SourceFile)::FileFingerprint = file.fingerprint
Projects.source_item_path(file::SourceFile)::String = file.filepath
Projects.source_item_timestamp(file::SourceFile)::Union{DateTime,Nothing} = file.timestamp

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

"""Index one physical source file without interpreting its data items."""
function index_source_file(path::AbstractString)::SourceFile
    normalized = normpath(abspath(expanduser(String(path))))
    filename = basename(normalized)
    return SourceFile(
        normalized,
        filename,
        parse_timestamp(filename),
        file_fingerprint(normalized),
    )
end

"""Package sidecar files that are never candidate source files."""
const _SIDECAR_FILENAMES = Set([
    "device_info.txt",          # collection metadata
    "tags.txt", "notes.txt", "layout.txt",   # annotations
    "measurementbrowser.toml",  # saved project view
])

"""
Return whether a directory entry is a candidate source file.

The package is data-agnostic, so every visible file is a candidate and each project's `detect`
decides which recipe (if any) claims it — files no recipe detects are simply skipped. Hidden files
and the package's own sidecars (collection metadata, annotations, saved view) are never candidates.
"""
function is_source_filename(name::AbstractString)::Bool
    startswith(name, ".") && return false
    return !(lowercase(name) in _SIDECAR_FILENAMES)
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
