using Dates

const DEFAULT_DIRECTORY_METADATA_FILE = "metadata.txt"

"""
Directory-backed source scanned into `SourceFile` source items.
"""
mutable struct DirectorySource <: AbstractDataSource
    root_path::String
    recursive::Bool
    metadata_file::Union{Nothing,String}
    collection_parameter_entries::Union{Nothing,Dict{Tuple{Vararg{String}},Dict{Symbol,Any}}}
end

function DirectorySource(
    root_path::AbstractString;
    recursive::Bool=true,
    metadata_file::Union{Nothing,AbstractString}=DEFAULT_DIRECTORY_METADATA_FILE,
)::DirectorySource
    normalized = normpath(abspath(expanduser(String(root_path))))
    return DirectorySource(
        normalized,
        recursive,
        metadata_file === nothing ? nothing : String(metadata_file),
        nothing,
    )
end

"""
Fingerprint used to decide whether a physical source file has changed.
"""
struct FileFingerprint
    path::String
    size_bytes::Int64
    mtime_ns::Int64
end

function Base.:(==)(left::FileFingerprint, right::FileFingerprint)::Bool
    return left.path == right.path &&
        left.size_bytes == right.size_bytes &&
        left.mtime_ns == right.mtime_ns
end

"""One physical source file discovered inside a directory source."""
struct SourceFile <: AbstractDataSourceItem
    filepath::String
    filename::String
    timestamp::Union{DateTime,Nothing}
    fingerprint::FileFingerprint
end

Projects.source_item_id(file::SourceFile)::String = file.filepath
Projects.source_item_label(file::SourceFile)::String = file.filename
Projects.fingerprint(file::SourceFile)::FileFingerprint = file.fingerprint
Projects.source_item_path(file::SourceFile)::String = file.filepath
Projects.source_item_timestamp(file::SourceFile)::Union{DateTime,Nothing} = file.timestamp

Projects.source_id(source::DirectorySource)::String = source.root_path
Projects.source_label(source::DirectorySource)::String = basename(source.root_path)

"""Extract a source-item timestamp from the supported filename conventions."""
function parse_timestamp(filename::AbstractString)::Union{DateTime,Nothing}
    if (result = match(
        r"; (\d{4}-\d{2}-\d{2}) (\d{2})_(\d{2})_(\d{2})\]",
        filename,
    )) !== nothing
        date, hour, minute, second = result.captures
        return try
            DateTime("$date $hour:$minute:$second", "yyyy-mm-dd HH:MM:SS")
        catch
            nothing
        end
    end
    result = match(r"_(\d{4})(\d{2})(\d{2})_(\d{2})(\d{2})(\d{2})_", filename)
    result === nothing && return nothing
    year, month, day, hour, minute, second = result.captures
    return try
        DateTime("$year-$month-$day $hour:$minute:$second", "yyyy-mm-dd HH:MM:SS")
    catch
        nothing
    end
end

"""Fingerprint an already-normalized absolute path (no redundant re-normalization)."""
function _file_fingerprint(normalized::String)::FileFingerprint
    stat_info = stat(normalized)
    return FileFingerprint(
        normalized,
        Int64(stat_info.size),
        Int64(stat_info.mtime * 1_000_000_000),
    )
end

file_fingerprint(path::AbstractString)::FileFingerprint =
    _file_fingerprint(normpath(abspath(expanduser(String(path)))))

function index_source_file(path::AbstractString)::SourceFile
    normalized = normpath(abspath(expanduser(String(path))))
    filename = basename(normalized)
    return SourceFile(
        normalized,
        filename,
        parse_timestamp(filename),
        _file_fingerprint(normalized),
    )
end

"""DirectorySource sidecars that are never candidate source files."""
const DIRECTORY_SOURCE_SIDECARS = Set([
    DEFAULT_DIRECTORY_METADATA_FILE,
    "tags.txt",
    "notes.txt",
    "layout.txt",
    "measurementbrowser.toml",
])

function is_source_filename(name::AbstractString)::Bool
    startswith(name, ".") && return false
    return !(lowercase(name) in DIRECTORY_SOURCE_SIDECARS)
end

function collect_source_files(source::DirectorySource)::Vector{SourceFile}
    source_files = SourceFile[]
    if source.recursive
        for (root, _, names) in walkdir(source.root_path)
            append_source_files!(source_files, root, names)
        end
    else
        names = readdir(source.root_path)
        append_source_files!(source_files, source.root_path, names)
    end
    return source_files
end

function append_source_files!(
    source_files::Vector{SourceFile},
    root::AbstractString,
    names::Vector{String},
)::Nothing
    for name in names
        check_cancel()
        is_source_filename(name) || continue
        path = joinpath(root, name)
        isfile(path) || continue
        push!(source_files, index_source_file(path))
    end
    return nothing
end

function parse_metadata_value(text::AbstractString)::Any
    value = strip(text)
    isempty(value) && return nothing
    lowercase_value = lowercase(value)
    lowercase_value in ("true", "t", "yes", "y", "1") && return true
    lowercase_value in ("false", "f", "no", "n", "0") && return false
    for type in (Int, Float64)
        parsed = tryparse(type, value)
        parsed === nothing || return parsed
    end
    for format in (dateformat"yyyy-mm-dd", dateformat"yyyy/mm/dd")
        parsed = tryparse(Date, value, format)
        parsed === nothing || return parsed
    end
    for format in (
        dateformat"yyyy-mm-dd HH:MM:SS",
        dateformat"yyyy-mm-ddTHH:MM:SS",
    )
        parsed = tryparse(DateTime, value, format)
        parsed === nothing || return parsed
    end
    return value
end

function load_collection_parameter_entries(
    path::AbstractString,
)::Dict{Tuple{Vararg{String}},Dict{Symbol,Any}}
    lines = readlines(path)
    isempty(lines) && return Dict{Tuple{Vararg{String}},Dict{Symbol,Any}}()
    header = strip.(split(first(lines), ','))
    length(header) >= 2 ||
        return Dict{Tuple{Vararg{String}},Dict{Symbol,Any}}()
    parameter_names = header[2:end]
    entries = Dict{Tuple{Vararg{String}},Dict{Symbol,Any}}()
    for line in Iterators.drop(lines, 1)
        isempty(strip(line)) && continue
        cells = split(line, ',')
        path_text = strip(first(cells))
        isempty(path_text) && continue
        parameters = Dict{Symbol,Any}()
        for (offset, name) in enumerate(parameter_names)
            index = offset + 1
            index > length(cells) && continue
            value = parse_metadata_value(cells[index])
            value === nothing || (parameters[Symbol(strip(name))] = value)
        end
        entries[Tuple(filter(!isempty, split(path_text, '/')))] = parameters
    end
    return entries
end

function collection_parameter_file_path(source::DirectorySource)::Union{Nothing,String}
    source.metadata_file === nothing && return nothing
    return joinpath(source.root_path, source.metadata_file)
end

function has_collection_parameter_file(source::DirectorySource)::Bool
    path = collection_parameter_file_path(source)
    return path !== nothing && isfile(path)
end

function load_scan_collection_parameters(
    source::DirectorySource,
)::Union{Nothing,Dict{Tuple{Vararg{String}},Dict{Symbol,Any}}}
    path = collection_parameter_file_path(source)
    return path !== nothing && isfile(path) ? load_collection_parameter_entries(path) : nothing
end

function matching_collection_parameters(
    entries::Dict{Tuple{Vararg{String}},Dict{Symbol,Any}},
    location::AbstractVector{<:AbstractString},
)::Dict{Symbol,Any}
    merged = Dict{Symbol,Any}()
    for width in eachindex(location)
        for start in 1:(length(location) - width + 1)
            parameters = get(
                entries,
                Tuple(location[start:(start + width - 1)]),
                nothing,
            )
            parameters === nothing || merge!(merged, parameters)
        end
    end
    return merged
end

function Projects.source_items(source::DirectorySource)::Vector{SourceFile}
    source.collection_parameter_entries = load_scan_collection_parameters(source)
    return collect_source_files(source)
end

function Projects.collection_parameters(
    source::DirectorySource,
    collection_path::AbstractVector{<:AbstractString},
)::Dict{Symbol,Any}
    source.collection_parameter_entries === nothing && return Dict{Symbol,Any}()
    return matching_collection_parameters(source.collection_parameter_entries, collection_path)
end

Projects.has_collection_parameters(source::DirectorySource)::Bool =
    source.collection_parameter_entries !== nothing

function Workspace.open_workspace(
    project::Project,
    root_path::AbstractString;
    recursive::Bool=true,
    metadata_file::Union{Nothing,AbstractString}=DEFAULT_DIRECTORY_METADATA_FILE,
    profile_internal::Bool=Profiling.environment_flag("MB_PROFILE_INTERNAL"),
    profile_cpu::Bool=Profiling.environment_flag("MB_PROFILE_CPU"),
    profile_output::Union{Nothing,AbstractString}=
        Profiling.environment_path("MB_PROFILE_OUTPUT"),
    crash_trace::Union{Nothing,AbstractString}=
        Profiling.environment_path("MB_CRASH_TRACE"),
)::Workspace.Workspace
    return Workspace.open_workspace(
        project,
        DirectorySource(root_path; recursive, metadata_file);
        profile_internal,
        profile_cpu,
        profile_output,
        crash_trace,
    )
end
