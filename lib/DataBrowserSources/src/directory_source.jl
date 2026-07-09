using Dates
using BetterFileWatching
using CancellationTokens: CancellationToken, CancellationTokenSource, OperationCanceledException, cancel, get_token, is_cancellation_requested

import DataBrowserAPI
using DataBrowserAPI:
    AbstractDataSource,
    AbstractDataSourceItem,
    Project,
    SourceChanges,
    SourceError
import DataBrowserAPI:
    close_source!,
    collection_metadata,
    fingerprint,
    has_collection_metadata,
    open_source,
    source_id,
    source_item_id,
    source_item_label,
    source_item_noun,
    source_item_path,
    source_item_timestamp,
    source_items,
    source_label,
    source_open_options,
    watch_source

const DEFAULT_DIRECTORY_METADATA_FILE = "metadata.txt"

"""
Fingerprint used to decide whether a physical source file has changed.
"""
struct FileFingerprint
    path::String
    size_bytes::Int64
    mtime_ns::Int64
end

"""
Directory-backed source scanned into `SourceFile` source items.
"""
mutable struct DirectorySource <: AbstractDataSource
    root_path::String
    recursive::Bool
    metadata_file::Union{Nothing,String}
    collection_metadata_entries::Dict{Tuple{Vararg{String}},Dict{Symbol,Any}}
    has_metadata::Bool
    metadata_lock::ReentrantLock
    watcher_task::Union{Nothing,Task}
    watcher_cancel::Union{Nothing,CancellationTokenSource}
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
        Dict{Tuple{Vararg{String}},Dict{Symbol,Any}}(),
        false,
        ReentrantLock(),
        nothing,
        nothing,
    )
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

source_item_id(file::SourceFile)::String = file.filepath
source_item_label(file::SourceFile)::String = file.filename
fingerprint(file::SourceFile)::FileFingerprint = file.fingerprint
source_item_path(file::SourceFile)::String = file.filepath
source_item_timestamp(file::SourceFile)::Union{DateTime,Nothing} = file.timestamp

source_id(source::DirectorySource)::String = source.root_path
source_label(source::DirectorySource)::String = basename(source.root_path)
source_item_noun(::DirectorySource)::String = "source files"

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

function file_fingerprint(path::AbstractString; normalized::Bool=false)::FileFingerprint
    normalized_path = normalized ? String(path) : normpath(abspath(expanduser(String(path))))
    stat_info = stat(normalized_path)
    return FileFingerprint(
        normalized_path,
        Int64(stat_info.size),
        Int64(stat_info.mtime * 1_000_000_000),
    )
end

function index_source_file(path::AbstractString)::SourceFile
    normalized = normpath(abspath(expanduser(String(path))))
    filename = basename(normalized)
    return SourceFile(
        normalized,
        filename,
        parse_timestamp(filename),
        file_fingerprint(normalized; normalized=true),
    )
end

"""DirectorySource sidecars that are never candidate source files."""
const DIRECTORY_SOURCE_SIDECARS = Set([
    DEFAULT_DIRECTORY_METADATA_FILE,
    "tags.txt",
    "notes.txt",
    "layout.txt",
    "databrowser.toml",
])

function is_source_filename(name::AbstractString)::Bool
    startswith(name, ".") && return false
    return !(lowercase(name) in DIRECTORY_SOURCE_SIDECARS)
end

function collect_source_files(
    source::DirectorySource;
    on_progress::Union{Nothing,Function}=nothing,
    cancel_token::CancellationToken,
)::Vector{SourceFile}
    source_files = SourceFile[]
    metadata_path = collection_metadata_file_path(source)
    found = Ref(0)
    if source.recursive
        for (root, _, names) in walkdir(source.root_path)
            append_source_files!(
                source_files, root, names, metadata_path;
                on_progress, found, cancel_token,
            )
        end
    else
        names = readdir(source.root_path)
        append_source_files!(
            source_files, source.root_path, names, metadata_path;
            on_progress, found, cancel_token,
        )
    end
    on_progress === nothing || on_progress(found[])
    return source_files
end

function append_source_files!(
    source_files::Vector{SourceFile},
    root::AbstractString,
    names::Vector{String},
    metadata_path::Union{Nothing,String}=nothing,
    ;
    on_progress::Union{Nothing,Function}=nothing,
    found::Base.RefValue{Int}=Ref(0),
    cancel_token::CancellationToken,
)::Nothing
    for name in names
        is_cancellation_requested(cancel_token) &&
            throw(OperationCanceledException(cancel_token))
        is_source_filename(name) || continue
        path = joinpath(root, name)
        metadata_path !== nothing && normpath(path) == metadata_path && continue
        isfile(path) || continue
        push!(source_files, index_source_file(path))
        found[] += 1
        on_progress !== nothing && found[] % 256 == 0 && on_progress(found[])
    end
    return nothing
end

function parse_metadata_value(text::AbstractString)::Any
    value = strip(text)
    isempty(value) && return nothing
    lowercase_value = lowercase(value)
    lowercase_value in ("true", "t", "yes", "y") && return true
    lowercase_value in ("false", "f", "no", "n") && return false
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

function load_collection_metadata_entries(
    path::AbstractString,
)::Dict{Tuple{Vararg{String}},Dict{Symbol,Any}}
    lines = readlines(path)
    isempty(lines) && return Dict{Tuple{Vararg{String}},Dict{Symbol,Any}}()
    header = strip.(split(first(lines), ','))
    length(header) >= 2 ||
        throw(ArgumentError(
            "Invalid DirectorySource metadata file '$path': expected a path column and " *
            "at least one metadata column",
        ))
    metadata_names = header[2:end]
    entries = Dict{Tuple{Vararg{String}},Dict{Symbol,Any}}()
    for line in Iterators.drop(lines, 1)
        isempty(strip(line)) && continue
        cells = split(line, ',')
        path_text = strip(first(cells))
        isempty(path_text) && continue
        metadata = Dict{Symbol,Any}()
        for (offset, name) in enumerate(metadata_names)
            index = offset + 1
            index > length(cells) && continue
            value = parse_metadata_value(cells[index])
            value === nothing || (metadata[Symbol(strip(name))] = value)
        end
        entries[Tuple(filter(!isempty, split(path_text, '/')))] = metadata
    end
    return entries
end

function collection_metadata_file_path(source::DirectorySource)::Union{Nothing,String}
    source.metadata_file === nothing && return nothing
    return normpath(joinpath(source.root_path, source.metadata_file))
end

function load_collection_metadata!(source::DirectorySource)::Bool
    path = collection_metadata_file_path(source)
    present = path !== nothing && isfile(path)
    entries = present ?
        load_collection_metadata_entries(path::String) :
        Dict{Tuple{Vararg{String}},Dict{Symbol,Any}}()
    return lock(source.metadata_lock) do
        changed = source.has_metadata != present ||
            source.collection_metadata_entries != entries
        source.collection_metadata_entries = entries
        source.has_metadata = present
        changed
    end
end

function matching_collection_metadata(
    entries::Dict{Tuple{Vararg{String}},Dict{Symbol,Any}},
    location::AbstractVector{<:AbstractString},
)::Dict{Symbol,Any}
    merged = Dict{Symbol,Any}()
    for width in eachindex(location)
        for start in 1:(length(location) - width + 1)
            metadata = get(
                entries,
                Tuple(location[start:(start + width - 1)]),
                nothing,
            )
            metadata === nothing || merge!(merged, metadata)
        end
    end
    return merged
end

function source_items(
    source::DirectorySource;
    on_progress::Union{Nothing,Function}=nothing,
    cancel_token::CancellationToken,
)::Vector{SourceFile}
    return collect_source_files(source; on_progress, cancel_token)
end

function has_collection_metadata(source::DirectorySource)::Bool
    return lock(source.metadata_lock) do
        source.has_metadata
    end
end

function collection_metadata(
    source::DirectorySource,
    collection_path::AbstractVector{<:AbstractString},
)::Dict{Symbol,Any}
    return lock(source.metadata_lock) do
        matching_collection_metadata(source.collection_metadata_entries, collection_path)
    end
end

function open_source(source::DirectorySource)::DirectorySource
    isdir(source.root_path) || throw(ArgumentError(
        "Cannot open DirectorySource '$(source.root_path)': directory does not exist",
    ))
    load_collection_metadata!(source)
    return source
end

function watch_source(
    source::DirectorySource,
    on_change::Function;
    cancel_token::CancellationToken,
)::Union{Nothing,Task}
    source.watcher_task === nothing || error(
        "DirectorySource '$(source.root_path)' is already being watched",
    )
    cancel_source = CancellationTokenSource(cancel_token)
    watch_token = get_token(cancel_source)
    known = Dict(
        source_item_id(file) => fingerprint(file)
        for file in collect_source_files(source; cancel_token=watch_token)
    )
    ignore(rel::AbstractString)::Bool =
        ".git" in split(String(rel), '/') || (!source.recursive && occursin('/', String(rel)))
    task = Base.Threads.@spawn begin
        errored = false
        try
            watch_folder(source.root_path, watch_token; ignore) do _
                is_cancellation_requested(watch_token) && return
                metadata_changed = try
                    load_collection_metadata!(source)
                catch error
                    errored = true
                    on_change(SourceError(sprint(showerror, error)))
                    return
                end
                files = collect_source_files(source; cancel_token=watch_token)
                is_cancellation_requested(watch_token) && return
                current = Dict(source_item_id(file) => fingerprint(file) for file in files)
                upserts = SourceFile[
                    file for file in files
                    if !haskey(known, source_item_id(file)) ||
                       !isequal(known[source_item_id(file)], fingerprint(file))
                ]
                removals = String[id for id in keys(known) if !haskey(current, id)]
                isempty(upserts) && isempty(removals) && !metadata_changed && !errored && return
                known = current
                errored = false
                on_change(SourceChanges(upserts, removals; metadata_changed))
            end
        catch error
            error isa OperationCanceledException && return
            on_change(SourceError(sprint(showerror, error)))
        end
    end
    source.watcher_task = task
    source.watcher_cancel = cancel_source
    return task
end

source_open_options(source::DirectorySource)::NamedTuple =
    (; recursive=source.recursive, metadata_file=source.metadata_file)

function close_source!(source::DirectorySource)::Nothing
    source.watcher_cancel === nothing || cancel(source.watcher_cancel)
    task = source.watcher_task
    if task !== nothing
        wait(task)
    end
    source.watcher_task = nothing
    source.watcher_cancel = nothing
    return nothing
end
