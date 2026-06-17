module ItemIndex

using Dates

import ..Projects
import ..Projects:
    AbstractDataItem,
    Project,
    compute_and_add_item_stats!,
    collection_path_label,
    interpret_file,
    kind_label,
    project_name,
    reset_scan_profile!

"""
Failure produced while interpreting or analyzing one physical source file.
"""
struct ItemFailure
    filepath::String
    measurement_id::String
    message::String
end

"""
Cancellation raised by package-owned background work.
"""
struct JobCancelled <: Exception end

const CANCEL_CALLBACK_KEY = :MeasurementBrowser_cancel_requested

"""
Return whether an exception represents cancellation of all contained work.
"""
function is_job_cancelled(error::Exception)::Bool
    error isa JobCancelled && return true
    error isa CompositeException || return false
    return !isempty(error.exceptions) && all(is_job_cancelled, error.exceptions)
end

"""
Run work with a task-local cancellation callback.
"""
function with_cancel(
    work::Function,
    cancel_requested::Union{Nothing,Function},
)::Any
    cancel_requested === nothing && return work()
    return task_local_storage(work, CANCEL_CALLBACK_KEY, cancel_requested)
end

"""
Stop the current package job when its task-local cancellation callback requests it.
"""
function check_cancel()::Nothing
    cancel_requested = get(task_local_storage(), CANCEL_CALLBACK_KEY, nothing)
    cancel_requested !== nothing && cancel_requested() && throw(JobCancelled())
    return nothing
end

collection_path_label(::Project, collection::AbstractVector{<:AbstractString})::String =
    join(collection, "_")

collection_path_key(collection::AbstractVector{<:AbstractString})::String = join(collection, "/")

"""
Parse a stored slash-separated collection key into the tuple used by the hierarchy index.
"""
function collection_path_tuple(key::AbstractString)::Tuple{Vararg{String}}
    stripped = strip(String(key))
    isempty(stripped) && error("Collection path key cannot be empty")
    segments = split(stripped, '/')
    any(isempty, segments) && error("Invalid collection path key '$key'")
    return Tuple(String.(segments))
end

"""
The internal metadata record for one logical item discovered inside a physical source file.

`collection` is the item's canonical placement in the tree (nested collection names);
`collection_metadata` is per-collection metadata merged from source-root metadata sources.
`parameters` describe acquisition settings known while interpreting the file; `stats` contain values
computed after the required context is available.
"""
struct ItemRecord
    unique_id::String
    filename::String
    filepath::String
    clean_title::String
    kind::Symbol
    timestamp::Union{DateTime,Nothing}
    collection::Vector{String}
    collection_metadata::Dict{Symbol,Any}
    parameters::Dict{Symbol,Any}
    stats::Dict{Symbol,Any}
end

"""
Construct an item record while normalizing its string fields.
"""
function ItemRecord(;
    filepath::AbstractString,
    kind::Symbol,
    collection::AbstractVector{<:AbstractString},
    clean_title::AbstractString,
    filename::AbstractString=basename(filepath),
    unique_id::AbstractString=filepath,
    timestamp::Union{DateTime,Nothing}=nothing,
    collection_metadata::Dict{Symbol,Any}=Dict{Symbol,Any}(),
    parameters::Dict{Symbol,Any}=Dict{Symbol,Any}(),
    stats::Dict{Symbol,Any}=Dict{Symbol,Any}(),
)::ItemRecord
    return ItemRecord(
        String(unique_id),
        String(filename),
        String(filepath),
        String(clean_title),
        kind,
        timestamp,
        String[String(segment) for segment in collection],
        collection_metadata,
        parameters,
        stats,
    )
end

"""
Copy an item record while replacing selected fields.
"""
function ItemRecord(
    record::ItemRecord;
    filepath::AbstractString=record.filepath,
    filename::AbstractString=record.filename,
    unique_id::AbstractString=record.unique_id,
    clean_title::AbstractString=record.clean_title,
    kind::Symbol=record.kind,
    timestamp::Union{DateTime,Nothing}=record.timestamp,
    collection::AbstractVector{<:AbstractString}=copy(record.collection),
    collection_metadata::Dict{Symbol,Any}=deepcopy(record.collection_metadata),
    parameters::Dict{Symbol,Any}=deepcopy(record.parameters),
    stats::Dict{Symbol,Any}=deepcopy(record.stats),
)::ItemRecord
    return ItemRecord(;
        filepath,
        filename,
        unique_id,
        clean_title,
        kind,
        timestamp,
        collection,
        collection_metadata,
        parameters,
        stats,
    )
end

"""
The normal concrete item the package ships and `register_item!` produces: a loaded, data-bearing
value handed to plot/view callbacks. It answers the `AbstractDataItem` contract from its own fields
and carries its payload as `item.data`. A project that needs more subtypes `AbstractDataItem`
directly instead, side by side with this type.

The engine materializes a `DataItem` from an internal `ItemRecord` plus a data payload at view time,
sharing the record's `collection`/`parameters`/`stats` so engine-computed metadata is already
present. `ItemRecord` is never a field of a `DataItem`.
"""
struct DataItem <: AbstractDataItem
    unique_id::String
    label::String
    kind::Symbol
    collection::Vector{String}
    parameters::Dict{Symbol,Any}
    stats::Dict{Symbol,Any}
    data::Any
end

"""Materialize a loaded item from an internal record and its (processed) payload."""
DataItem(record::ItemRecord, data)::DataItem = DataItem(
    record.unique_id,
    record.clean_title,
    record.kind,
    record.collection,
    record.parameters,
    record.stats,
    data,
)

"""
Construct a `DataItem` from an `entries` callback — the recipe API's per-item entry.

The engine fills `filepath`/`timestamp` and mints `unique_id` from the `SourceFile`, so a recipe
supplies only the metadata it knows: `kind`, `collection`, and optionally `label`/`parameters`. `data`
is omitted at scan (the engine resolves it later via `read`/`process`).
"""
function DataItem(;
    kind::Symbol,
    collection::AbstractVector{<:AbstractString},
    label::AbstractString="",
    parameters::Dict{Symbol,Any}=Dict{Symbol,Any}(),
    stats::Dict{Symbol,Any}=Dict{Symbol,Any}(),
    data=nothing,
    unique_id::AbstractString="",
)::DataItem
    return DataItem(
        String(unique_id),
        String(label),
        kind,
        String[String(segment) for segment in collection],
        parameters,
        stats,
        data,
    )
end

"""
Derive the internal `ItemRecord` from any item via the `AbstractDataItem` contract.

The contract carries no file identity, so the engine supplies `filepath`/`filename`/`timestamp` and
the minted `unique_id` from the scan-time `SourceFile`; the rest comes from the contract.
`collection_metadata` starts empty and is filled later from source-root metadata.
"""
function ItemRecord(
    item::AbstractDataItem;
    filepath::AbstractString,
    filename::AbstractString=basename(filepath),
    timestamp::Union{DateTime,Nothing}=nothing,
    kind::Symbol=Projects.kind(item),
    unique_id::AbstractString=Projects.item_id(item),
)::ItemRecord
    label = Projects.item_label(item)
    title = isempty(label) ?
        strip(join(filter(!isnothing, Any[timestamp, string(kind)]), " ")) :
        String(label)
    return ItemRecord(
        String(unique_id),
        String(filename),
        String(filepath),
        title,
        kind,
        timestamp,
        String[String(segment) for segment in Projects.collection(item)],
        Dict{Symbol,Any}(),
        Projects.parameters(item),
        Projects.stats(item),
    )
end

Projects.item_id(item::DataItem)::String = item.unique_id
Projects.item_label(item::DataItem)::String = item.label
Projects.kind(item::DataItem)::Symbol = item.kind
Projects.collection(item::DataItem)::Vector{String} = item.collection
Projects.parameters(item::DataItem)::Dict{Symbol,Any} = item.parameters
Projects.stats(item::DataItem)::Dict{Symbol,Any} = item.stats
Projects.item_data(item::DataItem) = item.data
Projects.read_data(item::DataItem) = item.data

"""
Extract a measurement timestamp from the supported filename conventions.
"""
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

include("ItemIndex/SourceFiles.jl")

"""
One node in the device hierarchy.
"""
struct HierarchyNode
    name::String
    kind::Symbol
    children::Vector{HierarchyNode}
    measurements::Vector{ItemRecord}
end

HierarchyNode(name::String, kind::Symbol) =
    HierarchyNode(name, kind, HierarchyNode[], ItemRecord[])

"""
The complete device tree and its indexes for one source root.
"""
struct Hierarchy
    root::HierarchyNode
    all_measurements::Vector{ItemRecord}
    root_path::String
    index::Dict{Tuple{Vararg{String}},HierarchyNode}
    has_collection_metadata::Bool
    project::Project
    skipped_count::Int
end

"""
The authoritative result of one completed filesystem scan.
"""
struct SourceScan
    root_path::String
    project::Project
    files::Vector{SourceFile}
    hierarchy::Hierarchy
    analysis_failures::Vector{ItemFailure}
end

"""Construct a successful scan with no recorded analysis failures."""
function SourceScan(
    root_path::String,
    project::Project,
    files::Vector{SourceFile},
    hierarchy::Hierarchy,
)::SourceScan
    return SourceScan(root_path, project, files, hierarchy, ItemFailure[])
end

"""
Convert a Roman-numeral path segment to its integer value.
"""
function roman_value(text::AbstractString)::Union{Nothing,Int}
    values = Dict('I' => 1, 'V' => 5, 'X' => 10, 'L' => 50)
    isempty(text) && return nothing
    total = 0
    previous = 0
    for character in reverse(uppercase(text))
        value = get(values, character, 0)
        value == 0 && return nothing
        total += value < previous ? -value : value
        previous = max(previous, value)
    end
    return total
end

"""
Return a tuple that sorts text segments and embedded integers naturally.
"""
function natural_key(text::AbstractString)::Tuple
    parts = Any[]
    for result in eachmatch(r"\d+|\D+", String(text))
        segment = result.match
        push!(parts, all(isdigit, segment) ?
            (1, parse(Int, segment)) :
            (0, lowercase(segment)))
    end
    return Tuple(parts)
end

item_timestamp_key(measurement::ItemRecord)::DateTime =
    measurement.timestamp === nothing ?
    DateTime(Dates.year(typemax(Date))) :
    measurement.timestamp

"""
Sort one hierarchy node recursively using natural names and measurement time.
"""
function Base.sort!(node::HierarchyNode)::HierarchyNode
    foreach(sort!, node.children)
    roman = !isempty(node.children) &&
        all(child -> roman_value(child.name) !== nothing, node.children)
    sort!(
        node.children;
        by=roman ?
            child -> roman_value(child.name) :
            child -> natural_key(child.name),
    )
    foreach(
        child -> sort!(child.measurements; by=item_timestamp_key),
        node.children,
    )
    return node
end

Base.sort!(hierarchy::Hierarchy)::Hierarchy = (
    sort!(hierarchy.root);
    hierarchy
)

"""
Create an empty hierarchy ready for progressive scan results.
"""
function Hierarchy(
    root_path::String,
    has_collection_metadata::Bool,
    project::Project,
    skipped_count::Int=0,
)::Hierarchy
    return Hierarchy(
        HierarchyNode("/", :root),
        ItemRecord[],
        root_path,
        Dict{Tuple{Vararg{String}},HierarchyNode}(),
        has_collection_metadata,
        project,
        skipped_count,
    )
end

"""
Insert one measurement into an existing hierarchy.
"""
function insert_item!(
    hierarchy::Hierarchy,
    measurement::ItemRecord,
)::ItemRecord
    parent = hierarchy.root
    for (depth, segment) in enumerate(measurement.collection)
        path = Tuple(measurement.collection[1:depth])
        child = get(hierarchy.index, path, nothing)
        if child === nothing
            kind = depth == length(measurement.collection) ? :leaf : :level
            child = HierarchyNode(segment, kind)
            push!(parent.children, child)
            hierarchy.index[path] = child
        end
        parent = child
    end
    push!(parent.measurements, measurement)
    push!(hierarchy.all_measurements, measurement)
    return measurement
end

"""
Build a complete device tree from a flat measurement list.
"""
function Hierarchy(
    measurements::Vector{ItemRecord},
    root_path::String,
    has_collection_metadata::Bool,
    project::Project,
    skipped_count::Int=0,
)::Hierarchy
    hierarchy = Hierarchy(
        root_path,
        has_collection_metadata,
        project,
        skipped_count,
    )
    foreach(measurement -> insert_item!(hierarchy, measurement), measurements)
    return sort!(hierarchy)
end

children(node::HierarchyNode)::Vector{HierarchyNode} = node.children

include("ItemIndex/Scanning.jl")

end
