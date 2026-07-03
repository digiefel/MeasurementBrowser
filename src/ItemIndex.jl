module ItemIndex

using Dates
using DataFrames: AbstractDataFrame

import ..Profiling

import ..Projects
import ..Projects:
    AbstractDataSource,
    AbstractDataSourceItem,
    AbstractDataItem,
    Project,
    collection,
    collection_parameters,
    collection_path_label,
    data_items,
    fingerprint,
    finish_source_profile!,
    has_collection_parameters,
    id,
    item_data,
    kind,
    kind_label,
    parameters,
    project_name,
    reset_scan_profile!,
    source_id,
    source_item_id,
    source_item_label,
    source_item_path,
    source_item_timestamp,
    source_items,
    source_label

"""
Failure produced while interpreting or analyzing one source item.
"""
struct ItemFailure
    source_item_id::String
    id::String
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
The value types a `parameters`/`stats` entry may hold.

Kept deliberately narrow: these dicts describe items and collections with flat, queryable scalars
(exactly what `parse_metadata_value` emits) plus `Symbol` tags, `Missing` for absent values, and
homogeneous numeric/string vectors. The cache stores them as proper typed columns, so anything
outside this union is rejected at index time rather than silently blobbed.
"""
const MetadataValue = Union{
    Bool, Int64, Float64, String, Symbol, Date, DateTime, Missing,
    Vector{Bool}, Vector{Int64}, Vector{Float64}, Vector{String},
}

"""The typed dict the engine stores for item/collection `parameters` and `stats`."""
const MetadataDict = Dict{Symbol,MetadataValue}

"""
Coerce one value into a [`MetadataValue`](@ref), normalizing common near-misses (any `Integer` to
`Int64`, any `AbstractFloat` to `Float64`, `AbstractString` to `String`, and the corresponding
vectors). Throws `ArgumentError` for anything else.
"""
function metadata_value(value)::MetadataValue
    value isa MetadataValue && return value
    value isa Integer && return Int64(value)            # Bool already matched above
    value isa AbstractFloat && return Float64(value)
    value isa AbstractString && return String(value)
    value isa AbstractVector{Bool} && return Vector{Bool}(value)   # before the Integer branch
    value isa AbstractVector{<:Integer} && return Int64.(value)
    value isa AbstractVector{<:AbstractFloat} && return Float64.(value)
    value isa AbstractVector{<:AbstractString} && return String.(value)
    throw(ArgumentError(
        "unsupported metadata value of type $(typeof(value)): $(repr(value)); " *
        "parameters/stats must be scalars (Bool/Int64/Float64/String/Symbol/Date/DateTime/Missing) " *
        "or homogeneous Bool/Int64/Float64/String vectors",
    ))
end

"""
Convert any `Symbol`-keyed dict into a validated [`MetadataDict`](@ref). A `MetadataDict` passes
through unchanged; for anything else each value is checked via [`metadata_value`](@ref), with the
offending key named on failure.
"""
metadata_dict(dict::MetadataDict)::MetadataDict = dict
function metadata_dict(dict::AbstractDict)::MetadataDict
    out = MetadataDict()
    for (key, value) in dict
        sym = Symbol(key)
        out[sym] = try
            metadata_value(value)
        catch err
            err isa ArgumentError || rethrow()
            throw(ArgumentError("metadata key :$sym: $(err.msg)"))
        end
    end
    return out
end

"""
The internal metadata record for one logical item discovered inside one source item.
"""
struct ItemRecord
    id::String
    source_item_id::String
    source_item_fingerprint::Any
    source_item_path::Union{Nothing,String}
    source_item_timestamp::Union{DateTime,Nothing}
    item_label::String
    kind::Symbol
    collection::Vector{String}
    parameters::MetadataDict
    stats::MetadataDict
    item_fingerprint::Any
end

"""
Construct an item record while normalizing its string fields.
"""
function ItemRecord(;
    id::AbstractString,
    source_item_id::AbstractString,
    source_item_fingerprint=nothing,
    source_item_path::Union{Nothing,AbstractString}=nothing,
    source_item_timestamp::Union{DateTime,Nothing}=nothing,
    item_label::AbstractString,
    kind::Symbol,
    collection::AbstractVector{<:AbstractString},
    parameters::AbstractDict=MetadataDict(),
    stats::AbstractDict=MetadataDict(),
    item_fingerprint=nothing,
)::ItemRecord
    record_id = String(id)
    isempty(record_id) && error("ItemRecord id cannot be empty")
    return ItemRecord(
        record_id,
        String(source_item_id),
        source_item_fingerprint,
        source_item_path === nothing ? nothing : String(source_item_path),
        source_item_timestamp,
        String(item_label),
        kind,
        String[String(segment) for segment in collection],
        metadata_dict(parameters),
        metadata_dict(stats),
        item_fingerprint,
    )
end

"""
Copy an item record while replacing selected fields.
"""
function ItemRecord(
    record::ItemRecord;
    id::AbstractString=record.id,
    source_item_id::AbstractString=record.source_item_id,
    source_item_fingerprint=record.source_item_fingerprint,
    source_item_path::Union{Nothing,AbstractString}=record.source_item_path,
    source_item_timestamp::Union{DateTime,Nothing}=record.source_item_timestamp,
    item_label::AbstractString=record.item_label,
    kind::Symbol=record.kind,
    collection::AbstractVector{<:AbstractString}=copy(record.collection),
    parameters::AbstractDict=deepcopy(record.parameters),
    stats::AbstractDict=deepcopy(record.stats),
    item_fingerprint=record.item_fingerprint,
)::ItemRecord
    return ItemRecord(;
        id,
        source_item_id,
        source_item_fingerprint,
        source_item_path,
        source_item_timestamp,
        item_label,
        kind,
        collection,
        parameters,
        stats,
        item_fingerprint,
    )
end

"""
The normal concrete item the package ships and `register_item!` produces: a loaded, data-bearing
value handed to plot/view callbacks. It answers the `AbstractDataItem` contract from its own fields
and carries its loaded data as `item.data`. A project that needs more subtypes `AbstractDataItem`
directly instead, side by side with this type.

The engine materializes a `DataItem` from an internal `ItemRecord` plus loaded data at view time,
sharing the record's `collection`/`parameters`/`stats` so engine-computed metadata is already
present. `ItemRecord` is never a field of a `DataItem`.
"""
struct DataItem <: AbstractDataItem
    id::String
    label::String
    kind::Symbol
    collection::Vector{String}
    parameters::Dict{Symbol,Any}
    stats::Dict{Symbol,Any}
    data::Any
end

"""Materialize a loaded item from an internal record and its processed data."""
DataItem(record::ItemRecord, data)::DataItem = DataItem(
    record.id,
    record.item_label,
    record.kind,
    record.collection,
    record.parameters,
    record.stats,
    data,
)

"""Copy one `DataItem`, replacing only its data."""
DataItem(item::DataItem, data)::DataItem = DataItem(
    item.id,
    item.label,
    item.kind,
    item.collection,
    item.parameters,
    item.stats,
    data,
)

"""
Construct a `DataItem` from an `entries` callback — the recipe API's per-item entry.

A recipe supplies the metadata it knows: `kind`, `collection`, and optionally
`label`/`parameters`/`id`. `data` carries the raw per-item data; an optional `process` callback
can return another item before views receive it.
"""
function DataItem(;
    kind::Symbol,
    collection::AbstractVector{<:AbstractString},
    label::AbstractString="",
    parameters::Dict{Symbol,Any}=Dict{Symbol,Any}(),
    stats::Dict{Symbol,Any}=Dict{Symbol,Any}(),
    data=nothing,
    id::AbstractString="",
)::DataItem
    return DataItem(
        String(id),
        String(label),
        kind,
        String[String(segment) for segment in collection],
        parameters,
        stats,
        data,
    )
end

"""
Derive the internal `ItemRecord` from any item via the source and item contracts.
"""
function ItemRecord(
    item::AbstractDataItem;
    source_item::AbstractDataSourceItem,
    kind::Symbol=Projects.kind(item),
    parameters=Projects.parameters(item),
)::ItemRecord
    label = Projects.item_label(item)
    title = isempty(label) ?
        strip(join(filter(!isnothing, Any[source_item_timestamp(source_item), string(kind)]), " ")) :
        String(label)
    return ItemRecord(;
        source_item_id=source_item_id(source_item),
        source_item_fingerprint=fingerprint(source_item),
        source_item_path=source_item_path(source_item),
        source_item_timestamp=source_item_timestamp(source_item),
        id=Projects.id(item),
        item_label=title,
        kind,
        collection=Projects.collection(item),
        parameters,
        stats=Projects.stats(item),
        item_fingerprint=fingerprint(item),
    )
end

Projects.id(item::DataItem)::String = item.id
Projects.item_label(item::DataItem)::String = item.label
Projects.kind(item::DataItem)::Symbol = item.kind
Projects.collection(item::DataItem)::Vector{String} = item.collection
Projects.parameters(item::DataItem)::Dict{Symbol,Any} = item.parameters
Projects.stats(item::DataItem)::Dict{Symbol,Any} = item.stats
Projects.item_data(item::DataItem) = item.data
Projects.fingerprint(item::DataItem) = nothing

# The built-in item opts DataFrame values and views into the native columnar cache. Other data types
# remain source-backed until they receive their own native storage method. Type-API items opt in
# themselves.
Projects.cacheable(item::DataItem)::Bool = item.data isa AbstractDataFrame

"""
One node in the collection hierarchy.
"""
struct HierarchyNode
    name::String
    kind::Symbol
    parameters::MetadataDict
    children::Vector{HierarchyNode}
    items::Vector{ItemRecord}
    stats::MetadataDict
end

HierarchyNode(name::String, kind::Symbol) =
    HierarchyNode(
        name,
        kind,
        MetadataDict(),
        HierarchyNode[],
        ItemRecord[],
        MetadataDict(),
    )

"""
The complete collection tree and its indexes for one source.
"""
struct Hierarchy
    root::HierarchyNode
    all_items::Vector{ItemRecord}
    source_id::String
    index::Dict{Tuple{Vararg{String}},HierarchyNode}
    has_collection_parameters::Bool
    skipped_count::Int
end

"""
The authoritative result of one completed source scan.
"""
struct SourceScan
    source_id::String
    source_label::String
    source_item_fingerprints::Dict{String,Any}
    hierarchy::Hierarchy
    analysis_failures::Vector{ItemFailure}
end

"""Construct a successful scan with no recorded analysis failures."""
function SourceScan(
    source::AbstractDataSource,
    source_item_fingerprints::Dict{String,Any},
    hierarchy::Hierarchy,
)::SourceScan
    return SourceScan(
        source_id(source),
        source_label(source),
        source_item_fingerprints,
        hierarchy,
        ItemFailure[],
    )
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

item_timestamp_key(item::ItemRecord)::DateTime =
    item.source_item_timestamp === nothing ?
    DateTime(Dates.year(typemax(Date))) :
    item.source_item_timestamp

"""
Sort one hierarchy node recursively using natural names and source-item time.
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
        child -> sort!(child.items; by=item_timestamp_key),
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
    source_id::String,
    has_collection_parameters::Bool,
    skipped_count::Int=0,
)::Hierarchy
    return Hierarchy(
        HierarchyNode("/", :root),
        ItemRecord[],
        source_id,
        Dict{Tuple{Vararg{String}},HierarchyNode}(),
        has_collection_parameters,
        skipped_count,
    )
end

"""
Insert one item into an existing hierarchy.
"""
function insert_item!(
    hierarchy::Hierarchy,
    item::ItemRecord,
)::ItemRecord
    parent = hierarchy.root
    for (depth, segment) in enumerate(item.collection)
        path = Tuple(item.collection[1:depth])
        child = get(hierarchy.index, path, nothing)
        if child === nothing
            kind = depth == length(item.collection) ? :leaf : :level
            child = HierarchyNode(segment, kind)
            push!(parent.children, child)
            hierarchy.index[path] = child
        end
        parent = child
    end
    push!(parent.items, item)
    push!(hierarchy.all_items, item)
    return item
end

function _effective_parameters(
    source::AbstractDataSource,
    collection_path::AbstractVector{<:AbstractString},
    local_parameters::AbstractDict,
)::MetadataDict
    parameters = metadata_dict(collection_parameters(source, collection_path))
    merge!(parameters, metadata_dict(local_parameters))
    return parameters
end

"""Apply source collection parameters to hierarchy nodes."""
function apply_collection_parameters!(
    hierarchy::Hierarchy,
    source::AbstractDataSource,
)::Nothing
    empty!(hierarchy.root.parameters)

    function visit!(
        node::HierarchyNode,
        path::Vector{String},
        inherited::MetadataDict,
    )::Nothing
        effective = copy(inherited)
        isempty(path) || merge!(effective, metadata_dict(collection_parameters(source, path)))
        empty!(node.parameters)
        merge!(node.parameters, effective)
        for child in node.children
            visit!(child, [path; child.name], effective)
        end
        return nothing
    end

    for child in hierarchy.root.children
        visit!(child, String[child.name], hierarchy.root.parameters)
    end
    return nothing
end

"""Return one record's source-inherited and item-local parameters."""
function effective_parameters(hierarchy::Hierarchy, record::ItemRecord)::MetadataDict
    node = get(hierarchy.index, Tuple(record.collection), nothing)
    effective = node === nothing ? MetadataDict() : copy(node.parameters)
    merge!(effective, record.parameters)
    return effective
end

"""Materialize a record with the effective parameters consumed by project callbacks."""
function effective_record(hierarchy::Hierarchy, record::ItemRecord)::ItemRecord
    return ItemRecord(record; parameters=effective_parameters(hierarchy, record))
end

"""
Build a complete collection tree from a flat item list.
"""
function Hierarchy(
    items::Vector{ItemRecord},
    source::AbstractDataSource,
    skipped_count::Int=0,
)::Hierarchy
    hierarchy = Hierarchy(
        source_id(source),
        has_collection_parameters(source),
        skipped_count,
    )
    foreach(item -> insert_item!(hierarchy, item), items)
    apply_collection_parameters!(hierarchy, source)
    return sort!(hierarchy)
end

children(node::HierarchyNode)::Vector{HierarchyNode} = node.children

include("ItemIndex/Scanning.jl")

end
