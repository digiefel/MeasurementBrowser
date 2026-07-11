module ItemIndex

using Dates

import ..DataBrowserAPI
import ..DataBrowserAPI:
    AbstractDataSource,
    AbstractDataSourceItem,
    AbstractDataItem,
    MetadataDict,
    MetadataValue,
    Project,
    collection,
    collection_metadata,
    collection_path_label,
    fingerprint,
    has_collection_metadata,
    id,
    item_data,
    item_label,
    kind,
    metadata,
    source_id,
    source_item_id,
    source_item_path,
    source_item_timestamp,
    source_label,
    cacheable,
    cacheable_data

"""
Failure produced while interpreting or analyzing one source item.
"""
struct ItemFailure
    source_item_id::String
    id::String
    message::String
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
        "metadata must be scalars (Bool/Int64/Float64/String/Symbol/Date/DateTime/Missing) " *
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
    source_item_path::Union{Nothing,String}
    source_item_timestamp::Union{DateTime,Nothing}
    item_label::String
    kind::Symbol
    collection::Vector{String}
    metadata::MetadataDict
    item_fingerprint::Any
end

"""
Construct an item record while normalizing its string fields.
"""
function ItemRecord(;
    id::AbstractString,
    source_item_id::AbstractString,
    source_item_path::Union{Nothing,AbstractString}=nothing,
    source_item_timestamp::Union{DateTime,Nothing}=nothing,
    item_label::AbstractString,
    kind::Symbol,
    collection::AbstractVector{<:AbstractString},
    metadata::AbstractDict=MetadataDict(),
    item_fingerprint=nothing,
)::ItemRecord
    record_id = String(id)
    isempty(record_id) && error("ItemRecord id cannot be empty")
    return ItemRecord(
        record_id,
        String(source_item_id),
        source_item_path === nothing ? nothing : String(source_item_path),
        source_item_timestamp,
        String(item_label),
        kind,
        String[String(segment) for segment in collection],
        metadata_dict(metadata),
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
    source_item_path::Union{Nothing,AbstractString}=record.source_item_path,
    source_item_timestamp::Union{DateTime,Nothing}=record.source_item_timestamp,
    item_label::AbstractString=record.item_label,
    kind::Symbol=record.kind,
    collection::AbstractVector{<:AbstractString}=copy(record.collection),
    metadata::AbstractDict=deepcopy(record.metadata),
    item_fingerprint=record.item_fingerprint,
)::ItemRecord
    return ItemRecord(;
        id,
        source_item_id,
        source_item_path,
        source_item_timestamp,
        item_label,
        kind,
        collection,
        metadata,
        item_fingerprint,
    )
end

"""
The normal concrete item the package ships and `register_item!` produces: a loaded, data-bearing
value handed to plot/view callbacks. It answers the `AbstractDataItem` contract from its own fields
and carries its loaded data as `item.data`. A project that needs more subtypes `AbstractDataItem`
directly instead, side by side with this type.

The engine materializes a `DataItem` from an internal `ItemRecord` plus loaded data at view time,
sharing the record's `collection`/`metadata` so engine-computed metadata is already present.
`ItemRecord` is never a field of a `DataItem`.
"""
struct DataItem <: AbstractDataItem
    id::String
    label::String
    kind::Symbol
    collection::Vector{String}
    metadata::Dict{Symbol,Any}
    data::Any
end

"""Materialize a loaded item from an internal record and its processed data."""
DataItem(record::ItemRecord, data)::DataItem = DataItem(
    record.id,
    record.item_label,
    record.kind,
    record.collection,
    record.metadata,
    data,
)

"""Copy one `DataItem`, replacing only its data."""
DataItem(item::DataItem, data)::DataItem = DataItem(
    item.id,
    item.label,
    item.kind,
    item.collection,
    item.metadata,
    data,
)

"""
Construct a `DataItem` from an `entries` callback — the recipe API's per-item entry.

A recipe supplies the metadata it knows: `kind`, `collection`, and optionally
`label`/`metadata`/`id`. `data` carries the raw per-item data; an optional `process` callback
can return another item before views receive it.
"""
function DataItem(;
    kind::Symbol,
    collection::AbstractVector{<:AbstractString},
    label::AbstractString="",
    metadata::Dict{Symbol,Any}=Dict{Symbol,Any}(),
    data=nothing,
    id::AbstractString="",
)::DataItem
    return DataItem(
        String(id),
        String(label),
        kind,
        String[String(segment) for segment in collection],
        metadata,
        data,
    )
end

"""
Derive the internal `ItemRecord` from any item via the source and item contracts.
"""
function ItemRecord(
    item::AbstractDataItem;
    source_item::AbstractDataSourceItem,
    kind::Symbol=kind(item),
    metadata=metadata(item),
)::ItemRecord
    label = item_label(item)
    title = isempty(label) ?
        strip(join(filter(!isnothing, Any[source_item_timestamp(source_item), string(kind)]), " ")) :
        String(label)
    return ItemRecord(;
        source_item_id=source_item_id(source_item),
        source_item_path=source_item_path(source_item),
        source_item_timestamp=source_item_timestamp(source_item),
        id=id(item),
        item_label=title,
        kind,
        collection=collection(item),
        metadata,
        item_fingerprint=fingerprint(item),
    )
end

id(item::DataItem)::String = item.id
item_label(item::DataItem)::String = item.label
kind(item::DataItem)::Symbol = item.kind
collection(item::DataItem)::Vector{String} = item.collection
metadata(item::DataItem)::Dict{Symbol,Any} = item.metadata
item_data(item::DataItem) = item.data
fingerprint(item::DataItem) = nothing

# The built-in item delegates to the payload trait: anything tabular is cacheable by default.
# Type-API items opt in themselves.
cacheable(item::DataItem)::Bool = cacheable_data(item.data)

"""
One node in the collection hierarchy.
"""
struct HierarchyNode
    name::String
    kind::Symbol
    metadata::MetadataDict
    children::Vector{HierarchyNode}
    items::Vector{ItemRecord}
    analysis::MetadataDict
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

const _AllItemsCache = Base.RefValue{Union{Nothing,Vector{ItemRecord}}}
_lazy_all_items()::_AllItemsCache = _AllItemsCache(nothing)

"""
The complete collection tree and its indexes for one source.

`all_items` is a flat view over every record in the tree. It is materialized lazily and cached
(`all_items_cache`) rather than rebuilt on every edit: the tree nodes own the records, and a
progressive scan publishes thousands of single-item edits, so eagerly rebuilding an O(total)
vector per publish is quadratic and floods the GC. Readers go through [`all_items`](@ref); the
cache is left empty by edits and filled on first read.
"""
struct Hierarchy
    root::HierarchyNode
    all_items_cache::_AllItemsCache
    source_id::String
    index::Dict{Tuple{Vararg{String}},HierarchyNode}
    has_collection_metadata::Bool
    skipped_count::Int
end

"""Push every record under one node into `acc`, in tree order."""
function _collect_items!(acc::Vector{ItemRecord}, node::HierarchyNode)::Nothing
    append!(acc, node.items)
    for child in node.children
        _collect_items!(acc, child)
    end
    return nothing
end

"""
Every record in the hierarchy as one flat vector, materialized on first read and cached.

Published hierarchies never mutate their node `items` in place (edits clone touched nodes), so the
cache stays valid for the life of the hierarchy object.
"""
function all_items(hierarchy::Hierarchy)::Vector{ItemRecord}
    cached = hierarchy.all_items_cache[]
    cached === nothing || return cached
    items = ItemRecord[]
    _collect_items!(items, hierarchy.root)
    hierarchy.all_items_cache[] = items
    return items
end

"""
The authoritative result of one completed source scan.
"""
struct SourceScan
    source_id::String
    source_label::String
    hierarchy::Hierarchy
    analysis_failures::Vector{ItemFailure}
end

"""Construct a successful scan with no recorded analysis failures."""
function SourceScan(
    source::AbstractDataSource,
    hierarchy::Hierarchy,
)::SourceScan
    return SourceScan(
        source_id(source),
        source_label(source),
        hierarchy,
        ItemFailure[],
    )
end

"""
Send one structured progress update when a callback is present.
"""
function emit_progress(
    on_progress::Union{Nothing,Function};
    phase::Symbol,
    total_source_items::Int,
    processed_source_items::Int,
    loaded_items::Int,
    skipped_source_items::Int,
    current_source_item::String="",
)::Nothing
    on_progress === nothing && return nothing
    on_progress((
        phase=phase,
        total_source_items=total_source_items,
        processed_source_items=processed_source_items,
        loaded_items=loaded_items,
        skipped_source_items=skipped_source_items,
        current_source_item=current_source_item,
    ))
    return nothing
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

"""Sort one node's children using Roman-numeral or natural name order."""
function sort_children!(node::HierarchyNode)::HierarchyNode
    roman = !isempty(node.children) &&
        all(child -> roman_value(child.name) !== nothing, node.children)
    sort!(
        node.children;
        by=roman ?
            child -> roman_value(child.name) :
            child -> natural_key(child.name),
    )
    return node
end

"""Sort one node's own items by source-item time."""
sort_items!(node::HierarchyNode)::HierarchyNode =
    (sort!(node.items; by=item_timestamp_key); node)

"""
Sort one hierarchy node recursively using natural names and source-item time.
"""
function Base.sort!(node::HierarchyNode)::HierarchyNode
    foreach(sort!, node.children)
    sort_children!(node)
    foreach(sort_items!, node.children)
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
    has_collection_metadata::Bool,
    skipped_count::Int=0,
)::Hierarchy
    return Hierarchy(
        HierarchyNode("/", :root),
        _lazy_all_items(),
        source_id,
        Dict{Tuple{Vararg{String}},HierarchyNode}(),
        has_collection_metadata,
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
    hierarchy.all_items_cache[] = nothing
    return item
end

function _effective_metadata(
    source::AbstractDataSource,
    collection_path::AbstractVector{<:AbstractString},
    local_metadata::AbstractDict,
)::MetadataDict
    metadata = metadata_dict(collection_metadata(source, collection_path))
    merge!(metadata, metadata_dict(local_metadata))
    return metadata
end

"""
Set one node's metadata from its parent's already-effective metadata plus the source's own
collection metadata for the path. Parents must be refreshed before children (depth order).
"""
function refresh_node_metadata!(
    hierarchy::Hierarchy,
    source::AbstractDataSource,
    path::Tuple{Vararg{String}},
)::Nothing
    node = get(hierarchy.index, path, nothing)
    node === nothing && return nothing
    inherited = length(path) == 1 ?
        hierarchy.root.metadata :
        hierarchy.index[path[1:end-1]].metadata
    effective = copy(inherited)
    merge!(effective, metadata_dict(collection_metadata(source, collect(path))))
    empty!(node.metadata)
    merge!(node.metadata, effective)
    return nothing
end

"""Apply source collection metadata to hierarchy nodes."""
function apply_collection_metadata!(
    hierarchy::Hierarchy,
    source::AbstractDataSource,
)::Nothing
    empty!(hierarchy.root.metadata)
    for path in sort!(collect(keys(hierarchy.index)); by=length)
        refresh_node_metadata!(hierarchy, source, path)
    end
    return nothing
end

"""Return one record's source-inherited and item-local metadata."""
function effective_metadata(hierarchy::Hierarchy, record::ItemRecord)::MetadataDict
    node = get(hierarchy.index, Tuple(record.collection), nothing)
    effective = node === nothing ? MetadataDict() : copy(node.metadata)
    merge!(effective, record.metadata)
    return effective
end

"""Materialize a record with the effective metadata consumed by project callbacks."""
function effective_record(hierarchy::Hierarchy, record::ItemRecord)::ItemRecord
    return ItemRecord(record; metadata=effective_metadata(hierarchy, record))
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
        has_collection_metadata(source),
        skipped_count,
    )
    foreach(item -> insert_item!(hierarchy, item), items)
    apply_collection_metadata!(hierarchy, source)
    return sort!(hierarchy)
end

children(node::HierarchyNode)::Vector{HierarchyNode} = node.children

# ------------------------------- Copy-on-write hierarchy edits -----------------------------------
#
# Publishing swaps `workspace.index.hierarchy` as one reference assignment so concurrent readers
# always traverse a complete tree. A `HierarchyEdit` keeps that contract at incremental cost: it
# clones only the nodes it mutates (plus their ancestor chain), shares every untouched subtree with
# the original, and re-sorts/re-parameterizes only what changed, using the same per-node primitives
# as the full build.

const CollectionPath = Tuple{Vararg{String}}

"""One in-progress copy-on-write edit of a hierarchy."""
mutable struct HierarchyEdit
    hierarchy::Hierarchy
    cloned::Set{CollectionPath}
    items_dirty::Set{CollectionPath}
    children_dirty::Set{CollectionPath}
    created::Vector{CollectionPath}
    pruned::Bool
end

_clone_node(node::HierarchyNode)::HierarchyNode = HierarchyNode(
    node.name,
    node.kind,
    copy(node.metadata),
    copy(node.children),
    copy(node.items),
    copy(node.analysis),
)

"""
Begin a copy-on-write edit; the original hierarchy and its nodes are never mutated.

The edit clones only touched nodes and their ancestors; `all_items` is left lazy, so neither the
original nor the result pays to rebuild the flat item view — it is derived from the tree on demand.
"""
function edit_hierarchy(hierarchy::Hierarchy)::HierarchyEdit
    return HierarchyEdit(
        Hierarchy(
            _clone_node(hierarchy.root),
            _lazy_all_items(),
            hierarchy.source_id,
            copy(hierarchy.index),
            hierarchy.has_collection_metadata,
            hierarchy.skipped_count,
        ),
        Set{CollectionPath}([()]),
        Set{CollectionPath}(),
        Set{CollectionPath}(),
        CollectionPath[],
        false,
    )
end

"""Whether an edit created or pruned collection nodes (as opposed to only moving items)."""
edit_changed_structure(edit::HierarchyEdit)::Bool = !isempty(edit.created) || edit.pruned

_node_at(hierarchy::Hierarchy, path::CollectionPath)::Union{Nothing,HierarchyNode} =
    path === () ? hierarchy.root : get(hierarchy.index, path, nothing)

"""Return the node at `path` cloned for mutation, cloning uncloned ancestors as needed."""
function editable_node!(edit::HierarchyEdit, path::CollectionPath)::HierarchyNode
    hierarchy = edit.hierarchy
    node = hierarchy.root
    for depth in eachindex(path)
        prefix = path[1:depth]
        child = hierarchy.index[prefix]
        if !(prefix in edit.cloned)
            child = _clone_node(child)
            position = findfirst(existing -> existing.name == child.name, node.children)
            node.children[position] = child
            hierarchy.index[prefix] = child
            push!(edit.cloned, prefix)
        else
            child = hierarchy.index[prefix]
        end
        node = child
    end
    return node
end

"""Insert one record, creating (and metadata-marking) any missing collection nodes."""
function insert_record!(edit::HierarchyEdit, record::ItemRecord)::Nothing
    hierarchy = edit.hierarchy
    leaf_path = Tuple(record.collection)
    for depth in eachindex(leaf_path)
        prefix = leaf_path[1:depth]
        haskey(hierarchy.index, prefix) && continue
        parent = editable_node!(edit, prefix[1:end-1])
        kind = depth == length(leaf_path) ? :leaf : :level
        child = HierarchyNode(String(leaf_path[depth]), kind)
        push!(parent.children, child)
        hierarchy.index[prefix] = child
        push!(edit.cloned, prefix)
        push!(edit.created, prefix)
        push!(edit.children_dirty, prefix[1:end-1])
    end
    node = editable_node!(edit, leaf_path)
    push!(node.items, record)
    push!(edit.items_dirty, leaf_path)
    return nothing
end

"""Remove records by id, pruning collection nodes the removal leaves empty."""
function remove_records!(edit::HierarchyEdit, records::Vector{ItemRecord})::Nothing
    isempty(records) && return nothing
    hierarchy = edit.hierarchy
    ids = Set(record.id for record in records)
    for path in unique(Tuple(record.collection) for record in records)
        node = editable_node!(edit, path)
        filter!(item -> !(item.id in ids), node.items)
        while path !== () && isempty(node.items) && isempty(node.children)
            parent_path = path[1:end-1]
            parent = editable_node!(edit, parent_path)
            position = findfirst(child -> child.name == node.name, parent.children)
            deleteat!(parent.children, position)
            delete!(hierarchy.index, path)
            delete!(edit.items_dirty, path)
            delete!(edit.children_dirty, path)
            edit.pruned = true
            path, node = parent_path, parent
        end
    end
    return nothing
end

"""Clear the analysis of one collection node if it still exists."""
function clear_node_analysis!(edit::HierarchyEdit, path::CollectionPath)::Nothing
    _node_at(edit.hierarchy, path) === nothing && return nothing
    empty!(editable_node!(edit, path).analysis)
    return nothing
end

"""
Finish an edit: metadata-mark created nodes, re-sort what changed, and return the new hierarchy,
ready to swap in with one reference assignment.
"""
function finish_edit!(edit::HierarchyEdit, source::AbstractDataSource)::Hierarchy
    hierarchy = edit.hierarchy
    for path in sort!(edit.created; by=length)
        refresh_node_metadata!(hierarchy, source, path)
    end
    for path in edit.items_dirty
        node = _node_at(hierarchy, path)
        node === nothing || sort_items!(node)
    end
    for path in edit.children_dirty
        node = _node_at(hierarchy, path)
        node === nothing || sort_children!(node)
    end
    # The tree now reflects every insert and prune; the flat item view is derived from it lazily,
    # so an incremental publish no longer rebuilds an O(total) vector.
    return hierarchy
end

end
