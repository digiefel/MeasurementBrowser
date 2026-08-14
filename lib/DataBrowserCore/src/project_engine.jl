"""
Source-item interpretation: run `read` then `entries` for one source item and normalize the result.

This is the only place the source and the item pipeline meet. Identity minting and collection
placement happen here, once, for every project alike — the stages themselves never see the source
beyond `read`, and nothing in this file knows any particular project dialect exists.
"""

using DataBrowserAPI
using DataBrowserAPI: AbstractCollection, AbstractDataItem, AbstractDataSource, AbstractDataSourceItem, AbstractProject, annotate_collection_path, default_collection_path, source_id, source_item_path, source_item_timestamp
using DataBrowserSources: DirectorySource, index_source_file
import DataBrowserAPI:
    attach_record,
    collection,
    entries,
    id,
    label,
    metadata,
    read
import DataBrowserAPI.ItemIndex:
    CollectionIndex,
    CollectionInput,
    ItemFailure,
    ItemRecord,
    MetadataDict,
    collection_inputs,
    effective_record,
    resolve_collection_path!,
    metadata_dict

"""
Resolve where one interpreted item is placed, with the source in hand.

An item declaring no path of its own takes the source's default placement; every path is then
offered to the source, which may attach its own metadata to each level. This is the one point at
which placement sees the source, so the stages before it stay pure functions of values.
"""
function _placed_collection_path(
    source::AbstractDataSource,
    source_item::AbstractDataSourceItem,
    item::AbstractDataItem,
)::Vector{AbstractCollection}
    path = _collection_path(item)
    isempty(path) && (path = default_collection_path(source, source_item))
    isempty(path) && return path
    return AbstractCollection[segment for segment in annotate_collection_path(source, path)]
end

"""Validate one item's `collection(item)` contract result and return its concrete segments."""
function _collection_path(item::AbstractDataItem)::Vector{AbstractCollection}
    value = collection(item)
    value isa AbstractVector || throw(ArgumentError(
        "collection(::$(typeof(item))) must return a vector of AbstractCollection values; " *
        "got $(typeof(value))",
    ))
    all(segment -> segment isa AbstractCollection, value) || throw(ArgumentError(
        "collection(::$(typeof(item))) must return only AbstractCollection values; got $(repr(value))",
    ))
    return AbstractCollection[segment for segment in value]
end

"""
One completed source-item pass.

`records` are retained by the index. `interpreted_items` carry effective item data only on direct
interpretation paths; workspace workers put that data in the memory cache before publishing a
lightweight completion.
"""
struct SourceItemInterpretation
    records::Vector{ItemRecord}
    collection_paths::Vector{Vector{CollectionInput}}
    interpreted_items::Vector{AbstractDataItem}
    failures::Vector{ItemFailure}
    # The source item's display label, resolved once here so the UI and cache never rerun
    # project label code.
    source_item_label::String
end

"""
Interpret every logical data item produced by one source item.

The source is touched only by `read`; `entries` expands its result into items without going back to
the origin. Collection placement is applied here, where the source is legitimately in hand: an item
declaring no path of its own takes the source's default, and every path is offered to the source to
annotate. Items are otherwise passed through untouched, whatever their type. Processing and
analysis belong to the workspace work graph. `source_item_key` is the workspace-minted surrogate
stamped on every record; transient interpretations outside a workspace (such as `items_for_file`)
leave it at 0.
"""
function interpret_source_item(
    project::AbstractProject,
    source::AbstractDataSource,
    source_item::AbstractDataSourceItem;
    source_item_key::Int64=Int64(0),
)::SourceItemInterpretation
    source_item_id_value = id(source_item)
    source_item_path_value = source_item_path(source_item)
    source_item_label_value = label(source_item)
    loaded = read(project, source, source_item)
    handles = entries(project, source_item, loaded)
    handles isa AbstractVector || error(
        "entries(::$(typeof(project)), ::$(typeof(source_item)), ::$(typeof(loaded))) must " *
        "return a vector of items; got $(typeof(handles))",
    )
    item_count = length(handles)
    records = Vector{ItemRecord}(undef, item_count)
    collection_paths = Vector{Vector{CollectionInput}}(undef, item_count)
    interpreted_items = Vector{AbstractDataItem}(undef, item_count)
    source_metadata = metadata_dict(metadata(source_item))
    for (index, handle) in pairs(handles)
        record_metadata = merge(copy(source_metadata), metadata_dict(metadata(handle)))
        item_label_value = label(handle)
        record = ItemRecord(;
            source_item_key,
            source_item_path=source_item_path_value,
            source_item_timestamp=source_item_timestamp(source_item),
            id=id(handle),
            label=isempty(item_label_value) ? label(source_item) : item_label_value,
            type=typeof(handle),
            collection_key=nothing,
            metadata=record_metadata,
        )
        records[index] = record
        collection_paths[index] = collection_inputs(
            _placed_collection_path(source, source_item, handle))
        interpreted_items[index] = attach_record(handle, record)
    end
    return SourceItemInterpretation(
        records, collection_paths, interpreted_items, ItemFailure[],
        String(source_item_label_value))
end

"""
Interpret one physical file into the data items produced by project code.

Every item answers `collection(item)` with its concrete collection segments; registered items
carry the normalized segments adapted from their registration callback's string path.
"""
function items_for_file(
    project::AbstractProject,
    filepath::AbstractString;
    meta::Union{Nothing,Dict{Tuple{Vararg{String}},Dict{Symbol,Any}}}=nothing,
)::Vector{AbstractDataItem}
    source = DirectorySource(dirname(filepath); metadata_file=nothing)
    if meta !== nothing
        source.collection_metadata_entries = meta
        source.has_metadata = true
    end
    interpretation = interpret_source_item(
        project,
        source,
        index_source_file(filepath, dirname(filepath)),
    )
    collections = CollectionIndex(source_id(source))
    records = ItemRecord[
        ItemRecord(record; collection_key=resolve_collection_path!(collections, path))
        for (record, path) in zip(
            interpretation.records,
            interpretation.collection_paths,
        )
    ]
    return AbstractDataItem[
        attach_record(item, effective_record(collections, record))
        for (record, item) in zip(records, interpretation.interpreted_items)
    ]
end

