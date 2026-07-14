module ItemIndex

using Dates

import ..DataBrowserAPI:
    AbstractDataSource,
    AbstractDataSourceItem,
    AbstractDataItem,
    AbstractCollection,
    MetadataDict,
    MetadataValue,
    Project,
    cacheable,
    cacheable_data,
    collection,
    collection_record_id,
    collection_path_label,
    fingerprint,
    id,
    item_data,
    item_label,
    kind,
    label,
    metadata,
    source_id,
    source_item_id,
    source_item_label,
    source_item_path,
    source_item_timestamp,
    source_label

"""Failure produced while interpreting or analyzing one source item."""
struct ItemFailure
    source_item_id::String
    id::String
    message::String
end

collection_path_label(::Project, path::AbstractVector{<:AbstractCollection})::String =
    join(label.(path), "_")

"""Normalize one value into the supported metadata value union."""
function metadata_value(value)::MetadataValue
    value isa MetadataValue && return value
    value isa Integer && return Int64(value)
    value isa AbstractFloat && return Float64(value)
    value isa AbstractString && return String(value)
    value isa AbstractVector{Bool} && return Vector{Bool}(value)
    value isa AbstractVector{<:Integer} && return Int64.(value)
    value isa AbstractVector{<:AbstractFloat} && return Float64.(value)
    value isa AbstractVector{<:AbstractString} && return String.(value)
    throw(ArgumentError(
        "unsupported metadata value of type $(typeof(value)): $(repr(value)); " *
        "metadata must be scalars (Bool/Int64/Float64/String/Symbol/Date/DateTime/Missing) " *
        "or homogeneous Bool/Int64/Float64/String vectors",
    ))
end

"""Convert a symbol-keyed dictionary into a validated [`MetadataDict`](@ref)."""
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

"""Package-owned collection value created by the registration string-path adapter."""
struct RegisteredCollection <: AbstractCollection
    name::String
    metadata::MetadataDict
end

RegisteredCollection(name::AbstractString; metadata::AbstractDict=MetadataDict()) =
    RegisteredCollection(String(name), metadata_dict(metadata))

label(collection::RegisteredCollection)::String = collection.name
metadata(collection::RegisteredCollection)::MetadataDict = collection.metadata
id(collection::RegisteredCollection)::String = collection.name
Base.:(==)(left::RegisteredCollection, right::RegisteredCollection)::Bool =
    left.name == right.name
Base.isequal(left::RegisteredCollection, right::RegisteredCollection)::Bool =
    isequal(left.name, right.name)
Base.hash(collection::RegisteredCollection, seed::UInt)::UInt = hash(collection.name, seed)

"""
One transient normalized collection level produced during interpretation.

This contains only package-owned projections of a live `AbstractCollection`; it never retains the
user value itself.
"""
struct CollectionInput
    id::String
    label::String
    own_metadata::MetadataDict
    registration_name::Union{Nothing,String}
end

"""Resolve a live collection path into package-owned inputs without retaining user values."""
function collection_inputs(path::AbstractVector{<:AbstractCollection})::Vector{CollectionInput}
    inputs = CollectionInput[]
    parent_id = ""
    for value in path
        collection_id = collection_record_id(parent_id, value)
        push!(inputs, CollectionInput(
            collection_id,
            String(label(value)),
            metadata_dict(metadata(value)),
            value isa RegisteredCollection ? value.name : nothing,
        ))
        parent_id = collection_id
    end
    return inputs
end

"""Normalize a registration string path into package-owned collection inputs."""
collection_inputs(
    ::AbstractDataSource,
    names::AbstractVector{<:AbstractString},
)::Vector{CollectionInput} = collection_inputs(
    AbstractCollection[RegisteredCollection(name) for name in names],
)

"""
One package-owned indexed collection occurrence.

`id` is the final deterministic occurrence ID derived from the parent occurrence ID, the concrete
collection type, and the value returned by `id(::AbstractCollection)`. `key` is a compact
workspace/cache-local integer. `label` is display text. These are distinct contracts.
"""
struct CollectionRecord
    key::Int64
    id::String
    parent_key::Union{Nothing,Int64}
    label::String
    own_metadata::MetadataDict
    registration_name::Union{Nothing,String}
    analysis::MetadataDict
end

label(collection_record::CollectionRecord)::String = collection_record.label

"""
Flat collection records plus their parent, child, and item-membership indexes.

There is no root record and no tree object graph. `nothing` is the parent/membership key for items
at the source root.
"""
mutable struct CollectionIndex
    source_id::String
    records::Dict{Int64,CollectionRecord}
    key_by_id::Dict{String,Int64}
    children_by_parent::Dict{Union{Nothing,Int64},Vector{Int64}}
    item_ids_by_collection::Dict{Union{Nothing,Int64},Vector{String}}
    next_key::Int64
    skipped_count::Int
end

CollectionIndex(source_id::AbstractString, skipped_count::Integer=0) = CollectionIndex(
    String(source_id),
    Dict{Int64,CollectionRecord}(),
    Dict{String,Int64}(),
    Dict{Union{Nothing,Int64},Vector{Int64}}(),
    Dict{Union{Nothing,Int64},Vector{String}}(),
    Int64(1),
    Int(skipped_count),
)

"""Copy a collection index deeply enough for independent mutation."""
function Base.copy(index::CollectionIndex)::CollectionIndex
    return CollectionIndex(
        index.source_id,
        copy(index.records),
        copy(index.key_by_id),
        Dict(key => copy(values) for (key, values) in index.children_by_parent),
        Dict(key => copy(values) for (key, values) in index.item_ids_by_collection),
        index.next_key,
        index.skipped_count,
    )
end

"""Register one persisted collection record, validating keys, parents, and IDs."""
function register_collection!(index::CollectionIndex, collection_record::CollectionRecord)::Nothing
    collection_record.key > 0 || error("Collection record keys must be positive")
    collection_record.parent_key === nothing ||
        haskey(index.records, collection_record.parent_key) || error(
        "Collection '$(collection_record.label)' refers to missing parent key " *
        "$(collection_record.parent_key)",
    )
    existing_key = get(index.key_by_id, collection_record.id, nothing)
    if existing_key !== nothing
        existing_collection_record = index.records[existing_key]
        existing_collection_record.key == collection_record.key || error(
            "Collection id $(collection_record.id) is assigned to keys " *
            "$(existing_collection_record.key) and $(collection_record.key)",
        )
        return nothing
    end
    haskey(index.records, collection_record.key) && error(
        "Collection key $(collection_record.key) is assigned to multiple collection IDs",
    )
    index.records[collection_record.key] = collection_record
    index.key_by_id[collection_record.id] = collection_record.key
    children = get!(() -> Int64[], index.children_by_parent, collection_record.parent_key)
    collection_record.key in children || push!(children, collection_record.key)
    index.next_key = max(index.next_key, collection_record.key + 1)
    return nothing
end

"""Resolve one interpreted path into collection records and return its leaf key."""
function resolve_collection_path!(
    index::CollectionIndex,
    inputs::Vector{CollectionInput},
)::Union{Nothing,Int64}
    parent_key::Union{Nothing,Int64} = nothing
    labels = String[]
    for input in inputs
        push!(labels, input.label)
        key = get(index.key_by_id, input.id, nothing)
        if key === nothing
            key = index.next_key
            register_collection!(index, CollectionRecord(
                key,
                input.id,
                parent_key,
                input.label,
                copy(input.own_metadata),
                input.registration_name,
                MetadataDict(),
            ))
        else
            existing = index.records[key]
            existing.parent_key == parent_key || error(
                "Collection '$(join(labels, " / "))' resolved below inconsistent parents",
            )
            if existing.label != input.label ||
               existing.own_metadata != input.own_metadata ||
               existing.registration_name != input.registration_name
                index.records[key] = CollectionRecord(
                    existing.key,
                    existing.id,
                    existing.parent_key,
                    input.label,
                    copy(input.own_metadata),
                    input.registration_name,
                    existing.analysis,
                )
            end
        end
        parent_key = key
    end
    return parent_key
end

"""Return one record's ancestor-to-self collection keys."""
function collection_path_keys(index::CollectionIndex, key::Int64)::Vector{Int64}
    path = Int64[]
    current::Union{Nothing,Int64} = key
    while current !== nothing
        push!(path, current)
        current = index.records[current].parent_key
    end
    reverse!(path)
    return path
end

collection_path_keys(::CollectionIndex, ::Nothing)::Vector{Int64} = Int64[]

"""Return one collection's ancestor-to-self final occurrence IDs."""
collection_id_path(index::CollectionIndex, key::Int64)::Vector{String} =
    String[index.records[path_key].id for path_key in collection_path_keys(index, key)]
collection_id_path(::CollectionIndex, ::Nothing)::Vector{String} = String[]

"""Return display labels for one collection path."""
collection_location(index::CollectionIndex, key::Int64)::Vector{String} =
    String[index.records[path_key].label for path_key in collection_path_keys(index, key)]

"""Return stored registration names, or `nothing` when the indexed path is typed."""
function registration_names(
    index::CollectionIndex,
    key::Union{Nothing,Int64},
)::Union{Nothing,Vector{String}}
    key === nothing && return String[]
    names = String[]
    for path_key in collection_path_keys(index, key)
        collection_record = index.records[path_key]
        collection_record.registration_name === nothing && return nothing
        push!(names, collection_record.registration_name)
    end
    return names
end

"""Return child keys sorted by their resolved display labels."""
function sorted_child_keys(
    index::CollectionIndex,
    parent_key::Union{Nothing,Int64},
)::Vector{Int64}
    keys = copy(get(index.children_by_parent, parent_key, Int64[]))
    names = String[index.records[key].label for key in keys]
    roman = !isempty(names) && all(name -> roman_value(name) !== nothing, names)
    sort!(keys; by=roman ?
        key -> roman_value(index.records[key].label) :
        key -> natural_key(index.records[key].label))
    return keys
end

"""Attach one item id to its direct collection membership."""
function insert_item!(index::CollectionIndex, item_id::AbstractString, key)::Nothing
    ids = get!(() -> String[], index.item_ids_by_collection, key)
    stable_id = String(item_id)
    stable_id in ids || push!(ids, stable_id)
    return nothing
end

"""Remove one item membership and prune empty collection records toward the root."""
function remove_item!(index::CollectionIndex, item_id::AbstractString, key)::Nothing
    ids = get(index.item_ids_by_collection, key, nothing)
    ids === nothing || filter!(id -> id != item_id, ids)
    ids !== nothing && isempty(ids) && delete!(index.item_ids_by_collection, key)
    current = key
    while current !== nothing
        !isempty(get(index.item_ids_by_collection, current, String[])) && break
        !isempty(get(index.children_by_parent, current, Int64[])) && break
        collection_record = index.records[current]
        siblings = get(index.children_by_parent, collection_record.parent_key, nothing)
        siblings === nothing || filter!(child -> child != current, siblings)
        siblings !== nothing && isempty(siblings) &&
            delete!(index.children_by_parent, collection_record.parent_key)
        delete!(index.key_by_id, collection_record.id)
        delete!(index.records, current)
        current = collection_record.parent_key
    end
    return nothing
end

"""Return direct member ids for one collection key."""
collection_item_ids(index::CollectionIndex, key)::Vector{String} =
    get(index.item_ids_by_collection, key, String[])

"""Return whether one record has child collection records."""
has_collection_children(index::CollectionIndex, key::Int64)::Bool =
    !isempty(get(index.children_by_parent, key, Int64[]))

"""Return metadata inherited from collection parents through the direct collection."""
function collection_metadata(index::CollectionIndex, key)::MetadataDict
    effective = MetadataDict()
    key === nothing && return effective
    for path_key in collection_path_keys(index, key)
        merge!(effective, index.records[path_key].own_metadata)
    end
    return effective
end

"""Clear one collection record's analysis while preserving its resolved facts."""
function clear_collection_analysis!(index::CollectionIndex, key::Int64)::Nothing
    collection_record = index.records[key]
    index.records[key] = CollectionRecord(
        collection_record.key,
        collection_record.id,
        collection_record.parent_key,
        collection_record.label,
        collection_record.own_metadata,
        collection_record.registration_name,
        MetadataDict(),
    )
    return nothing
end

"""Replace one collection record's analysis result."""
function set_collection_analysis!(
    index::CollectionIndex,
    key::Int64,
    analysis::AbstractDict,
)::Nothing
    collection_record = index.records[key]
    index.records[key] = CollectionRecord(
        collection_record.key,
        collection_record.id,
        collection_record.parent_key,
        collection_record.label,
        collection_record.own_metadata,
        collection_record.registration_name,
        metadata_dict(analysis),
    )
    return nothing
end

"""The internal metadata record for one logical item discovered inside one source item."""
struct ItemRecord
    id::String
    source_item_id::String
    source_item_path::Union{Nothing,String}
    source_item_timestamp::Union{DateTime,Nothing}
    item_label::String
    kind::Symbol
    collection_key::Union{Nothing,Int64}
    metadata::MetadataDict
    item_fingerprint::Any
end

"""Construct an item record while normalizing its fields."""
function ItemRecord(;
    id::AbstractString,
    source_item_id::AbstractString,
    source_item_path::Union{Nothing,AbstractString}=nothing,
    source_item_timestamp::Union{DateTime,Nothing}=nothing,
    item_label::AbstractString,
    kind::Symbol,
    collection_key::Union{Nothing,Integer}=nothing,
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
        collection_key === nothing ? nothing : Int64(collection_key),
        metadata_dict(metadata),
        item_fingerprint,
    )
end

"""Copy an item record while replacing selected fields."""
function ItemRecord(
    record::ItemRecord;
    id::AbstractString=record.id,
    source_item_id::AbstractString=record.source_item_id,
    source_item_path::Union{Nothing,AbstractString}=record.source_item_path,
    source_item_timestamp::Union{DateTime,Nothing}=record.source_item_timestamp,
    item_label::AbstractString=record.item_label,
    kind::Symbol=record.kind,
    collection_key::Union{Nothing,Integer}=record.collection_key,
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
        collection_key,
        metadata,
        item_fingerprint,
    )
end

"""
Private carrier for ordinary data produced by `register_item!`.

Its collection path stays a `Vector{String}` at registration-facing boundaries. Interpretation
adapts those names separately when constructing package-owned collection records.
"""
struct RegisteredDataItem{D} <: AbstractDataItem
    id::String
    label::String
    registration::Symbol
    collection::Vector{String}
    data::D
    metadata::MetadataDict
end

"""Reconstruct registered data from a record, payload, and its registration string path."""
RegisteredDataItem(
    record::ItemRecord,
    data,
    path::Vector{String}=String[],
)::RegisteredDataItem = RegisteredDataItem(
    record.id,
    record.item_label,
    record.kind,
    path,
    data,
    record.metadata,
)

"""Copy registered data while replacing only its payload."""
RegisteredDataItem(item::RegisteredDataItem, data)::RegisteredDataItem = RegisteredDataItem(
    item.id,
    item.label,
    item.registration,
    item.collection,
    data,
    item.metadata,
)

"""Derive an unresolved item record from any item; collection membership is assigned on publish."""
function ItemRecord(
    item::AbstractDataItem;
    source_item::AbstractDataSourceItem,
    id::AbstractString=id(item),
    kind::Symbol=kind(item),
    metadata=metadata(item),
)::ItemRecord
    item_title = item_label(item)
    title = isempty(item_title) ? source_item_label(source_item) : String(item_title)
    return ItemRecord(;
        source_item_id=source_item_id(source_item),
        source_item_path=source_item_path(source_item),
        source_item_timestamp=source_item_timestamp(source_item),
        id,
        item_label=title,
        kind,
        collection_key=nothing,
        metadata,
        item_fingerprint=fingerprint(item),
    )
end

id(item::RegisteredDataItem)::String = item.id
item_label(item::RegisteredDataItem)::String = item.label
kind(item::RegisteredDataItem)::Symbol = item.registration
collection(item::RegisteredDataItem)::Vector{String} = item.collection
metadata(item::RegisteredDataItem)::MetadataDict = item.metadata
item_data(item::RegisteredDataItem) = item.data
fingerprint(item::RegisteredDataItem) = nothing
cacheable(item::RegisteredDataItem)::Bool = cacheable_data(item.data)

"""Return one item record's inherited collection metadata plus its own entries layer."""
function effective_metadata(index::CollectionIndex, record::ItemRecord)::MetadataDict
    effective = collection_metadata(index, record.collection_key)
    merge!(effective, record.metadata)
    return effective
end

"""Materialize a record with its inherited collection metadata."""
effective_record(index::CollectionIndex, record::ItemRecord)::ItemRecord =
    ItemRecord(record; metadata=effective_metadata(index, record))

"""The authoritative cached result of one completed source scan."""
struct SourceScan
    source_id::String
    source_label::String
    collections::CollectionIndex
    items::Vector{ItemRecord}
    analysis_failures::Vector{ItemFailure}
end

"""Construct a successful scan with no recorded analysis failures."""
SourceScan(
    source::AbstractDataSource,
    collections::CollectionIndex,
    items::Vector{ItemRecord},
) = SourceScan(source_id(source), source_label(source), collections, items, ItemFailure[])

"""Send one structured progress update when a callback is present."""
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
    on_progress((;
        phase,
        total_source_items,
        processed_source_items,
        loaded_items,
        skipped_source_items,
        current_source_item,
    ))
    return nothing
end

"""Convert a Roman-numeral path segment to its integer value."""
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

"""Return a tuple that sorts text segments and embedded integers naturally."""
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

end
