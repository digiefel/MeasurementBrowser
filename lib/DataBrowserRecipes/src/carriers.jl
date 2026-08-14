"""
Package-owned carriers for the registration dialect.

A registration name is a type parameter, not a field, so every stage after `read` selects its
recipe by dispatch and the engine can specialize the whole call behind one function barrier.
`RegisteredDataItem{:pund,DataFrame}` is a distinct concrete type, which is what lets registered
items behave like any other typed item rather than like a tagged bag.
"""

"""
Private carrier handing one registration's `read` output to its `entries` stage.

`detect` runs inside the dialect's `read`, which is why the `{K}` tag first appears here.
"""
struct RegisteredReadResult{K,D}
    data::D
    metadata::MetadataDict
end

RegisteredReadResult{K}(data::D, metadata::MetadataDict) where {K,D} =
    RegisteredReadResult{K,D}(data, metadata)

"""Result of a `read` whose source item matched no registration; `entries` yields no items."""
struct NoMatch end

"""
Private carrier for ordinary data produced by `register_item!`.

Adaptation converts the registration callback's collection strings into normalized
`AbstractCollection` segments, so the carrier answers the generic item contract directly. Before
interpretation normalizes the item, `id` holds only the callback-supplied sibling key (or `""`);
the carrier delivered by interpretation is rebuilt on its record and carries the final minted id.
"""
struct RegisteredDataItem{K,D} <: AbstractDataItem
    id::String
    label::String
    collection::Vector{AbstractCollection}
    data::D
    metadata::MetadataDict
end

RegisteredDataItem{K}(
    id::AbstractString,
    label::AbstractString,
    collection::Vector{AbstractCollection},
    data::D,
    metadata::MetadataDict,
) where {K,D} = RegisteredDataItem{K,D}(String(id), String(label), collection, data, metadata)

"""Copy registered data while replacing only its payload."""
RegisteredDataItem(item::RegisteredDataItem{K}, data) where {K} = RegisteredDataItem{K}(
    item.id, item.label, item.collection, data, item.metadata)

id(item::RegisteredDataItem)::String = item.id
label(item::RegisteredDataItem)::String = item.label
label(::Type{<:RegisteredDataItem{K}}) where {K} = K

collection(item::RegisteredDataItem)::Vector{AbstractCollection} = item.collection
metadata(item::RegisteredDataItem)::MetadataDict = item.metadata
item_data(item::RegisteredDataItem) = item.data

"""
A registered carrier adopts its normalized record wholesale, and the index's collection path when
one is supplied. Its own segments survive interpretation, where the index has nothing yet.

The record's registration is validated against the carrier's tag: a mismatch means a cached record
was paired with the wrong payload, which must fail loudly rather than silently produce an item of
the wrong registration.
"""
function attach_record(
    item::RegisteredDataItem{K},
    record::ItemRecord,
    path::AbstractVector=AbstractCollection[],
) where {K}
    label(record.type) === K || error(
        "Cached record '$(record.id)' is registered as :$(label(record.type)) but its payload was " *
        "rebuilt as :$K",
    )
    return RegisteredDataItem{K}(
        record.id,
        record.label,
        isempty(path) ? item.collection : AbstractCollection[s for s in path],
        item.data,
        record.metadata,
    )
end

"""
Rebuild a registered carrier from a cached payload.

The stored id is handed back verbatim, and `attach_record` restores the label and collection path
from the record and index straight afterwards. The registration comes from the type, so it survives
the round trip without ever entering the payload.
"""
# `Type{<:RegisteredDataItem{K}}` accepts both the seeded UnionAll `RegisteredDataItem{kind}` and
# a concrete `RegisteredDataItem{kind,D}` remembered from a live item.
reconstruct(::Type{<:RegisteredDataItem{K}}, id::AbstractString, data, metadata::Dict) where {K} =
    RegisteredDataItem{K}(String(id), "", AbstractCollection[], data, metadata_dict(metadata))

"""
One collection level of a registered path, identified by its name.

`metadata` carries whatever the source attached to this level — `metadata.txt` entries, say — which
is why it survives the round trip through `reconstruct` rather than being rebuilt from the name.
"""
struct NamedCollection <: AbstractCollection
    name::String
    metadata::MetadataDict
end

NamedCollection(name::AbstractString; metadata::AbstractDict=MetadataDict()) =
    NamedCollection(String(name), metadata_dict(metadata))

id(level::NamedCollection)::String = level.name
label(level::NamedCollection)::String = level.name
metadata(level::NamedCollection)::MetadataDict = level.metadata
reconstruct(::Type{NamedCollection}, identity::AbstractString, metadata::Dict) =
    NamedCollection(String(identity); metadata)
Base.:(==)(left::NamedCollection, right::NamedCollection)::Bool = left.name == right.name
Base.isequal(left::NamedCollection, right::NamedCollection)::Bool =
    isequal(left.name, right.name)
Base.hash(level::NamedCollection, seed::UInt)::UInt = hash(level.name, seed)

"""Wrap a registered path's names in collection levels."""
named_collection_path(names::AbstractVector{<:AbstractString})::Vector{AbstractCollection} =
    AbstractCollection[NamedCollection(name) for name in names]
