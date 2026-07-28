abstract type AbstractDataItem end

"""One concrete level in an item's collection hierarchy."""
abstract type AbstractCollection end

"""
    id(collection)::String

What identifies one collection level, and the only thing besides its metadata that survives to
rebuild it. There is no default: a type whose job is to identify a grouping has to say what
identifies it.

It is combined with the parent occurrence ID and the concrete type into a one-way digest, so the
digest cannot give it back — `reconstruct(::Type{T}, id, metadata)` receives this value verbatim.
"""
id(collection::AbstractCollection)::String = error(
    "Collection type $(typeof(collection)) must implement id(::$(typeof(collection)))::String")

"""Human-readable label for one collection level."""
label(collection::AbstractCollection)::String = string(collection)

"""Return metadata supplied directly by one collection level."""
metadata(::AbstractCollection)::Dict = Dict()

"""
    id(item)::String

What identifies one item. Stored, displayed in messages, and used for selection, annotation, and
cache lookup exactly as returned — nothing wraps or namespaces it.

There is no default, and no uniqueness is inferred: two items answering the same `id` collide, and
that is a project error the engine reports. Source items and collections answer the same contract
the same way.
"""
id(item::AbstractDataItem)::String = error(
    "Item type $(typeof(item)) must implement id(::$(typeof(item)))::String")

"""Human-readable label for an item. An empty value uses a source-derived label."""
label(::AbstractDataItem)::String = ""

"""Internal item category. Custom item types default to their type name."""
kind(item::AbstractDataItem)::Symbol = Symbol(nameof(typeof(item)))

"""Return an item's complete root-to-leaf path of concrete collection values."""
collection(::AbstractDataItem)::Vector{AbstractCollection} = AbstractCollection[]

"""Return metadata supplied directly by a data item. The default is an empty `Dict`."""
metadata(::AbstractDataItem)::Dict = Dict()

"""The data represented by an item. A custom item is its own data by default."""
item_data(item::AbstractDataItem) = item

"""Process an item. Optional; default identity."""
process(item::AbstractDataItem) = item

"""Analyze a processed item into additional metadata. Optional; default empty `Dict`."""
analyze(::AbstractDataItem)::Dict = Dict()

"""
Whether a payload value can be stored natively by the data cache. Tables are first-class: by
default anything implementing the Tables.jl interface is cacheable, and the cache still requires
storable column types at write time. A type can opt out (or a non-tabular type opt in) by dispatch.
"""
cacheable_data(data)::Bool = Tables.istable(data)

# ---------------------------------------------------------------------------
# Internal workspace hooks
# ---------------------------------------------------------------------------

"""
Let an item adopt the normalized record interpretation produced for it, and the collection path the
index holds for it. Internal workspace hook; the default keeps the item unchanged, because a typed
item derives its own path from its own state. Package-owned carriers adopt both.
"""
attach_record(item::AbstractDataItem, record, path::AbstractVector=AbstractCollection[]) = item
