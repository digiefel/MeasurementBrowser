abstract type AbstractDataItem end

"""One concrete level in an item's collection hierarchy."""
abstract type AbstractCollection end

"""
Value used to derive one collection level's deterministic occurrence ID.

The default uses the complete concrete collection value. Override this only when the value contains
state that is deliberately not part of the collection ID or cannot be canonically encoded.
"""
id(collection::AbstractCollection) = collection

"""Human-readable label for one collection level."""
label(collection::AbstractCollection)::String = string(collection)

"""Return metadata supplied directly by one collection level."""
metadata(::AbstractCollection)::Dict = Dict()

"""Stable identity of an item. An empty value lets DataBrowser assign one."""
id(::AbstractDataItem) = ""

"""Human-readable label for an item. An empty value uses a source-derived label."""
item_label(::AbstractDataItem)::String = ""

"""Internal item category. Custom item types default to their type name."""
kind(item::AbstractDataItem)::Symbol = Symbol(nameof(typeof(item)))

"""Return an item's collection path."""
collection(::AbstractDataItem) = String[]

"""Return metadata supplied directly by a data item. The default is an empty `Dict`."""
metadata(::AbstractDataItem)::Dict = Dict()

"""The data represented by an item. A custom item is its own data by default."""
item_data(item::AbstractDataItem) = item

"""Process an item. Optional; default identity."""
process(item::AbstractDataItem) = item

"""Analyze a processed item into additional metadata. Optional; default empty `Dict`."""
analyze(::AbstractDataItem)::Dict = Dict()

"""Whether an item's data should be persisted by the data cache. Optional; default `false`."""
cacheable(::AbstractDataItem)::Bool = false

"""
Whether a payload value can be stored natively by the data cache. Tables are first-class: by
default anything implementing the Tables.jl interface is cacheable, and the cache still requires
storable column types at write time. A type can opt out (or a non-tabular type opt in) by dispatch.
"""
cacheable_data(data)::Bool = Tables.istable(data)

fingerprint(::AbstractDataItem) = nothing

"""Collection metadata contributed by a source for one collection path."""
collection_metadata(::AbstractDataSource, ::AbstractVector{<:AbstractString})::Dict{Symbol,Any} =
    Dict{Symbol,Any}()

"""Whether a source supplied collection metadata for this scan."""
has_collection_metadata(::AbstractDataSource)::Bool = false

# ---------------------------------------------------------------------------
# Internal workspace hooks
# ---------------------------------------------------------------------------

"""Optional rewrite of a collection's members (one output per input). Internal workspace hook."""
function _process_collection end

"""Optional fold over a collection's members into collection-node metadata. Internal workspace hook."""
function _analyze_collection end

"""Per-item analysis metadata computed after indexing. Internal workspace hook."""
function _analyze_item end

"""Whether one item kind has a registered collection `process` stage."""
function _has_collection_process end

"""Whether one item kind has a registered collection `analyze` stage."""
function _has_collection_analysis end
