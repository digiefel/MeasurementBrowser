abstract type AbstractDataItem end

"""Nestable container an item is placed in. Future home for collection behaviour."""
abstract type AbstractCollection end

"""Stable string identity of an item."""
function id end

"""Human-readable label for an item."""
function item_label end

"""The kind symbol the plot registry is keyed on."""
function kind end

"""Canonical tree placement of an item, as nested collection names (`Vector{String}`)."""
function collection end

"""Metadata of an item (`Dict{Symbol,Any}`): parsed parameters, computed values, provenance merged."""
function metadata end

"""The materialized data carried by an item (also reachable as `item.data`)."""
function item_data end

"""Process an item into the item a view consumes. Optional; default identity."""
process(item::AbstractDataItem) = item

"""Whether an item's data should be persisted by the data cache. Optional; default `false`."""
cacheable(::AbstractDataItem)::Bool = false

"""
Whether a payload value can be stored natively by the data cache. Core claims the first-class
payload types (`AbstractDataFrame`); extension packages that provide support for a type add its
method alongside that support. Everything else stays source-backed. The built-in `DataItem`
answers `cacheable` through this trait.
"""
cacheable_data(::Any)::Bool = false

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
