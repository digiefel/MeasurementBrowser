"""
Typed metadata values stored in the item index and DuckDB cache.

Scalars and homogeneous vectors only — the cache persists these as typed columns.
"""
const MetadataValue = Union{
    Bool, Int64, Float64, String, Symbol, Date, DateTime, Missing,
    Vector{Bool}, Vector{Int64}, Vector{Float64}, Vector{String},
}

"""The typed dict the engine stores for item/collection `metadata`."""
const MetadataDict = Dict{Symbol,MetadataValue}
