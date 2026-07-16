using DataBrowserAPI: AbstractDataItem, id

"""Stable storage key for annotations — item id, registration path, or typed node identity."""
const AnnotationKey = String

"""Annotation key for one item, from [`id`](@ref)."""
item_annotation_key(item::AbstractDataItem) = String(id(item))

"""Annotation key for a collection path given as nested segment names."""
collection_annotation_key(segments::AbstractVector{<:AbstractString}) = join(segments, '/')

collection_annotation_key(key::AbstractString) = String(key)

"""Parent collection keys for inheritance lookups on slash-joined paths."""
function ancestor_annotation_keys(path::AnnotationKey)::Vector{AnnotationKey}
    parts = split(path, '/')
    length(parts) <= 1 && return AnnotationKey[]
    out = AnnotationKey[]
    for i in 1:(length(parts) - 1)
        push!(out, join(parts[1:i], '/'))
    end
    return out
end
