"""
Tags, notes, coordinates, and layout persistence keyed on item and collection identity.
"""
module DataBrowserAnnotations

using DataBrowserAPI: AbstractDataItem

include("identity.jl")
include("coords.jl")
include("layout.jl")
include("tags.jl")
include("notes.jl")

export Coords, Layout, Tags, Notes
export AnnotationKey, item_annotation_key, collection_annotation_key, ancestor_annotation_keys

end # module
