module Annotations

include("Coords.jl")
include("Layout.jl")
include("Tags.jl")
include("Notes.jl")

using .Coords
using .Layout
using .Tags
using .Notes

export Coords, Layout, Tags, Notes

end # module
