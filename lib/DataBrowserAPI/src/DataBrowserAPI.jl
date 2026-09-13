"""
Payload-agnostic project and item contracts for DataBrowser.

The project contract, the typed pipeline stage contract, the source and item contracts, and the
shared item data model (`ItemIndex`: item records, the collection hierarchy, and source scans) live
here. The registration dialect lives in `DataBrowserRecipes`, plot registration and rendering in
`DataBrowserPlots`; both are ordinary clients of what this package declares.
"""
module DataBrowserAPI

using Dates
using SHA

# `read` is the pipeline's source-touching stage, and that is exactly what `Base.read` means.
# Extending it keeps `using DataBrowser` from shadowing a function every Julia user already knows.
import Base: read

include("project_contract.jl")
include("metadata_types.jl")
include("source_contract.jl")
include("item_contract.jl")
include("collection_id.jl")
include("stage_contract.jl")
include("ItemIndex.jl")
include("timing_debug.jl")

end
