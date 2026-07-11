"""
Payload-agnostic project and item contracts for DataBrowser.

Types, generic function declarations, project construction (`define_project`, `register_*`),
plot-kind identities, and the shared item data model (`ItemIndex`: item records, the collection
hierarchy, and source scans) live here.
"""
module DataBrowserAPI

using Dates
using InteractiveUtils: subtypes
import Tables

include("project_types.jl")
include("metadata_types.jl")
include("source_contract.jl")
include("item_contract.jl")
include("interface.jl")
include("construction.jl")
include("plots.jl")
include("ItemIndex.jl")

end
