"""
Payload-agnostic project and item contracts for DataBrowser.

Types, generic function declarations, project construction (`define_project`, `register_item!`,
`register_collection_analysis!`), and the shared item data model (`ItemIndex`: item records, the
collection hierarchy, and source scans) live here. Plot registration and rendering live in
`DataBrowserPlots`.
"""
module DataBrowserAPI

using Dates
using SHA
import Tables

include("project_types.jl")
include("metadata_types.jl")
include("source_contract.jl")
include("item_contract.jl")
include("collection_id.jl")
include("interface.jl")
include("construction.jl")
include("ItemIndex.jl")

end
