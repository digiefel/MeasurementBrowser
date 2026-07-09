"""
Payload-agnostic project and item contracts for DataBrowser.

Types, generic function declarations, project construction (`define_project`, `register_*`), and
plot-kind identities live here. Engine methods that drive callbacks and touch concrete payloads are
defined in the main package until `DataBrowserCore` is extracted.
"""
module DataBrowserAPI

using Dates
using InteractiveUtils: subtypes

include("project_types.jl")
include("metadata_types.jl")
include("source_contract.jl")
include("item_contract.jl")
include("interface.jl")
include("construction.jl")
include("plots.jl")

end
