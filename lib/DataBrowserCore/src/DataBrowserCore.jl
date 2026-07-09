"""Headless workspace engine: item index, cache, work graph, and workspace lifecycle."""
module DataBrowserCore

using DataBrowserAPI
using DataBrowserAPI:
    CollectionRecipe,
    ItemRecipe,
    PlotRecipe,
    Project,
    SourceItemProfile,
    data_items,
    detect_kind,
    kind_label,
    source_item_id,
    source_item_label,
    source_item_path
using DataBrowserAnnotations
using DataBrowserSources
import DataBrowserProfiling as Profiling

include("ItemIndex.jl")
include("project_engine.jl")
include("WorkGraph.jl")
include("Cache.jl")
include("Workspace.jl")

end
