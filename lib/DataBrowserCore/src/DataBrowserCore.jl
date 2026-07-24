"""
Headless workspace engine: work graph and workspace lifecycle over the `DataBrowserAPI` item index
and the `DataBrowserCache` store.
"""
module DataBrowserCore

using DataBrowserAPI
using DataBrowserAPI:
    CollectionRecipe,
    ItemRecipe,
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
using DataBrowserAPI: @timed_dbg

include("project_engine.jl")
include("WorkGraph.jl")
include("Workspace.jl")
include("TableModel.jl")

export InspectorTable, merge_item_tables

end
