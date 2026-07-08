"""Engine: item index and project callback methods (phase 1 — before cache/workspace)."""
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
using DataBrowserSources
using DataBrowserProfiling

include("ItemIndex.jl")
include("project_engine.jl")

end
