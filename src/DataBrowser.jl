"""DataBrowser umbrella: public API re-exports and default wiring."""
module DataBrowser

using PrecompileTools
using DataBrowserAPI
using DataBrowserSources
using DataBrowserCore
using DataBrowserPlots
using DataBrowserGUI

using DataBrowserAPI:
    AbstractDataItem,
    AbstractCollection,
    PlotKind,
    RegisteredPlot,
    cacheable,
    collection,
    define_project,
    display_label,
    fingerprint,
    id,
    item_data,
    item_label,
    kind,
    metadata,
    plot_kind_from_name,
    plot_kind_label,
    plot_kind_name,
    process,
    register_collection_analysis!,
    register_item!,
    register_plot!,
    registered_plot_kinds,
    close_source!,
    open_source,
    scan_profile_summary,
    scan_source_profile,
    source_items
import DataBrowserAPI: plot_data!, setup_plot

using DataBrowserSources: DirectorySource, SourceFile, inspect_table
using DataBrowserCore.ItemIndex: DataItem
using DataBrowserCore: items_for_file
using DataBrowserCore.Workspace:
    close_workspace!,
    open_workspace,
    query_items,
    read_item_data,
    select_items!
using DataBrowserGUI: open_browser

export open_browser,
    open_workspace,
    close_workspace!,
    select_items!,
    query_items,
    read_item_data,
    define_project,
    register_item!,
    register_collection_analysis!,
    register_plot!,
    setup_plot,
    plot_data!,
    RegisteredPlot,
    AbstractDataItem,
    AbstractCollection,
    id,
    item_label,
    kind,
    collection,
    metadata,
    item_data,
    process,
    cacheable,
    fingerprint,
    Project,
    DirectorySource,
    SourceFile,
    DataItem,
    inspect_table,
    items_for_file,
    PlotKind,
    display_label,
    plot_kind_from_name,
    plot_kind_label,
    plot_kind_name,
    registered_plot_kinds,
    open_source,
    close_source!,
    source_items,
    scan_profile_summary,
    scan_source_profile

using DataBrowserAPI: Project, project_name, registered_plot_kinds
import DataBrowserCore.Cache as Cache
import DataBrowserCore.Workspace as Workspace
import DataBrowserGUI: Browser

include("Precompile.jl")

end
