"""DataBrowser umbrella: public API re-exports and default wiring."""
module DataBrowser

using PrecompileTools
using DataBrowserAPI
using DataBrowserRecipes
using DataBrowserSources
using DataBrowserCore
using DataBrowserPlots
using DataBrowserGUI

using DataBrowserAPI:
    AbstractDataItem,
    AbstractDataSource,
    AbstractDataSourceItem,
    AbstractCollection,
    AbstractProject,
    analyze,
    collection,
    entries,
    display_label,
    fingerprint,
    id,
    item_data,
    label,
    metadata,
    process,
    close_source!,
    open_source,
    source_items,
    source_id,
    source_label,
    source_item_path,
    source_item_timestamp,
    watch_source,
    reconstruct,
    read
using DataBrowserRecipes:
    Project,
    define_project,
    register_collection_analysis!,
    register_csv!,
    register_item!
using DataBrowserPlots:
    PlotKind,
    RegisteredPlot,
    plot_data!,
    plot_kind_from_name,
    plot_kind_label,
    plot_kind_name,
    register_plot!,
    registered_plot_kinds,
    setup_plot

using DataBrowserSources: DirectorySource, SourceFile
using DataBrowserCore: items_for_file
using DataBrowserCore.Workspace:
    close_workspace!,
    materialize_items,
    modify_workspace!,
    open_workspace,
    query_items,
    read_item_data,
    select_items!,
    wait_workspace_idle!,
    workspace_status
using DataBrowserGUI: open_browser, close_browser!, BrowserSession, gui_timings, reset_timings!
using DataBrowserGUI: wait_browser_ready
export wait_browser_ready

export open_browser,
    close_browser!,
    BrowserSession,
    gui_timings,
    reset_timings!,
    open_workspace,
    modify_workspace!,
    close_workspace!,
    select_items!,
    query_items,
    materialize_items,
    read_item_data,
    wait_workspace_idle!,
    workspace_status,
    define_project,
    register_csv!,
    register_item!,
    register_collection_analysis!,
    register_plot!,
    setup_plot,
    plot_data!,
    RegisteredPlot,
    AbstractDataItem,
    AbstractDataSource,
    AbstractDataSourceItem,
    AbstractCollection,
    id,
    label,
    collection,
    metadata,
    item_data,
    process,
    analyze,
    fingerprint,
    Project,
    DirectorySource,
    SourceFile,
    items_for_file,
    PlotKind,
    display_label,
    plot_kind_from_name,
    plot_kind_label,
    plot_kind_name,
    registered_plot_kinds,
    open_source,
    close_source!,
    source_id,
    source_label,
    source_items,
    source_item_path,
    source_item_timestamp,
    watch_source,
    entries,
    reconstruct,
    AbstractProject

using DataBrowserAPI: project_name
import DataBrowserCache as Cache
import DataBrowserCore.Workspace as Workspace
import DataBrowserGUI: Browser

include("Precompile.jl")

end
