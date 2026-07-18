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
    AbstractDataSource,
    AbstractDataSourceItem,
    AbstractCollection,
    analyze,
    cacheable,
    collection,
    define_project,
    data_items,
    display_label,
    fingerprint,
    id,
    item_data,
    item_label,
    kind,
    label,
    metadata,
    process,
    register_collection_analysis!,
    register_item!,
    close_source!,
    open_source,
    scan_profile_summary,
    scan_source_profile,
    source_items,
    source_id,
    source_label,
    source_item_id,
    source_item_label,
    source_item_path,
    source_item_timestamp,
    source_open_options,
    watch_source
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

using DataBrowserSources: DirectorySource, SourceFile, inspect_table
using DataBrowserCore: items_for_file
using DataBrowserCore.Workspace:
    close_workspace!,
    materialize_items,
    open_workspace,
    query_items,
    read_item_data,
    select_items!,
    wait_workspace_idle!,
    workspace_status
using DataBrowserGUI: open_browser, close_browser!, BrowserSession

export open_browser,
    close_browser!,
    BrowserSession,
    open_workspace,
    close_workspace!,
    select_items!,
    query_items,
    materialize_items,
    read_item_data,
    wait_workspace_idle!,
    workspace_status,
    define_project,
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
    item_label,
    kind,
    label,
    collection,
    metadata,
    item_data,
    process,
    analyze,
    cacheable,
    fingerprint,
    Project,
    DirectorySource,
    SourceFile,
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
    source_id,
    source_label,
    source_items,
    source_item_id,
    source_item_label,
    source_item_path,
    source_item_timestamp,
    source_open_options,
    watch_source,
    data_items,
    scan_profile_summary,
    scan_source_profile

using DataBrowserAPI: Project, project_name
import DataBrowserCache as Cache
import DataBrowserCore.Workspace as Workspace
import DataBrowserGUI: Browser

include("Precompile.jl")

end
