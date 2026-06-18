module MeasurementBrowser

using PrecompileTools
using Annotations

include("Projects.jl")
using .Projects:
    AbstractDataSource,
    AbstractDataSourceItem,
    AbstractDataItem,
    Collection,
    CollectionStatRecipe,
    DEFAULT_PROJECT,
    ItemRecipe,
    KindProfile,
    PROJECTS,
    PlotRecipe,
    Project,
    cacheable,
    close_source!,
    collection,
    collection_stats,
    data_items,
    item_data,
    item_fingerprint,
    item_id,
    item_label,
    load_data_item,
    open_source,
    parameters,
    source_fingerprint,
    source_id,
    source_item_fingerprint,
    source_item_id,
    source_item_label,
    source_item_path,
    source_item_timestamp,
    source_items,
    source_label,
    stats
import .Projects:
    detect_kind,
    collection_path_label,
    display_label,
    kind,
    kind_label,
    process,
    project_description,
    project_name,
    reset_scan_profile!,
    scan_profile_summary

include("ItemIndex.jl")
using .ItemIndex:
    CANCEL_CALLBACK_KEY,
    DataItem,
    HierarchyNode,
    JobCancelled,
    ItemFailure,
    Hierarchy,
    ItemRecord,
    SourceScan,
    check_cancel,
    children,
    collection_path_key,
    collection_path_tuple,
    emit_progress,
    insert_item!,
    item_record_key,
    is_job_cancelled,
    item_timestamp_key,
    scan_source,
    with_cancel

include("Cache.jl")

include("Workspace.jl")
using .Workspace:
    close_workspace!,
    open_workspace,
    read_item_data,
    select_items!

include("DataSources/DirectorySource.jl")

include("Visualization.jl")
using .Visualization:
    PlotKind,
    RegisteredPlot,
    plot_kind_from_name,
    plot_kind_label,
    plot_kind_name,
    plot_kind_symbol,
    plot_kinds,
    registered_plot_kinds
import .Visualization:
    debug_plot,
    plot_data!,
    setup_plot

include("TableInspector.jl")
using .TableInspector: TablePreview, inspect_table

include("Project.jl")
include("Precompile.jl")
include("Browser.jl")
using .Browser: open_browser

# Public API — exactly the contract. Everything else (engine internals: ItemRecord, the plot
# dispatch, scan/data-access plumbing, table inspector, project introspection) stays un-exported and
# is reachable only as `MeasurementBrowser.name` when genuinely needed.
#
# Run the app
export open_browser, open_workspace, close_workspace!, select_items!
# Build a project
export define_project, register_item!, register_collection_stat!, register_plot!
# The AbstractDataItem contract — implement these for a custom item type
export AbstractDataItem, Collection
export item_id, item_label, kind, collection, parameters, stats, item_data, process, cacheable
# Types you name
export Project, DirectorySource, SourceFile, DataItem

end # module MeasurementBrowser
