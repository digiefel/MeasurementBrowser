module MeasurementBrowser

using PrecompileTools
using Annotations

include("Projects.jl")
using .Projects:
    AbstractDataItem,
    Collection,
    DEFAULT_PROJECT,
    DeviceStatRecipe,
    KindProfile,
    MeasurementRecipe,
    PROJECTS,
    PlotRecipe,
    Project,
    cacheable,
    collection,
    item_data,
    item_id,
    item_label,
    parameters,
    read_data,
    stats
import .Projects:
    compute_and_add_item_stats!,
    detect_kind,
    collection_path_label,
    display_label,
    interpret_file,
    kind,
    kind_label,
    load_source_data,
    parse_collection_info,
    process,
    process_item_data,
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
    SourceFile,
    SourceScan,
    build_clean_title,
    check_cancel,
    children,
    collect_source_files,
    collection_metadata_path,
    collection_path_key,
    collection_path_tuple,
    emit_progress,
    file_fingerprint,
    has_collection_metadata,
    index_source_file,
    insert_item!,
    interpret_items,
    is_job_cancelled,
    load_scan_metadata,
    item_timestamp_key,
    items_for_file,
    parse_timestamp,
    scan_source,
    with_cancel

include("Cache.jl")

include("Workspace.jl")
using .Workspace:
    close_workspace!,
    open_workspace,
    read_item_data,
    select_measurements!

include("Visualization.jl")
using .Visualization:
    PlotKind,
    RegistryPlot,
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

export open_browser
export SourceFile, ItemRecord, ItemFailure
export AbstractDataItem, DataItem, Collection
export item_id, item_label, kind, collection, parameters, stats, item_data, read_data, process, cacheable
export PlotKind, plot_kinds
export parse_collection_info, detect_kind, kind_label, display_label, interpret_file
export items_for_file
export load_source_data, read_item_data, process_item_data, setup_plot, plot_data!
export compute_and_add_item_stats!
export debug_plot
export TablePreview, inspect_table
export project_name, project_description
export open_workspace, close_workspace!, select_measurements!
export define_project, register_measurement!, register_device_stat!, register_plot!
export Project

end # module MeasurementBrowser
