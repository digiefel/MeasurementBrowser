module MeasurementBrowser

using PrecompileTools
using Annotations

include("Profiling.jl")

include("Projects.jl")
using .Projects:
    AbstractDataSource,
    AbstractDataSourceItem,
    AbstractDataItem,
    Collection,
    CollectionRecipe,
    DEFAULT_PROJECT,
    ItemRecipe,
    KindProfileRow,
    PROJECTS,
    PlotRecipe,
    Project,
    SourceChanges,
    SourceError,
    SourceItemProfile,
    SourceProfileRow,
    cacheable,
    close_source!,
    collection,
    data_items,
    fingerprint,
    id,
    item_data,
    item_label,
    metadata,
    open_source,
    source_id,
    source_item_id,
    source_item_label,
    source_item_path,
    source_item_timestamp,
    source_items,
    source_label,
    watch_source
import .Projects:
    detect_kind,
    collection_path_label,
    display_label,
    kind,
    kind_label,
    process,
    project_description,
    project_name,
    record_scan_phase!,
    reset_scan_profile!,
    scan_profile_summary,
    scan_source_profile,
    finish_source_profile!

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
    is_job_cancelled,
    item_timestamp_key,
    with_cancel

include("Cache.jl")

include("Workspace.jl")
using .Workspace:
    close_workspace!,
    open_workspace,
    query_items,
    read_item_data,
    select_items!

include("DataSources/DirectorySource.jl")

const ENGINE_ONLY_BENCHMARK_LOAD = get(ENV, "MB_BENCH_ENGINE_ONLY", "0") == "1"

if !ENGINE_ONLY_BENCHMARK_LOAD
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
        plot_data!,
        setup_plot
end

include("TableInspector.jl")
using .TableInspector: TablePreview, inspect_table

include("Project.jl")
if ENGINE_ONLY_BENCHMARK_LOAD
    open_browser(args...; kwargs...) =
        error("open_browser is unavailable when MB_BENCH_ENGINE_ONLY=1")
else
    include("Browser.jl")
    using .Browser: open_browser
    include("Precompile.jl")
end

# Public API — exactly the contract. Everything else (engine internals: ItemRecord, the plot
# dispatch, scan/data-access plumbing, table inspector, project introspection) stays un-exported and
# is reachable only as `MeasurementBrowser.name` when genuinely needed.
#
# Run the app
export open_browser, open_workspace, close_workspace!, select_items!, query_items
# Build a project
export define_project, register_item!, register_collection_analysis!, register_plot!
# The AbstractDataItem contract — implement these for a custom item type
export AbstractDataItem, Collection
export id, item_label, kind, collection, metadata, item_data, process, cacheable, fingerprint
# Types you name
export Project, DirectorySource, SourceFile, DataItem

end # module MeasurementBrowser
