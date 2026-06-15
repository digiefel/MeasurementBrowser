module MeasurementBrowser

using PrecompileTools
using Annotations

include("Projects.jl")
using .Projects:
    AbstractProject,
    DEFAULT_PROJECT,
    PROJECTS
import .Projects:
    compute_and_add_measurement_stats!,
    detect_kind,
    device_path_label,
    display_label,
    interpret_file,
    kind_label,
    load_source_data,
    parse_device_info,
    process_measurement_data,
    project_description,
    project_name,
    reset_scan_profile!,
    scan_profile_summary

include("MeasurementIndex.jl")
using .MeasurementIndex:
    CANCEL_CALLBACK_KEY,
    DeviceInfo,
    HierarchyNode,
    JobCancelled,
    MeasurementAnalysisFailure,
    MeasurementHierarchy,
    MeasurementInfo,
    SourceFile,
    SourceScan,
    build_clean_title,
    check_cancel,
    children,
    collect_source_files,
    device_info_path,
    device_path_key,
    device_path_tuple,
    emit_progress,
    file_fingerprint,
    has_device_metadata,
    index_source_file,
    insert_measurement!,
    interpret_measurements,
    is_job_cancelled,
    load_scan_metadata,
    measurement_timestamp_key,
    measurements_for_file,
    parse_timestamp,
    scan_source,
    with_cancel

include("Cache.jl")

include("Workspace.jl")
using .Workspace:
    close_workspace!,
    open_workspace,
    read_measurement_data,
    select_measurements!

include("Visualization.jl")
using .Visualization:
    PlotKind,
    RegistryPlot,
    plot_kind_from_name,
    plot_kind_name,
    plot_kind_symbol,
    plot_kinds
import .Visualization:
    debug_plot,
    plot_data!,
    setup_plot

include("TableInspector.jl")
using .TableInspector: TablePreview, inspect_table

include("projects/TASEProject.jl")
include("Project.jl")
include("FigureScripts.jl")
include("Precompile.jl")
include("Browser.jl")
using .Browser: open_browser

export open_browser
export AbstractProject
export SourceFile, MeasurementInfo, DeviceInfo, MeasurementAnalysisFailure
export PlotKind, plot_kinds
export parse_device_info, detect_kind, kind_label, display_label, interpret_file
export measurements_for_file
export load_source_data, read_measurement_data, process_measurement_data, setup_plot, plot_data!
export compute_and_add_measurement_stats!
export debug_plot
export TablePreview, inspect_table
export project_name, project_description
export open_workspace, close_workspace!, select_measurements!
export define_project, register_measurement!, register_device_stat!, register_plot!
export Project

end # module MeasurementBrowser
