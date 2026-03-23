module MeasurementBrowser

using PrecompileTools

# Core modules
include("scan/Types.jl")
include("scan/Indexer.jl")
include("AnalysisResult.jl")
include("DeviceParser.jl")
include("BadRegistry.jl")
include("projects/RuO2Project.jl")   # defines methods + registers RUO2_PROJECT
include("projects/TASEProject.jl")   # defines methods + registers TASE_PROJECT
include("Precompile.jl")
include("Gui.jl")

export start_browser, scan_directory
export IndexedCsvFile, MeasurementItem, index_csv_file
export AnalysisResult
export MeasurementHierarchy, HierarchyNode, MeasurementInfo, DeviceInfo
export AbstractProject, RuO2Project, TASEProject, RUO2_PROJECT, TASE_PROJECT
export parse_device_info, detect_kind, kind_label, display_label
export expand_measurement, load_plot_for_file, analyze_plot_for_file, draw_plot_for_file
export load_plot_for_files, analyze_plot_for_files, draw_plot_for_files
export available_analyses, run_analysis, draw_analysis_view
export combined_plot_types, compatible_kinds
export project_name, project_description

end # module MeasurementBrowser
