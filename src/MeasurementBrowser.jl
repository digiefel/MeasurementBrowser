module MeasurementBrowser

using PrecompileTools
using Annotations

# Core modules
include("scan/Types.jl")
include("scan/Indexer.jl")
include("AnalysisResult.jl")
include("DeviceParser.jl")
include("FigureScripts.jl")
include("projects/RuO2Project.jl")   # defines methods + registers RUO2_PROJECT
include("projects/TASEProject.jl")   # defines methods + registers TASE_PROJECT
include("ProjectCache.jl")
include("projects/RuO2/Cache.jl")
include("Precompile.jl")
include("PlotJobs.jl")
include("Gui.jl")

export start_browser, scan_source
export SourceFile, SourceScan, FileFingerprint, index_source_file
export AnalysisResult
export MeasurementHierarchy, HierarchyNode, MeasurementInfo, DeviceInfo
export AbstractProject, RuO2Project, TASEProject, RUO2_PROJECT, TASE_PROJECT
export MeasurementFilterClause, MeasurementGroupFilter, NamedMeasurementGroup
export FigureMeasurement, FigureScriptData, prepare_figure_script_data, infer_measurement_group
export parse_device_info, detect_kind, kind_label, display_label
export measurements_for_file
export load_plot_for_file, analyze_plot_for_file, draw_plot_for_file
export load_plot_for_files, analyze_plot_for_files, draw_plot_for_files
export debug_plot
export available_analyses, run_analysis, draw_analysis_view
export combined_plot_types, compatible_kinds
export project_name, project_description
export ProjectCacheIdentity, ProjectCacheSnapshot, ProjectCacheStatus, ProjectCacheFileError
export ProjectCacheUnsupportedError, ProjectCacheMissingError, ProjectCacheInvalidError
export ProjectCacheBuildError
export write_project_cache!, load_project_cache, cache_status
export project_cache_identity, project_cache_path, project_cache_id, new_project_cache_id

end # module MeasurementBrowser
