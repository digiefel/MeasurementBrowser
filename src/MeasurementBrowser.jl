module MeasurementBrowser

using PrecompileTools
using Annotations

# Core modules
include("scan/Types.jl")
include("scan/Indexer.jl")
include("AnalysisResult.jl")
include("DeviceParser.jl")
include("FigureScripts.jl")
include("DataAccess.jl")
include("projects/RuO2Project.jl")   # defines methods + registers RUO2_PROJECT
include("projects/TASEProject.jl")   # defines methods + registers TASE_PROJECT
include("ProjectCache.jl")
include("Precompile.jl")
include("PlotJobs.jl")
include("Gui.jl")

export start_browser, scan_source
export SourceFile, SourceScan, FileFingerprint, MeasurementAnalysisFailure, index_source_file
export AnalysisResult
export MeasurementHierarchy, HierarchyNode, MeasurementInfo, DeviceInfo
export AbstractProject, RuO2Project, TASEProject, RUO2_PROJECT, TASE_PROJECT
export PlotKind, TASEFourTerminalIVPlot
export RuO2PUNDPlot, RuO2IVSweepPlot, RuO2TLM4PointPlot, RuO2CVSweepPlot
export RuO2TLMAnalysisPlot, RuO2TLMTemperaturePlot, RuO2PUNDFatiguePlot
export MeasurementFilterClause, MeasurementGroupFilter, NamedMeasurementGroup
export FigureMeasurement, FigureScriptData, prepare_figure_script_data, infer_measurement_group
export parse_device_info, detect_kind, kind_label, display_label
export measurements_for_file
export load_source_data, data_of_measurements, setup_plot, plot_data!
export available_plot_kinds, default_plot_kind, plot_kind_label, plot_kind_description
export plot_kind_measurement_kinds, plot_kind_min_measurements, supports_plot_kind
export debug_plot
export available_analyses, run_analysis, draw_analysis_view
export project_name, project_description
export ProjectCacheIdentity, ProjectCacheSnapshot, ProjectCacheStatus, ProjectCacheFileError
export ProjectCacheMissingError, ProjectCacheInvalidError
export write_project_cache!, load_project_cache, cache_status
export project_cache_identity, project_cache_path, project_cache_id, new_project_cache_id

end # module MeasurementBrowser
