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
include("Gui.jl")

export start_browser, scan_source
export AbstractProject
export SourceFile, MeasurementInfo, DeviceInfo
export PlotKind, plot_kinds
export parse_device_info, detect_kind, kind_label, display_label, interpret_file
export measurements_for_file
export load_source_data, read_measurement_data, process_measurement_data, setup_plot, plot_data!
export debug_plot
export project_name, project_description

end # module MeasurementBrowser
