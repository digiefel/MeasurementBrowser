module MeasurementBrowser

# Core modules
include("scan/Types.jl")
include("scan/Indexer.jl")
include("DeviceParser.jl")
include("projects/RuO2Project.jl")   # defines methods + registers RUO2_PROJECT
include("projects/TASEProject.jl")   # defines methods + registers TASE_PROJECT
include("Gui.jl")

export start_browser, scan_directory
export IndexedCsvFile, MeasurementItem, index_csv_file
export MeasurementHierarchy, HierarchyNode, MeasurementInfo, DeviceInfo
export AbstractProject, RuO2Project, TASEProject, RUO2_PROJECT, TASE_PROJECT
export parse_device_info, detect_kind, kind_label, display_label
export expand_measurement, figure_for_file, figure_for_files
export combined_plot_types, compatible_kinds
export project_name, project_description

end # module MeasurementBrowser
