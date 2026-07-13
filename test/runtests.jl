using Test
using DataBrowserAPI
using DataBrowserCache
using DataBrowserCore
using DataBrowserGUI
using DataBrowserSources

const TEST_CACHE_DEPOT = mktempdir()
pushfirst!(DEPOT_PATH, TEST_CACHE_DEPOT)
atexit(() -> rm(TEST_CACHE_DEPOT; force=true, recursive=true))

include("test_project.jl")

@testset "DataBrowser Tests" begin
    include("test_tase_analysis.jl")
    include("test_project_view_state.jl")
    include("test_directory_source.jl")
    include("test_work_graph.jl")
    include("test_hierarchy_edit.jl")
    include("test_open_options.jl")
    include("test_scan_profile.jl")
    include("test_wide_cache.jl")
    include("test_metadata_pipeline.jl")
    include("test_metadata_contract.jl")
    include("test_registration_contract.jl")
    include("test_event_driven_workspace.jl")
    include("test_profiling.jl")
    include("test_registry_plot.jl")
    include("test_type_api.jl")
    include("test_table_inspector.jl")
    include("test_performance_window.jl")
    include("test_debug_tools.jl")
    include("test_gui_extensions.jl")
    include("test_plots_extension_boundary.jl")
    include("test_annotations.jl")
    get(ENV, "MB_GLFW_SCAN_STRESS", "") == "1" &&
        include("test_glfw_scan_startup_stress.jl")
end
