using Test

const TEST_CACHE_DEPOT = mktempdir()
pushfirst!(DEPOT_PATH, TEST_CACHE_DEPOT)
atexit(() -> rm(TEST_CACHE_DEPOT; force=true, recursive=true))

include("test_project.jl")

@testset "MeasurementBrowser Tests" begin
    include("test_tase_analysis.jl")
    include("test_project_view_state.jl")
    include("test_directory_source.jl")
    include("test_work_graph.jl")
    include("test_scan_profile.jl")
    include("test_event_driven_workspace.jl")
    include("test_profiling.jl")
    include("test_registry_plot.jl")
    include("test_type_api.jl")
    include("test_table_inspector.jl")
    include("test_performance_window.jl")
    include("test_annotations.jl")
    get(ENV, "MB_GLFW_SCAN_STRESS", "") == "1" &&
        include("test_glfw_scan_startup_stress.jl")
end
