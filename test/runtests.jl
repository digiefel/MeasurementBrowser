using Test

const TEST_CACHE_DEPOT = mktempdir()
pushfirst!(DEPOT_PATH, TEST_CACHE_DEPOT)
atexit(() -> rm(TEST_CACHE_DEPOT; force=true, recursive=true))

include("test_project.jl")

@testset "MeasurementBrowser Tests" begin
    include("test_tase_analysis.jl")
    include("test_project_view_state.jl")
    include("test_scan_directory_progress.jl")
    include("test_scan_profile.jl")
    include("test_profiling.jl")
    include("test_registry_plot.jl")
    include("test_type_api.jl")
    include("test_project_cache.jl")
    include("test_table_inspector.jl")
    include("test_annotations.jl")
    if Threads.nthreads() >= 4
        include("test_threaded_crash_regressions.jl")
    else
        @info "Skipping threaded workspace smoke test. Run `julia --project --threads=4 -e 'using Pkg; Pkg.test()'` to exercise threaded scanning."
    end
    get(ENV, "MB_GLFW_SCAN_STRESS", "") == "1" &&
        include("test_glfw_scan_startup_stress.jl")
end
