using Test

include("test_project.jl")

@testset "MeasurementBrowser Tests" begin
    include("test_tase_analysis.jl")
    include("test_project_view_state.jl")
    include("test_scan_directory_progress.jl")
    include("test_scan_profile.jl")
    include("test_registry_plot.jl")
    include("test_project_cache.jl")
    include("test_table_inspector.jl")
    include("test_annotations.jl")
end
