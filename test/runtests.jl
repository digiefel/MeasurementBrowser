using Test

@testset "MeasurementBrowser Tests" begin
    include("test_analysis_core.jl")
    include("test_tase_analysis.jl")
    include("test_wakeup.jl")
    include("test_pund_fatigue.jl")
    include("test_ruo2_interpretation.jl")
    include("test_cvsweep.jl")
    include("test_ruo2_simple_plot_api.jl")
    include("test_scan_directory_progress.jl")
    include("test_project_cache.jl")
    include("test_annotations.jl")
end
