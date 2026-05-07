using Test

@testset "MeasurementBrowser Tests" begin
    include("test_analysis_core.jl")
    include("test_tase_analysis.jl")
    include("test_wakeup.jl")
    include("test_pund_fatigue.jl")
    include("test_cvsweep.jl")
    include("test_scan_directory_progress.jl")
    include("test_project_cache.jl")
    include("test_plot_job_controller.jl")
    include("test_annotations.jl")
    include("test_figure_scripts.jl")
    include("test_gui_helpers.jl")
end
