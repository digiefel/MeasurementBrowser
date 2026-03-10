using Test

@testset "MeasurementBrowser Tests" begin
    include("test_analysis_core.jl")
    include("test_wakeup.jl")
    include("test_pund_fatigue.jl")
    include("test_scan_directory_progress.jl")
    include("test_scan_job_controller.jl")
    include("test_plot_job_controller.jl")
    include("test_gui_helpers.jl")
end
