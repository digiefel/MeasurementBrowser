using Test

@testset "MeasurementBrowser Tests" begin
    include("test_wakeup.jl")
    include("test_scan_directory_progress.jl")
    include("test_scan_job_controller.jl")
end
