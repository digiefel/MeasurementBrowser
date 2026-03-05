using Test

@testset "MeasurementBrowser Tests" begin
    include("test_wakeup.jl")
    include("test_scan_directory_progress.jl")
end
