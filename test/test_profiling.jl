import DataBrowserProfiling as Profiling
using Test

const PROFILE = Profiling

@testset "process memory diagnostics" begin
    @test PROFILE.process_rss_bytes() > 0
end

@testset "CPU sampling diagnostics" begin
    PROFILE.start_sampling!()
    result = try
        deadline = time() + 0.05
        value = 0.0
        while time() < deadline
            for index in 1:10_000
                value += sin(index)
            end
        end
        @test isfinite(value)
        PROFILE.stop_sampling!()
    finally
        PROFILE.cancel_sampling!()
    end

    @test result.total_samples > 0
    @test result.delay_seconds > 0
    @test !isempty(result.rows)
    @test all(row -> 0 <= row.self_samples <= row.samples, result.rows)
    @test issorted(result.rows; by=row -> (row.self_samples, row.samples), rev=true)
end
