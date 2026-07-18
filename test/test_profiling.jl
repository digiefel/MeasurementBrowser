import DataBrowserProfiling as Profiling
using Test

const PROFILE = Profiling

module DebugTimingFixture
import DataBrowserProfiling as Profiling

timed_identity(value) = Profiling.@time_dbg identity(value)
end

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

@testset "explicit debug timings" begin
    @test DebugTimingFixture.timed_identity(9) == 9

    timings = PROFILE.DebugTimings(start_ns=time_ns())
    value = PROFILE.with_debug_timings(timings) do
        DebugTimingFixture.timed_identity(42)
    end
    @test value == 42
    rows = PROFILE.debug_timing_rows(timings)
    @test only(rows).label == "identity"
    @test only(rows).calls == 1
    @test timings.stop_ns >= timings.start_ns

    PROFILE.reset_debug_timings!(timings)
    PROFILE.with_debug_timings(timings) do
        DebugTimingFixture.timed_identity(1)
        DebugTimingFixture.timed_identity(2)
    end
    rows = PROFILE.debug_timing_rows(timings)
    @test only(rows).calls == 2

    requests = Channel{Int}(1)
    responses = Channel{Int}(1)
    worker = @async begin
        for request in requests
            put!(responses, DebugTimingFixture.timed_identity(request))
        end
    end
    PROFILE.reset_debug_timings!(timings)
    PROFILE.with_debug_timings(timings) do
        put!(requests, 7)
        @test take!(responses) == 7
    end
    close(requests)
    wait(worker)
    rows = PROFILE.debug_timing_rows(timings)
    @test only(rows).calls == 1

    output = mktempdir()
    PROFILE.write_debug_timings(output, timings)
    @test isfile(joinpath(output, "debug_timings.txt"))
    @test isfile(joinpath(output, "debug_timings.csv"))
    text = read(joinpath(output, "debug_timings.txt"), String)
    @test occursin("Measurement wall time", text)
    @test occursin("identity", text)
    @test !occursin("% measured", text)
end
