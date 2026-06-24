const PERFORMANCE_BROWSER = MeasurementBrowser.Browser

@testset "performance window live figures" begin
    live = PERFORMANCE_BROWSER.LivePlotsState(capacity=3)
    timings = Dict{Symbol,Vector{Float64}}(
        :plot_load => [1.0, 2.0],
        :plot_setup => [3.0, 4.0],
        :plot_data => [5.0, 6.0],
        :plot_draw => [7.0, 8.0, 9.0],
    )

    @test PERFORMANCE_BROWSER._make_timings_figure(live) !== nothing
    PERFORMANCE_BROWSER._update_live_timings!(live, timings)
    @test live.load_x[] == Float32[1, 2]
    @test live.total_x[] == Float32[1, 2, 3]
    @test length(live.load_x[]) == length(live.load_obs[])
    @test length(live.total_x[]) == length(live.total_obs[])

    @test PERFORMANCE_BROWSER._make_build_figure(live) !== nothing
    PERFORMANCE_BROWSER._sample_build_progress!(live, nothing)
    @test length(live.elapsed_obs[]) == 1
    @test length(live.elapsed_obs[]) == length(live.throughput_obs[])
    @test length(live.elapsed_obs[]) == length(live.rss_obs[])

    profiler = MeasurementBrowser.Profiling.ProfileSession(true, false, nothing, nothing)
    called = Ref(false)
    PERFORMANCE_BROWSER._profile_action!(profiler) do
        called[] = true
    end
    @test called[]
    @test profiler.state === :idle
    MeasurementBrowser.Profiling.close!(profiler)
end
