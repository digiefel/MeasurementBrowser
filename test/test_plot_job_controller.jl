using MeasurementBrowser
using Test
using Dates

function _test_plot_job(;
    job_key=:job,
    target=:main,
    target_id="main",
    plot_key=nothing,
    cache_version=nothing,
    debug=false,
)
    fixture = joinpath(@__DIR__, "fixtures", "RuO2", "3V FE PUND [RuO2test_A9_VI_D1(2) ; 2025-10-01 17_12_33].csv")
    measurement = MeasurementInfo(fixture, MeasurementBrowser.RUO2_PROJECT)
    params = merge(measurement.device_info.parameters, measurement.parameters)
    identity = MeasurementBrowser.project_cache_identity(
        "20260430_120007",
        MeasurementBrowser.RUO2_PROJECT,
        dirname(fixture),
    )
    return MeasurementBrowser.PlotJob(
        job_key,
        plot_key,
        MeasurementBrowser.RUO2_PROJECT,
        target,
        target_id,
        [measurement],
        measurement.measurement_kind,
        [params],
        identity,
        cache_version,
        debug,
    )
end

@testset "plot job controller guards" begin
    ui = Dict{Symbol,Any}()
    MeasurementBrowser._init_plot_state!(ui)
    @test ui[:plot_state] == :idle

    # Stale events must be ignored.
    ui[:active_plot_id] = 2
    ui[:plot_events] = Channel{NamedTuple}(8)
    ui[:active_plot_job] = _test_plot_job(plot_key=:fresh_key)
    put!(ui[:plot_events], (kind=:error, plot_id=1, error=ErrorException("stale"), bt=nothing))
    put!(ui[:plot_events], (kind=:data, plot_id=2, data=nothing))
    MeasurementBrowser._poll_plot_events!(ui)
    @test ui[:plot_state] == :error
    @test occursin("no figure", ui[:plot_error])
    @test !occursin("stale", ui[:plot_error])
    @test ui[:active_plot_job] === nothing

    # Job while active should queue and move to canceling state.
    ui2 = Dict{Symbol,Any}()
    MeasurementBrowser._init_plot_state!(ui2)
    ui2[:plot_state] = :loading
    ui2[:active_plot_job] = _test_plot_job(job_key=:active)
    MeasurementBrowser._queue_plot_job!(ui2, _test_plot_job(job_key=:pending))
    @test ui2[:plot_state] == :canceling
    @test ui2[:pending_plot_job].job_key == :pending

    # Duplicate pending job should not change queued key.
    MeasurementBrowser._queue_plot_job!(ui2, _test_plot_job(job_key=:pending))
    @test ui2[:pending_plot_job].job_key == :pending

    # Main target application should clear cached figure and remember failed key
    # when drawing returns no figure, preventing an immediate retry loop.
    ui3 = Dict{Symbol,Any}(:plot_figure => :old, :_last_plot_key => :old_key)
    failed = MeasurementBrowser._finish_plot_job!(
        ui3,
        _test_plot_job(plot_key=:new_key);
        fig=nothing,
    )
    @test failed
    @test !haskey(ui3, :plot_figure)
    @test ui3[:_last_plot_key] == :new_key
    @test occursin("no figure", ui3[:plot_error])

    # Extra-window target application should update only matching entry.
    cache_version = ("cache-id", "/tmp/cache.h5", 1, 0, 0, 0, 0)
    entry_a = Dict{Symbol,Any}(:target_id => "A", :figure => :figure_a)
    entry_b = Dict{Symbol,Any}(:target_id => "B", :figure => :figure_b)
    ui4 = Dict{Symbol,Any}(:open_plot_windows => [entry_a, entry_b])
    failed = MeasurementBrowser._finish_plot_job!(
        ui4,
        _test_plot_job(target=:extra, target_id="A", cache_version=cache_version);
        fig=nothing,
    )
    @test failed
    @test !haskey(entry_a, :figure)
    @test entry_a[:cache_version] == cache_version
    @test occursin("no figure", entry_a[:plot_error])
    @test entry_b[:figure] == :figure_b

    # Failed worker events should apply to the target instead of only logging.
    ui5 = Dict{Symbol,Any}(:plot_figure => :old, :_last_plot_key => :old_key)
    MeasurementBrowser._init_plot_state!(ui5)
    ui5[:plot_figure] = :old
    ui5[:_last_plot_key] = :old_key
    ui5[:active_plot_id] = 1
    ui5[:plot_events] = Channel{NamedTuple}(8)
    ui5[:active_plot_job] = _test_plot_job(plot_key=:failed_key)
    put!(ui5[:plot_events], (kind=:error, plot_id=1, error=ErrorException("boom"), bt=nothing))
    MeasurementBrowser._poll_plot_events!(ui5)
    @test ui5[:plot_state] == :error
    @test !haskey(ui5, :plot_figure)
    @test ui5[:_last_plot_key] == :failed_key
    @test occursin("boom", ui5[:plot_error])
    @test ui5[:active_plot_job] === nothing

    entry_c = Dict{Symbol,Any}(:target_id => "C", :figure => :figure_c)
    ui6 = Dict{Symbol,Any}(:open_plot_windows => [entry_c])
    MeasurementBrowser._init_plot_state!(ui6)
    ui6[:open_plot_windows] = [entry_c]
    ui6[:active_plot_id] = 1
    ui6[:plot_events] = Channel{NamedTuple}(8)
    ui6[:active_plot_job] = _test_plot_job(
        target=:extra,
        target_id="C",
        cache_version=cache_version,
    )
    put!(ui6[:plot_events], (kind=:error, plot_id=1, error=ErrorException("extra boom"), bt=nothing))
    MeasurementBrowser._poll_plot_events!(ui6)
    @test ui6[:plot_state] == :error
    @test !haskey(entry_c, :figure)
    @test entry_c[:cache_version] == cache_version
    @test occursin("extra boom", entry_c[:plot_error])
end

@testset "extra plot windows preserve logical measurement identity" begin
    fixture = joinpath(
        @__DIR__,
        "fixtures", "RuO2",
        "3V FE PUND [RuO2test_A9_VI_D1(2) ; 2025-10-01 17_12_33].csv",
    )
    expanded = measurements_for_file(MeasurementBrowser.RUO2_PROJECT, fixture)
    pund = only(expanded)

    entry = MeasurementBrowser._measurement_plot_window_entry(pund)
    @test entry[:target_id] == pund.id
    @test entry[:measurement_kind] == :pund

    cache_id = "20260430_120007"
    ui = Dict{Symbol,Any}(
        :cache_identity => MeasurementBrowser.project_cache_identity(
            cache_id,
            MeasurementBrowser.RUO2_PROJECT,
            dirname(fixture),
        ),
        :cache_status => MeasurementBrowser.ProjectCacheStatus(1, 1, 1, 0, 0, 0, 0),
    )
    job = MeasurementBrowser._plot_job(
        ui,
        MeasurementBrowser.RUO2_PROJECT,
        [pund],
        :pund,
        :extra;
        target_id=pund.id,
    )

    @test job.target == :extra
    @test job.target_id == pund.id
    @test only(job.measurements) == pund
    @test job.cache_identity.cache_id == cache_id
    @test job.plot_kind == :pund
end
