using MeasurementBrowser
using Test
using Dates

@testset "plot job controller guards" begin
    ui = Dict{Symbol,Any}()
    MeasurementBrowser._init_plot_state!(ui)
    @test ui[:plot_state] == :idle

    # Stale events must be ignored.
    ui[:active_plot_id] = 2
    ui[:plot_events] = Channel{NamedTuple}(8)
    ui[:active_plot_request] = Dict{Symbol,Any}(
        :kind => :single_file,
        :target => :main,
        :plot_key => :fresh_key,
        :project => MeasurementBrowser.RUO2_PROJECT,
        :debug_plot_mode => false,
        :measurement_kind => :iv,
        :device_params => Dict{Symbol,Any}(),
    )
    put!(ui[:plot_events], (kind=:error, plot_id=1, error=ErrorException("stale"), bt=nothing))
    put!(ui[:plot_events], (kind=:analyzed, plot_id=2, analyzed=nothing))
    MeasurementBrowser._poll_plot_events!(ui)
    @test ui[:plot_state] == :done
    @test ui[:plot_error] == ""
    @test ui[:active_plot_request] === nothing

    # Request while active should queue and move to canceling state.
    ui2 = Dict{Symbol,Any}()
    MeasurementBrowser._init_plot_state!(ui2)
    ui2[:plot_state] = :loading
    ui2[:active_plot_request] = Dict{Symbol,Any}(:job_key => :active)
    MeasurementBrowser._queue_plot_job!(ui2, Dict{Symbol,Any}(:job_key => :pending))
    @test ui2[:plot_state] == :canceling
    @test ui2[:pending_plot_request][:job_key] == :pending

    # Duplicate pending request should not change queued key.
    MeasurementBrowser._queue_plot_job!(ui2, Dict{Symbol,Any}(:job_key => :pending))
    @test ui2[:pending_plot_request][:job_key] == :pending

    # Main target application should clear cached figure when fig is nothing.
    ui3 = Dict{Symbol,Any}(:plot_figure => :old, :_last_plot_key => :old_key)
    MeasurementBrowser._apply_plot_result!(
        ui3,
        Dict{Symbol,Any}(:target => :main, :plot_key => :new_key),
        nothing,
    )
    @test !haskey(ui3, :plot_figure)
    @test !haskey(ui3, :_last_plot_key)

    # Extra-window target application should update only matching entry.
    cache_version = ("cache-id", "/tmp/cache.h5", 1, 0, 0, 0, 0)
    entry_a = Dict{Symbol,Any}(:target_id => "A", :figure => :figure_a)
    entry_b = Dict{Symbol,Any}(:target_id => "B", :figure => :figure_b)
    ui4 = Dict{Symbol,Any}(:open_plot_windows => [entry_a, entry_b])
    MeasurementBrowser._apply_plot_result!(
        ui4,
        Dict{Symbol,Any}(:target => :extra, :target_id => "A", :cache_version => cache_version),
        nothing,
    )
    @test !haskey(entry_a, :figure)
    @test entry_a[:cache_version] == cache_version
    @test entry_b[:figure] == :figure_b
end

@testset "extra plot windows preserve logical measurement identity" begin
    fixture = joinpath(
        @__DIR__,
        "3V PUND Fatigue [RuO2test_A9_VI_D1(2) ; 2025-10-01 17_12_33].csv",
    )
    meas = MeasurementInfo(fixture, MeasurementBrowser.RUO2_PROJECT)
    expanded = expand_measurement(MeasurementBrowser.RUO2_PROJECT, meas)
    cycle_two = only(filter(m -> m.parameters[:fatigue_cycle] == 2, expanded))

    entry = MeasurementBrowser._measurement_plot_window_entry(cycle_two)
    @test entry[:target_id] == cycle_two.id
    @test entry[:measurement_kind] == :pund
    @test entry[:params][:fatigue_cycle] == 2
    @test entry[:params][:voltage_V] == 3.0

    cache_id = "20260430_120007"
    ui = Dict{Symbol,Any}(
        :cache_identity => MeasurementBrowser.project_cache_identity(
            cache_id,
            MeasurementBrowser.RUO2_PROJECT,
            dirname(fixture),
        ),
        :cache_status => MeasurementBrowser.ProjectCacheStatus(1, 1, 1, 0, 0, 0, 0),
    )
    request = MeasurementBrowser._extra_plot_window_request(
        ui,
        MeasurementBrowser.RUO2_PROJECT,
        entry,
    )

    @test request[:target] == :extra
    @test request[:target_id] == cycle_two.id
    @test request[:measurement] == cycle_two
    @test request[:cache_identity].cache_id == cache_id
    @test request[:measurement_kind] == :pund
    @test request[:device_params][:fatigue_cycle] == 2
    @test request[:device_params][:voltage_V] == 3.0
end
