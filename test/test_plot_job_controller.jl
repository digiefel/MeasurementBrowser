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
        :target => :main,
        :plot_key => :fresh_key,
        :project => MeasurementBrowser.RUO2_PROJECT,
        :debug_plot_mode => false,
    )
    put!(ui[:plot_events], (kind=:error, plot_id=1, error=ErrorException("stale"), bt=nothing))
    put!(ui[:plot_events], (kind=:prepared, plot_id=2, payload=nothing))
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
    ts = DateTime(2026, 1, 1)
    entry_a = Dict{Symbol,Any}(:target_id => "A", :figure => :figure_a)
    entry_b = Dict{Symbol,Any}(:target_id => "B", :figure => :figure_b)
    ui4 = Dict{Symbol,Any}(:open_plot_windows => [entry_a, entry_b])
    MeasurementBrowser._apply_plot_result!(
        ui4,
        Dict{Symbol,Any}(:target => :extra, :target_id => "A", :mtime => ts),
        nothing,
    )
    @test !haskey(entry_a, :figure)
    @test entry_a[:mtime] == ts
    @test entry_b[:figure] == :figure_b
end
