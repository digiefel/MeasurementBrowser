# @timed and its internals live in the Browser submodule; imported at top level so the macro
# resolves at expansion time. TimerOutputs' binding is pulled through Browser.
using DataBrowserGUI.Browser: @timed, MAIN_TIMER, PerfHistory, PERF_HISTORY_CAP,
    _record_throughput!, TimerOutputs

@testset "always-on @timed main-task timer" begin
    inner_sum(n) = sum(1:n)
    run_outer(n) = @timed "outer" begin
        @timed inner_sum(n)   # no label → label is the expression text ("inner_sum(n)")
    end

    TimerOutputs.reset_timer!(MAIN_TIMER)
    @test isempty(MAIN_TIMER.inner_timers)

    run_outer(5)
    run_outer(5)

    outer = MAIN_TIMER.inner_timers["outer"]
    @test TimerOutputs.ncalls(outer) == 2                    # call count accumulates
    @test TimerOutputs.time(outer) > 0                        # time recorded
    @test haskey(outer.inner_timers, "inner_sum(n)")          # nested + label is the expression text
    @test TimerOutputs.ncalls(outer.inner_timers["inner_sum(n)"]) == 2

    TimerOutputs.reset_timer!(MAIN_TIMER)
    @test isempty(MAIN_TIMER.inner_timers)                    # reset clears the tree
end

@testset "throughput sampling (throttle + cap + rates)" begin
    h = PerfHistory()

    # First admitted call seeds baselines only — no data point, no spike.
    @test _record_throughput!(h, 100.0, 10, 2, 5, 3, 4) == true
    @test isempty(h.items_per_s)
    @test h.last_completed == 10

    # Too soon (dt < 0.25): throttled, nothing recorded.
    @test _record_throughput!(h, 100.1, 20, 2, 5, 3, 4) == false
    @test isempty(h.items_per_s)

    # dt = 1.0s: deltas become per-second rates; levels pass through.
    @test _record_throughput!(h, 101.0, 20, 3, 7, 8, 9) == true
    @test h.items_per_s[end] == 10.0f0    # (20-10)/1
    @test h.active[end] == 3.0f0
    @test h.pending_rows[end] == 7.0f0
    @test h.scan_per_s[end] == 5.0f0      # (8-3)/1
    @test h.cache_per_s[end] == 5.0f0     # (9-4)/1

    # Ring buffers stay capped at PERF_HISTORY_CAP.
    for i in 1:(PERF_HISTORY_CAP + 20)
        _record_throughput!(h, 101.0 + i, 20, 0, 0, 8, 9)
    end
    @test length(h.items_per_s) == PERF_HISTORY_CAP
end
