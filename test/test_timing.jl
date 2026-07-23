@testset "always-on @timed main-task timer" begin
    B = DataBrowserGUI.Browser
    TO = B.TimerOutputs

    inner_sum(n) = sum(1:n)
    run_outer(n) = B.@timed "outer" begin
        B.@timed inner_sum(n)   # no label → derived from the callee ("inner_sum")
    end

    TO.reset_timer!(B.MAIN_TIMER)
    @test isempty(B.MAIN_TIMER.inner_timers)

    run_outer(5)
    run_outer(5)

    outer = B.MAIN_TIMER.inner_timers["outer"]
    @test TO.ncalls(outer) == 2                       # call count accumulates
    @test TO.time(outer) > 0                           # time recorded
    @test haskey(outer.inner_timers, "inner_sum")      # nested + label derived from callee
    @test TO.ncalls(outer.inner_timers["inner_sum"]) == 2

    TO.reset_timer!(B.MAIN_TIMER)
    @test isempty(B.MAIN_TIMER.inner_timers)           # reset clears the tree
end
