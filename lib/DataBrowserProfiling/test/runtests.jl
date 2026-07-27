using Test
using DataBrowserProfiling
using DataBrowserProfiling: snapshot_debug_timings, take_debug_timings!,
    finish_debug_timings!, reset_debug_timings!, process_rss_bytes
using DataBrowserAPI: @timed_dbg
import DataBrowserAPI
import TimerOutputs
using TimerOutputs: TimerOutput

# Loading DataBrowserProfiling installs the timing hooks into DataBrowserAPI and
# raises its profiling level, so `@timed_dbg` sections now record.

@testset "DataBrowserProfiling" begin
    @testset "loading enables instrumentation" begin
        @test DataBrowserAPI.profile_level() == 1
    end

    @testset "nested sections nest, and each task merges once" begin
        reset_debug_timings!()
        op(n) = @timed_dbg "outer" begin
            @timed_dbg "inner" sum(1:n)
        end
        # Each spawned task owns its own TimerOutput and merges on outermost exit.
        foreach(fetch, [Threads.@spawn op(1000) for _ in 1:4])
        snap = snapshot_debug_timings()
        @test snap isa TimerOutput
        @test TimerOutputs.ncalls(snap["outer"]) == 4
        @test haskey(snap["outer"].inner_timers, "inner")
        @test TimerOutputs.ncalls(snap["outer"].inner_timers["inner"]) == 4
    end

    @testset "a snapshot is independent of later work" begin
        reset_debug_timings!()
        fetch(Threads.@spawn (@timed_dbg "a" sum(1:10)))
        snap1 = snapshot_debug_timings()
        n1 = TimerOutputs.ncalls(snap1["a"])
        fetch(Threads.@spawn (@timed_dbg "a" sum(1:10)))
        @test TimerOutputs.ncalls(snap1["a"]) == n1                       # unchanged
        @test TimerOutputs.ncalls(snapshot_debug_timings()["a"]) == n1 + 1
    end

    @testset "take! returns and resets the master" begin
        reset_debug_timings!()
        fetch(Threads.@spawn (@timed_dbg "x" sum(1:10)))
        taken = take_debug_timings!()
        @test haskey(taken.inner_timers, "x")
        @test isempty(snapshot_debug_timings().inner_timers)             # reset
    end

    @testset "an exception still closes and records the section" begin
        reset_debug_timings!()
        function boom()
            try
                @timed_dbg "boom" error("expected")
            catch
            end
        end
        fetch(Threads.@spawn boom())
        @test haskey(snapshot_debug_timings().inner_timers, "boom")
    end

    @testset "finish! stops recording" begin
        reset_debug_timings!()
        @test finish_debug_timings!() isa TimerOutput
        fetch(Threads.@spawn (@timed_dbg "after" sum(1:10)))
        @test isempty(snapshot_debug_timings().inner_timers)            # not recorded
        reset_debug_timings!()                                          # re-enable
    end

    @testset "process_rss_bytes reports a positive size" begin
        @test process_rss_bytes() > 0
    end
end
