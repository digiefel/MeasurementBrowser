using Test, DataBrowserProfiling
using DataBrowserAPI: @timed_dbg
using TimerOutputs: ncalls

@testset "concurrent sections merge once and snapshots remain independent" begin
    reset_debug_timings!()
    op() = @timed_dbg "outer" begin
        @timed_dbg "inner" sum(1:1000)
    end
    foreach(fetch, [Threads.@spawn op() for _ in 1:4])
    snapshot = snapshot_debug_timings()
    @test ncalls(snapshot["outer"]) == ncalls(snapshot["outer"]["inner"]) == 4
    fetch(Threads.@spawn op())
    @test ncalls(snapshot["outer"]) == 4
    @test ncalls(snapshot_debug_timings()["outer"]) == 5
end

@testset "take resets recorded work and finish stops recording" begin
    reset_debug_timings!()
    fetch(Threads.@spawn (@timed_dbg "sample" sum(1:10)))
    @test haskey(take_debug_timings!(), "sample")
    @test !haskey(snapshot_debug_timings(), "sample")
    finish_debug_timings!()
    fetch(Threads.@spawn (@timed_dbg "after" sum(1:10)))
    @test !haskey(snapshot_debug_timings(), "after")
    reset_debug_timings!()
end

@testset "exceptions still publish their timing section" begin
    reset_debug_timings!()
    function fail()
        try
            @timed_dbg "failure" error("fixture failed")
        catch
        end
    end
    fetch(Threads.@spawn fail())
    @test haskey(snapshot_debug_timings(), "failure")
end
