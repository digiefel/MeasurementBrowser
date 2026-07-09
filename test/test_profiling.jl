using JSON
using DataBrowser
import DataBrowserProfiling as Profiling
using CSV
using DataFrames: DataFrame, nrow
using Test

const PROFILE = Profiling

profile_snapshot() = (
    scan_done=0,
    scan_total=0,
    processing_done=0,
    processing_total=0,
    queue_depth=0,
)

@testset "structured profiler" begin
    @test PROFILE.environment_flag("MB_TEST_PROFILE_MISSING", true)
    withenv("MB_TEST_PROFILE_FLAG" => "1") do
        @test PROFILE.environment_flag("MB_TEST_PROFILE_FLAG")
    end
    withenv("MB_TEST_PROFILE_FLAG" => "0") do
        @test !PROFILE.environment_flag("MB_TEST_PROFILE_FLAG", true)
    end
    withenv("MB_TEST_PROFILE_FLAG" => "yes") do
        @test_throws ArgumentError PROFILE.environment_flag("MB_TEST_PROFILE_FLAG")
    end
    withenv("MB_TEST_PROFILE_PATH" => "") do
        @test_throws ArgumentError PROFILE.environment_path("MB_TEST_PROFILE_PATH")
    end
    @test_throws ArgumentError PROFILE.ProfileSession(false, true, nothing, nothing)
    @test_throws ArgumentError PROFILE.ProfileSession(false, false, "/tmp/out.json", nothing)
    @test PROFILE.process_rss_bytes() > 0

    for session in (
        PROFILE.ProfileSession(false, false, nothing, nothing),
        PROFILE.ProfileSession(true, false, nothing, nothing),
    )
        attributes_evaluated = Ref(false)
        result = PROFILE.@profile_span session :test :disabled begin
            attributes_evaluated[] = true
            PROFILE.ProfileAttributes()
        end begin
            42
        end
        @test result == 42
        @test !attributes_evaluated[]
        PROFILE.close!(session)
    end

    mktempdir() do dir
        session = PROFILE.ProfileSession(true, false, nothing, nothing)
        PROFILE.start!(session, profile_snapshot)
        PROFILE.@profile_span session :test :outer PROFILE.ProfileAttributes(items=2) begin
            task = Threads.@spawn PROFILE.@profile_span session :test :inner PROFILE.ProfileAttributes(
                batch_size=2,
            ) begin
                sleep(0.001)
            end
            wait(task)
        end
        concurrent = [Threads.@spawn begin
            PROFILE.@profile_span session :concurrent :work PROFILE.ProfileAttributes(
                item_id=string(index),
            ) begin
                yield()
            end
        end for index in 1:128]
        foreach(wait, concurrent)
        @test_throws ErrorException begin
            PROFILE.@profile_span session :test :failure PROFILE.ProfileAttributes() begin
                error("profiled failure")
            end
        end
        report = PROFILE.stop!(session)
        @test length(report.events) == 131
        @test PROFILE.event_count(session) == 131
        @test isempty(session.events)
        @test count(event -> event.category === :concurrent, report.events) == 128
        @test length(Set(event.id for event in report.events)) == length(report.events)
        outer = only(event for event in report.events if event.operation === :outer)
        inner = only(event for event in report.events if event.operation === :inner)
        failed = only(event for event in report.events if event.operation === :failure)
        @test inner.parent_id == outer.id
        @test failed.status === :error
        @test report.dropped_events == 0
        @test !isempty(report.summary)

        output = joinpath(dir, "profile.json")
        @test PROFILE.export!(session, output; numbered=true) == output
        second = PROFILE.export!(session, output; numbered=true)
        @test second == joinpath(dir, "profile-2.json")
        parsed = JSON.parsefile(output)
        @test parsed["displayTimeUnit"] == "ms"
        @test any(event -> event["ph"] == "X", parsed["traceEvents"])
        @test any(event -> event["ph"] == "C", parsed["traceEvents"])

        PROFILE.reset!(session)
        PROFILE.start!(session, profile_snapshot)
        resize!(session.events, 2)
        for index in 1:3
            PROFILE.@profile_span session :overflow Symbol("event_$index") PROFILE.ProfileAttributes() begin
                nothing
            end
        end
        overflow = PROFILE.stop!(session)
        @test length(overflow.events) == 2
        @test overflow.dropped_events == 1

        PROFILE.reset!(session)
        @test_throws ErrorException PROFILE.start!(session, () -> error("snapshot failed"))
        @test session.state === :error
        PROFILE.reset!(session)
        @test session.state === :idle
        PROFILE.close!(session)
    end
end

@testset "profile keyword precedence" begin
    mktempdir() do dir
        project = DataBrowser.define_project("ProfileFlags_$(basename(dir))")
        source = DataBrowser.DirectorySource(dir)
        withenv(
            "MB_PROFILE_INTERNAL" => "1",
            "MB_PROFILE_CPU" => "1",
            "MB_PROFILE_OUTPUT" => joinpath(dir, "ignored.json"),
            "MB_CRASH_TRACE" => nothing,
        ) do
            workspace = DataBrowser.open_workspace(
                project,
                source;
                profile_internal=false,
                profile_cpu=false,
                profile_output=nothing,
                crash_trace=nothing,
            )
            try
                @test workspace.profiler.state === :disabled
            finally
                DataBrowser.close_workspace!(workspace)
            end
        end
    end
end

@testset "crash trace survives process termination" begin
    mktempdir() do dir
        path = joinpath(dir, "crash.jsonl")
        project = dirname(@__DIR__)
        code = """
            import DataBrowserProfiling as Profiling
            session = Profiling.ProfileSession(false, false, nothing, $(repr(path)))
            Profiling.start_span!(session, :test, :unfinished)
            exit(17)
        """
        command = `$(Base.julia_cmd()) --project=$project -e $code`
        process = run(ignorestatus(command))
        @test process.exitcode == 17
        records = JSON.parse.(readlines(path))
        @test first(records)["type"] == "session_start"
        starts = [record for record in records if record["type"] == "span_start"]
        ends = [record for record in records if record["type"] == "span_end"]
        @test length(starts) == 1
        @test isempty(ends)
    end
end
