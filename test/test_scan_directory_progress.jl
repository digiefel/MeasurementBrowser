using MeasurementBrowser
using MeasurementBrowser: ItemRecord
using Test

@testset "scan_source progress and cancel" begin
    mktempdir() do dir
        write_test_source(joinpath(dir, "first.csv"))
        write_test_source(joinpath(dir, "second.csv"), 10)

        source = MeasurementBrowser.scan_source(dir; project=TEST_PROJECT)
        @test source isa MeasurementBrowser.SourceScan
        @test source.hierarchy isa MeasurementBrowser.Hierarchy
        @test source.hierarchy.skipped_count == 0

        events = NamedTuple[]
        MeasurementBrowser.scan_source(
            dir;
            project=TEST_PROJECT,
            on_progress=progress -> push!(events, progress),
            count_first=true,
        )
        @test any(event -> event.phase == :counting, events)
        @test any(event -> event.phase == :discovering, events)
        @test any(event -> event.phase == :scanning, events)
        @test any(event -> event.phase == :analyzing, events)
        scanning = filter(event -> event.phase == :scanning, events)
        @test maximum(event -> event.processed_csv, scanning) == 2
        @test all(event -> event.total_csv == 2, scanning)

        write(joinpath(dir, ".hidden.csv"), "")
        source = MeasurementBrowser.scan_source(dir; project=TEST_PROJECT)
        @test length(source.files) == 2

        streamed = Vector{ItemRecord}[]
        source = MeasurementBrowser.scan_source(
            dir;
            project=TEST_PROJECT,
            on_measurements=measurements -> push!(streamed, copy(measurements)),
        )
        @test sum(length, streamed) == length(source.hierarchy.all_measurements)
        @test Set(
            measurement.unique_id
            for batch in streamed
            for measurement in batch
        ) == Set(
            measurement.unique_id
            for measurement in source.hierarchy.all_measurements
        )

        write_test_source(joinpath(dir, "broken.csv"))
        source = MeasurementBrowser.scan_source(dir; project=TEST_PROJECT)
        @test length(source.files) == 3
        @test length(source.analysis_failures) == 1
        @test only(source.analysis_failures).filepath == joinpath(dir, "broken.csv")

        fired = Ref(false)
        @test_throws MeasurementBrowser.JobCancelled MeasurementBrowser.with_cancel(
            () -> begin
                fired[] && return true
                fired[] = true
                return false
            end,
        ) do
            MeasurementBrowser.scan_source(
                dir;
                project=TEST_PROJECT,
                count_first=true,
            )
        end
    end
end
