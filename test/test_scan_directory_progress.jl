using MeasurementBrowser
using MeasurementBrowser: ItemRecord
using Test

@testset "scan_source progress and cancel" begin
    mktempdir() do dir
        write_test_source(joinpath(dir, "first.csv"))
        write_test_source(joinpath(dir, "second.csv"), 10)

        source = scan_test_source(TEST_PROJECT, dir)
        @test source isa MeasurementBrowser.SourceScan
        @test source.hierarchy isa MeasurementBrowser.Hierarchy
        @test source.hierarchy.skipped_count == 0

        events = NamedTuple[]
        MeasurementBrowser.scan_source(
            test_source(TEST_PROJECT, dir);
            on_progress=progress -> push!(events, progress),
            count_first=true,
        )
        @test any(event -> event.phase == :counting, events)
        @test any(event -> event.phase == :discovering, events)
        @test any(event -> event.phase == :scanning, events)
        @test any(event -> event.phase == :analyzing, events)
        scanning = filter(event -> event.phase == :scanning, events)
        @test maximum(event -> event.processed_source_items, scanning) == 2
        @test all(event -> event.total_source_items == 2, scanning)

        write(joinpath(dir, ".hidden.csv"), "")
        source = scan_test_source(TEST_PROJECT, dir)
        @test length(source.source_item_fingerprints) == 2

        streamed = Vector{ItemRecord}[]
        source = MeasurementBrowser.scan_source(
            test_source(TEST_PROJECT, dir);
            on_items=items -> push!(streamed, copy(items)),
        )
        @test sum(length, streamed) == length(source.hierarchy.all_items)
        @test Set(
            MeasurementBrowser.item_record_key(item)
            for batch in streamed
            for item in batch
        ) == Set(
            MeasurementBrowser.item_record_key(item)
            for item in source.hierarchy.all_items
        )

        write_test_source(joinpath(dir, "broken.csv"))
        source = scan_test_source(TEST_PROJECT, dir)
        @test length(source.source_item_fingerprints) == 3
        @test length(source.analysis_failures) == 1
        @test only(source.analysis_failures).source_item_id == joinpath(dir, "broken.csv")

        fired = Ref(false)
        @test_throws MeasurementBrowser.JobCancelled MeasurementBrowser.with_cancel(
            () -> begin
                fired[] && return true
                fired[] = true
                return false
            end,
        ) do
            MeasurementBrowser.scan_source(
                test_source(TEST_PROJECT, dir);
                count_first=true,
            )
        end
    end
end
