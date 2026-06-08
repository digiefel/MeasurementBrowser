using MeasurementBrowser
using Test

function _copy_fixture(dst::String, filename::String; subdir::String="RuO2")
    src = joinpath(@__DIR__, "fixtures", subdir, filename)
    cp(src, joinpath(dst, filename); force=true)
end

@testset "scan_source progress and cancel" begin
    mktempdir() do dir
        _copy_fixture(dir, "TLM_4P [RuO2test_A9_VI_TLML100W2(12) ; 2025-10-01 16_21_45].csv")
        _copy_fixture(dir, "3V FE PUND [RuO2test_A9_VI_D1(2) ; 2025-10-01 17_12_33].csv")

        source = scan_source(dir)
        @test source isa SourceScan
        @test source.hierarchy isa MeasurementHierarchy
        @test source.hierarchy.skipped_count == 0

        # Progress callback shape
        events = NamedTuple[]
        scan_source(dir; on_progress=(p) -> push!(events, p), count_first=true)
        @test !isempty(events)
        @test any(e -> e.phase == :counting, events)
        @test any(e -> e.phase == :discovering, events)
        @test any(e -> e.phase == :scanning, events)
        @test any(e -> e.phase == :analyzing, events)
        last_discovering = last(filter(e -> e.phase == :discovering, events))
        @test last_discovering.processed_csv == 2
        @test last_discovering.total_csv == 0
        scanning_events = filter(e -> e.phase == :scanning, events)
        @test maximum(e -> e.processed_csv, scanning_events) == 2
        @test all(e -> e.total_csv == 2, scanning_events)

        # Measurement batches are emitted before scan_source returns its final SourceScan.
        streamed = Vector{MeasurementInfo}[]
        source = scan_source(dir; on_measurements=(measurements) -> push!(streamed, copy(measurements)))
        @test !isempty(streamed)
        @test sum(length, streamed) == length(source.hierarchy.all_measurements)
        @test Set(m.unique_id for batch in streamed for m in batch) ==
              Set(m.unique_id for m in source.hierarchy.all_measurements)

        # Cooperative cancellation
        fired = Ref(false)
        @test_throws MeasurementBrowser.JobCancelled MeasurementBrowser._with_cancel(
            () -> begin
                if fired[]
                    return true
                end
                fired[] = true
                return false
            end,
        ) do
            scan_source(dir; count_first=true)
        end
    end
end
