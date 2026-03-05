using MeasurementBrowser
using Test

function _copy_fixture(dst::String, filename::String)
    src = joinpath(@__DIR__, filename)
    cp(src, joinpath(dst, filename); force=true)
end

@testset "scan_directory progress and cancel" begin
    mktempdir() do dir
        _copy_fixture(dir, "Wakeup 3V [RuO2test_A9_VI_D1(2) ; 2025-10-01 17_10_48].csv")
        _copy_fixture(dir, "3V FE PUND [RuO2test_A9_VI_D1(2) ; 2025-10-01 17_12_33].csv")
        write(joinpath(dir, "not_a_measurement.csv"), "header\n1,2,3\n")

        # Backwards-compatible call shape
        h = scan_directory(dir)
        @test h isa MeasurementHierarchy
        @test h.skipped_count >= 1

        # Progress callback shape
        events = NamedTuple[]
        scan_directory(dir; on_progress=(p) -> push!(events, p), count_first=true)
        @test !isempty(events)
        @test any(e -> e.phase == :counting, events)
        @test any(e -> e.phase == :scanning, events)
        last_scanning = last(filter(e -> e.phase == :scanning, events))
        @test last_scanning.processed_csv == 3
        @test last_scanning.total_csv == 3

        # Cooperative cancellation
        fired = Ref(false)
        @test_throws MeasurementBrowser.ScanCancelled scan_directory(
            dir;
            should_cancel=() -> begin
                if fired[]
                    return true
                end
                fired[] = true
                return false
            end,
            count_first=true,
        )
    end
end
