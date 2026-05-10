using MeasurementBrowser
using Test

function _copy_fixture(dst::String, filename::String)
    src = joinpath(@__DIR__, filename)
    cp(src, joinpath(dst, filename); force=true)
end

function _write_stub_csv(path::String)
    mkpath(dirname(path))
    write(path, "Time_s,Current_A,Voltage_V\n0.0,0.0,0.0\n")
end

@testset "scan_source progress and cancel" begin
    mktempdir() do dir
        _copy_fixture(dir, "Wakeup 3V [RuO2test_A9_VI_D1(2) ; 2025-10-01 17_10_48].csv")
        _copy_fixture(dir, "3V FE PUND [RuO2test_A9_VI_D1(2) ; 2025-10-01 17_12_33].csv")
        write(joinpath(dir, "not_a_measurement.csv"), "header\n1,2,3\n")

        source = scan_source(dir)
        @test source isa SourceScan
        @test source.hierarchy isa MeasurementHierarchy
        @test source.hierarchy.skipped_count == 1

        # Progress callback shape
        events = NamedTuple[]
        scan_source(dir; on_progress=(p) -> push!(events, p), count_first=true)
        @test !isempty(events)
        @test any(e -> e.phase == :counting, events)
        @test any(e -> e.phase == :scanning, events)
        last_scanning = last(filter(e -> e.phase == :scanning, events))
        @test last_scanning.processed_csv == 3
        @test last_scanning.total_csv == 3

        # Cooperative cancellation
        fired = Ref(false)
        @test_throws MeasurementBrowser.ScanCancelled scan_source(
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

@testset "scan_source device metadata matches exact path units across RuO2 layouts" begin
    mktempdir() do dir
        legacy_dir = joinpath(dir, "RuO2test_A9")
        mkpath(legacy_dir)
        legacy_file = "3V FE PUND [RuO2test_A9_VI_D1(2) ; 2025-10-01 17_12_33].csv"
        _copy_fixture(legacy_dir, legacy_file)

        nested_pund = joinpath(
            dir,
            "RuO2test_A10",
            "VI",
            "25um",
            "D1",
            "RuO2test_A10_VI_25um_D1_20260126_194830_273K_PUND.csv",
        )
        _write_stub_csv(nested_pund)

        nested_tlm = joinpath(
            dir,
            "RuO2test_A10",
            "VI",
            "TLM",
            "L100W2",
            "RuO2test_A10_VI_TLM_L100W2_20251128_012549_FourTerminalIV.csv",
        )
        _write_stub_csv(nested_tlm)

        write(joinpath(dir, "device_info.txt"), """
device,radius_um,tag,should_not_match
D1,50.0,leaf_device,
A9,,legacy_site,
A10,,substring_bug,bad
RuO2test_A10,,new_chip,
RuO2test/A9/VI/D1,10.0,legacy_exact,
RuO2test_A10/VI/25um/D1,25.0,new_exact,
""")

        source = scan_source(dir; project=RUO2_PROJECT)
        hierarchy = source.hierarchy
        @test length(hierarchy.all_measurements) == 3
        by_location = Dict(Tuple(m.device_info.location) => m for m in hierarchy.all_measurements)

        legacy = by_location[("RuO2test", "A9", "VI", "D1")]
        @test legacy.device_info.parameters[:radius_um] == 10.0
        @test legacy.device_info.parameters[:tag] == "legacy_exact"

        nested = by_location[("RuO2test_A10", "VI", "25um", "D1")]
        @test nested.device_info.parameters[:radius_um] == 25.0
        @test nested.device_info.parameters[:tag] == "new_exact"
        @test !haskey(nested.device_info.parameters, :should_not_match)

        tlm = by_location[("RuO2test_A10", "VI", "TLM", "L100W2")]
        @test tlm.measurement_kind == :tlm4p

        prefixed_device = joinpath(
            dir,
            "RuO2test_A10",
            "VI",
            "FecapBD",
            "RuO2test_A1A2",
            "RuO2test_A10_VI_FecapBD_A1A2_20251128_175954_OxideBreakdown.csv",
        )
        _write_stub_csv(prefixed_device)

        prefixed_info = parse_device_info(RUO2_PROJECT, index_source_file(prefixed_device))
        @test prefixed_info.location == ["RuO2test_A10", "VI", "FecapBD", "A1A2"]
    end
end
