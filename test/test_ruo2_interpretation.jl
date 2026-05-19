using Dates
using Test

using MeasurementBrowser

const _RUO2_FIXTURE_DIR = joinpath(@__DIR__, "fixtures", "RuO2")
const _RUO2_BREAKDOWN_FIXTURE =
    "Break of oxide 15V [RuO2test_A9_VI_D1D2(1) ; 2025-10-01 16_27_53].csv"
const _RUO2_PUND_FIXTURE =
    "3V FE PUND [RuO2test_A9_VI_D1(2) ; 2025-10-01 17_12_33].csv"
const _RUO2_TLM_FIXTURE =
    "TLM_4P [RuO2test_A9_VI_TLML100W2(12) ; 2025-10-01 16_21_45].csv"
const _RUO2_IV_FIXTURE =
    "RuO2test_A11_XI_FeCapBD_A1A2_20260509_184021_IVSweep.csv"
const _RUO2_CV_FIXTURE =
    "RuO2test_A9_VI_FeCap_A4_20260320_181744_298K_CVSweep.csv"
const _RUO2_WAKEUP_FIXTURE =
    "RuO2test_A11_XI_FeCap_A9_20260509_231259_PUND_Wakeup.csv"
const _RUO2_FATIGUE_FIXTURE =
    "RuO2test_A12_XI_FeCap_D6_20260511_222203_PUND_Fatigue_V2.csv"

_ruo2_fixture_path(name::AbstractString) = joinpath(_RUO2_FIXTURE_DIR, name)

@testset "RuO2 fixture interpretation contracts" begin
    @testset "single-file parsing" begin
        cases = [
            (
                file=_RUO2_PUND_FIXTURE,
                kind=:pund,
                location=["RuO2test", "A9", "VI", "D1"],
                timestamp=DateTime(2025, 10, 1, 17, 12, 33),
                parameters=Dict(:count => 2),
            ),
            (
                file=_RUO2_TLM_FIXTURE,
                kind=:tlm4p,
                location=["RuO2test", "A9", "VI", "TLML100W2"],
                timestamp=DateTime(2025, 10, 1, 16, 21, 45),
                parameters=Dict(:count => 12),
            ),
            (
                file=_RUO2_IV_FIXTURE,
                kind=:iv,
                location=["RuO2test_A11", "XI", "FeCapBD", "A1A2"],
                timestamp=DateTime(2026, 5, 9, 18, 40, 21),
                parameters=Dict{Symbol,Any}(),
            ),
            (
                file=_RUO2_CV_FIXTURE,
                kind=:cvsweep,
                location=["RuO2test_A9", "VI", "FeCap", "A4"],
                timestamp=DateTime(2026, 3, 20, 18, 17, 44),
                parameters=Dict(:temperature_K => 298),
            ),
        ]

        for case in cases
            path = _ruo2_fixture_path(case.file)
            source = index_source_file(path)
            measurement = MeasurementInfo(path, RUO2_PROJECT)
            @test source.id == path
            @test measurement.id == path
            @test measurement.filepath == path
            @test measurement.measurement_kind === case.kind
            @test measurement.device_info.location == case.location
            @test measurement.timestamp == case.timestamp
            for (key, value) in case.parameters
                @test measurement.parameters[key] == value
            end
        end
    end

    @testset "expanded measurement IDs and source paths" begin
        breakdown = measurements_for_file(RUO2_PROJECT, _ruo2_fixture_path(_RUO2_BREAKDOWN_FIXTURE))
        @test length(breakdown) == 2
        @test Set(last(m.device_info.location) for m in breakdown) == Set(["D1", "D2"])
        @test all(m -> m.measurement_kind === :breakdown, breakdown)
        @test all(m -> m.filepath == _ruo2_fixture_path(_RUO2_BREAKDOWN_FIXTURE), breakdown)
        @test Set(m.id for m in breakdown) == Set(
            _ruo2_fixture_path(_RUO2_BREAKDOWN_FIXTURE) .* ["#split=D1", "#split=D2"],
        )

        fatigue = measurements_for_file(RUO2_PROJECT, _ruo2_fixture_path(_RUO2_FATIGUE_FIXTURE))
        cycles, _ = MeasurementBrowser._ruo2_scan_fatigue_file(
            _ruo2_fixture_path(_RUO2_FATIGUE_FIXTURE),
        )
        @test length(fatigue) == length(cycles)
        @test all(m -> m.measurement_kind === :pund, fatigue)
        @test Set(m.id for m in fatigue) == Set(
            [_ruo2_fixture_path(_RUO2_FATIGUE_FIXTURE) * "#cycle=$cycle" for cycle in cycles],
        )

        wakeup = measurements_for_file(RUO2_PROJECT, _ruo2_fixture_path(_RUO2_WAKEUP_FIXTURE))
        @test Set(m.measurement_kind for m in wakeup) == Set([:wakeup_pn, :wakeup_pund])
        @test length(unique(m.id for m in wakeup)) == length(wakeup)
        @test all(m -> startswith(m.id, _ruo2_fixture_path(_RUO2_WAKEUP_FIXTURE) * "#split="), wakeup)
        @test all(m -> m.filepath == _ruo2_fixture_path(_RUO2_WAKEUP_FIXTURE), wakeup)

        pund = only(measurements_for_file(RUO2_PROJECT, _ruo2_fixture_path(_RUO2_PUND_FIXTURE)))
        @test isequal((
            pund.parameters[:wakeup_count],
            pund.parameters[:fatigue_count],
            pund.parameters[:wakeup_f],
            pund.parameters[:wakeup_V],
            pund.parameters[:fatigue_f],
            pund.parameters[:fatigue_V],
        ), (0, 0, NaN, NaN, NaN, NaN))
        @test (
            pund.stats[:V_base],
            pund.stats[:V_min],
            pund.stats[:V_max],
            pund.stats[:V_amp],
        ) == (0.0, -3.0, 3.0, 3.0)
    end

    @testset "scan over copied fixtures" begin
        fixture_names = [
            _RUO2_BREAKDOWN_FIXTURE,
            _RUO2_PUND_FIXTURE,
            _RUO2_TLM_FIXTURE,
            _RUO2_IV_FIXTURE,
            _RUO2_CV_FIXTURE,
            _RUO2_WAKEUP_FIXTURE,
        ]

        mktempdir() do dir
            for name in fixture_names
                cp(_ruo2_fixture_path(name), joinpath(dir, name))
            end

            source = scan_source(dir; project=RUO2_PROJECT)
            measurements = source.hierarchy.all_measurements
            @test length(source.files) == length(fixture_names)
            @test source.hierarchy.skipped_count == 0
            @test all(m -> startswith(m.filepath, dir), measurements)
            @test count(m -> m.measurement_kind === :breakdown, measurements) == 2
            @test count(m -> m.measurement_kind === :pund, measurements) == 1
            @test count(m -> m.measurement_kind === :tlm4p, measurements) == 1
            @test count(m -> m.measurement_kind === :iv, measurements) == 1
            @test count(m -> m.measurement_kind === :cvsweep, measurements) == 1
            @test count(m -> m.measurement_kind === :wakeup_pn, measurements) == 3
            @test count(m -> m.measurement_kind === :wakeup_pund, measurements) == 3
            @test haskey(source.hierarchy.index, ("RuO2test", "A9", "VI", "D1"))
            @test haskey(source.hierarchy.index, ("RuO2test", "A9", "VI", "D2"))
            @test haskey(source.hierarchy.index, ("RuO2test_A11", "XI", "FeCap", "A9"))
        end
    end
end
