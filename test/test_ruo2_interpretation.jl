using Dates
using Test

using MeasurementBrowser

const RUO2_PROJECT = MeasurementBrowser.RUO2_PROJECT

const _RUO2_FIXTURE_DIR = joinpath(@__DIR__, "fixtures", "RuO2")
const _RUO2_BREAKDOWN_FIXTURE =
    "Break of oxide 15V [RuO2test_A9_VI_D1D2(1) ; 2025-10-01 16_27_53].csv"
const _RUO2_PUND_FIXTURE =
    "3V FE PUND [RuO2test_A9_VI_D1(2) ; 2025-10-01 17_12_33].csv"
const _RUO2_TLM_FIXTURE =
    "TLM_4P [RuO2test_A9_VI_TLML100W2(12) ; 2025-10-01 16_21_45].csv"
const _RUO2_IV_FIXTURE =
    "RuO2test_A11_XI_FeCapBD_A1A2_20260509_184021_IVSweep.csv"
const _RUO2_TLM_DEVICE_IV_FIXTURE =
    "RuO2test_A11_XI_TLM_L100W4_20260609_183048_14K_IVSweep.csv"
const _RUO2_FECAP_BD_FIXTURE =
    "RuO2test_A10_VI_FeCap_BD_A1A2_20251128_175954_OxideBreakdown.csv"
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
            ),
            (
                file=_RUO2_TLM_FIXTURE,
                kind=:tlm4p,
                location=["RuO2test", "A9", "VI", "TLML100W2"],
                timestamp=DateTime(2025, 10, 1, 16, 21, 45),
            ),
            (
                file=_RUO2_IV_FIXTURE,
                kind=:iv,
                location=["RuO2test_A11", "XI", "FeCapBD", "A1A2"],
                timestamp=DateTime(2026, 5, 9, 18, 40, 21),
            ),
            (
                file=_RUO2_TLM_DEVICE_IV_FIXTURE,
                kind=:iv,
                location=["RuO2test_A11", "XI", "TLM", "L100W4"],
                timestamp=DateTime(2026, 6, 9, 18, 30, 48),
            ),
            (
                file=_RUO2_FECAP_BD_FIXTURE,
                kind=:breakdown,
                location=["RuO2test_A10", "VI", "FeCapBD", "A1A2"],
                timestamp=DateTime(2025, 11, 28, 17, 59, 54),
            ),
            (
                file=_RUO2_CV_FIXTURE,
                kind=:cvsweep,
                location=["RuO2test_A9", "VI", "FeCap", "A4"],
                timestamp=DateTime(2026, 3, 20, 18, 17, 44),
            ),
        ]

        for case in cases
            path = _ruo2_fixture_path(case.file)
            source = MeasurementBrowser.index_source_file(path)
            @test source.unique_id == path
            @test detect_kind(RUO2_PROJECT, basename(path)) === case.kind
            @test parse_device_info(RUO2_PROJECT, source).location == case.location
            @test source.timestamp == case.timestamp
        end

        @test MeasurementBrowser._ruo2_location_from_filename(
            "RuO2testA10_VI_FeCap_BD_A1A2_20251128_175954_OxideBreakdown.csv",
        ) == ["RuO2testA10", "VI", "FeCapBD", "A1A2"]
    end

    @testset "partial measurement labels" begin
        path = _ruo2_fixture_path(_RUO2_PUND_FIXTURE)
        measurement = only(MeasurementBrowser.interpret_file(RUO2_PROJECT, MeasurementBrowser.index_source_file(path)))
        @test isempty(measurement.stats)
        @test display_label(RUO2_PROJECT, measurement) isa String
    end

    @testset "analysis failure stays local" begin
        measurement = MeasurementInfo(;
            filepath=joinpath(tempdir(), "missing_pund.csv"),
            measurement_kind=:pund,
            device_info=DeviceInfo(["RuO2test", "A0", "I", "D0"]),
            clean_title="Missing PUND",
        )
        failures = MeasurementBrowser.compute_and_add_measurement_stats!(
            RUO2_PROJECT,
            [measurement],
            SourceFile[],
        )
        @test length(failures) == 1
        @test failures[1].measurement_id == measurement.unique_id
        @test measurement.stats[:wakeup_count] == 0
        @test measurement.stats[:fatigue_count] == 0
    end

    @testset "expanded measurement IDs and source paths" begin
        breakdown = measurements_for_file(RUO2_PROJECT, _ruo2_fixture_path(_RUO2_BREAKDOWN_FIXTURE))
        @test length(breakdown) == 2
        @test Set(last(m.device_info.location) for m in breakdown) == Set(["D1", "D2"])
        @test all(m -> m.measurement_kind === :breakdown, breakdown)
        @test all(m -> m.filepath == _ruo2_fixture_path(_RUO2_BREAKDOWN_FIXTURE), breakdown)
        @test Set(m.unique_id for m in breakdown) == Set(
            _ruo2_fixture_path(_RUO2_BREAKDOWN_FIXTURE) .* ["#device=D1", "#device=D2"],
        )

        fatigue = measurements_for_file(RUO2_PROJECT, _ruo2_fixture_path(_RUO2_FATIGUE_FIXTURE))
        cycles = unique(MeasurementBrowser.read_pund_fatigue_file(
            _ruo2_fixture_path(_RUO2_FATIGUE_FIXTURE),
        ).cycle)
        @test length(fatigue) == length(cycles)
        @test all(m -> m.measurement_kind === :pund, fatigue)
        @test Set(m.unique_id for m in fatigue) == Set(
            [_ruo2_fixture_path(_RUO2_FATIGUE_FIXTURE) * "#fatigue_count=$cycle" for cycle in cycles],
        )

        wakeup = measurements_for_file(RUO2_PROJECT, _ruo2_fixture_path(_RUO2_WAKEUP_FIXTURE))
        @test Set(m.measurement_kind for m in wakeup) == Set([:wakeup_pn, :wakeup_pund])
        @test length(unique(m.unique_id for m in wakeup)) == length(wakeup)
        @test all(m -> startswith(m.unique_id, _ruo2_fixture_path(_RUO2_WAKEUP_FIXTURE) * "#wakeup_V="), wakeup)
        @test all(m -> m.filepath == _ruo2_fixture_path(_RUO2_WAKEUP_FIXTURE), wakeup)

        pund = only(measurements_for_file(RUO2_PROJECT, _ruo2_fixture_path(_RUO2_PUND_FIXTURE)))
        @test isequal((
            pund.stats[:wakeup_count],
            pund.stats[:fatigue_count],
        ), (0, 0))
        @test !haskey(pund.stats, :wakeup_f)
        @test !haskey(pund.stats, :wakeup_V)
        @test !haskey(pund.stats, :fatigue_f)
        @test !haskey(pund.stats, :fatigue_V)
        @test (
            pund.stats[:V_base],
            pund.stats[:V_min],
            pund.stats[:V_max],
            pund.stats[:V_amp],
        ) == (0.0, -3.0, 3.0, 3.0)

        cv = only(measurements_for_file(RUO2_PROJECT, _ruo2_fixture_path(_RUO2_CV_FIXTURE)))
        @test cv.parameters[:temperature_K] == 298.0
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

            source = MeasurementBrowser.scan_source(dir; project=RUO2_PROJECT)
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
