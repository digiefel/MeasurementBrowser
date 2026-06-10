using Test
using DataFrames: nrow
using MeasurementBrowser

const RUO2_PROJECT = MeasurementBrowser.RUO2_PROJECT

const _WAKEUP_FIXTURE = joinpath(
    @__DIR__,
    "fixtures", "RuO2",
    "RuO2test_A11_XI_FeCap_A9_20260509_231259_PUND_Wakeup.csv",
)
const _WAKEUP_ACCUM_FIRST = joinpath(
    @__DIR__,
    "fixtures", "RuO2",
    "RuO2test_A11_XI_FeCap_A2_20260509_223436_PUND_Wakeup.csv",
)
const _WAKEUP_ACCUM_SECOND = joinpath(
    @__DIR__,
    "fixtures", "RuO2",
    "RuO2test_A11_XI_FeCap_A2_20260509_223508_PUND_Wakeup.csv",
)

@testset "pund_wakeup detection and device parsing" begin
    @test isfile(_WAKEUP_FIXTURE)
    @test detect_kind(RUO2_PROJECT, basename(_WAKEUP_FIXTURE)) == :pund_wakeup

    source = MeasurementBrowser.index_source_file(_WAKEUP_FIXTURE)
    @test parse_device_info(RUO2_PROJECT, source).location == ["RuO2test_A11", "XI", "FeCap", "A9"]
end

@testset "pund_wakeup expansion" begin
    expanded = measurements_for_file(RUO2_PROJECT, _WAKEUP_FIXTURE)

    # 3 amplitudes × 2 segment kinds (pn + pund) = 6 items
    @test length(expanded) == 6

    pn_items   = filter(m -> m.measurement_kind === :wakeup_pn,   expanded)
    pund_items = filter(m -> m.measurement_kind === :wakeup_pund, expanded)
    @test length(pn_items)   == 3
    @test length(pund_items) == 3

    @test Set(m.stats[:wakeup_count] for m in expanded) == Set([1000.0])
    @test Set(m.stats[:wakeup_f] for m in expanded) == Set([1000.0])
    @test Set(m.parameters[:wakeup_count] for m in expanded) == Set([1000.0])
    @test Set(m.parameters[:wakeup_f] for m in expanded) == Set([1000.0])
    @test all(m -> m.stats[:fatigue_count] == 0, expanded)
    @test all(m -> !haskey(m.stats, :fatigue_f), expanded)
    @test all(m -> !haskey(m.stats, :fatigue_V), expanded)

    wakeup_voltages = Set(m.stats[:wakeup_V] for m in expanded)
    @test wakeup_voltages == Set([1.5, 2.5, 3.5])
    @test Set(m.parameters[:wakeup_V] for m in expanded) == Set([1.5, 2.5, 3.5])
    @test Set(
        (
            m.measurement_kind,
            m.stats[:wakeup_V],
            m.stats[:V_base],
            m.stats[:V_min],
            m.stats[:V_max],
            m.stats[:V_amp],
        )
        for m in expanded
    ) == Set([
        (:wakeup_pn, 1.5, 0.0, -1.5, 1.5, 1.5),
        (:wakeup_pn, 2.5, 0.0, -2.5, 2.5, 2.5),
        (:wakeup_pn, 3.5, 0.0, -3.5, 3.5, 3.5),
        (:wakeup_pund, 1.5, 0.0, -1.5, 1.5, 1.5),
        (:wakeup_pund, 2.5, 0.0, -2.5, 2.5, 2.5),
        (:wakeup_pund, 3.5, 0.0, -3.5, 3.5, 3.5),
    ])

    # Every expanded item shares the same source file
    @test all(m -> m.filepath == _WAKEUP_FIXTURE, expanded)
end

@testset "pund_wakeup count accumulation across files" begin
    @test isfile(_WAKEUP_ACCUM_FIRST)
    @test isfile(_WAKEUP_ACCUM_SECOND)

    source = MeasurementBrowser.scan_source(
        joinpath(@__DIR__, "fixtures", "RuO2");
        project=RUO2_PROJECT,
    )
    measurements = source.hierarchy.all_measurements

    first = only(filter(
        m -> m.filepath == _WAKEUP_ACCUM_FIRST &&
             m.measurement_kind === :wakeup_pund &&
             m.stats[:wakeup_V] == 3.5,
        measurements,
    ))
    second = only(filter(
        m -> m.filepath == _WAKEUP_ACCUM_SECOND &&
             m.measurement_kind === :wakeup_pund &&
             m.stats[:wakeup_V] == 3.5,
        measurements,
    ))

    @test first.stats[:wakeup_count] == 1000.0
    @test second.stats[:wakeup_count] == 2000.0
end

@testset "pund_wakeup expanded IDs and parameters" begin
    expanded = measurements_for_file(RUO2_PROJECT, _WAKEUP_FIXTURE)

    @test all(m -> startswith(m.unique_id, _WAKEUP_FIXTURE * "#wakeup_V="), expanded)
    @test length(unique(m.unique_id for m in expanded)) == length(expanded)
    @test all(m -> m.device_info.location == ["RuO2test_A11", "XI", "FeCap", "A9"], expanded)
end

@testset "pund_wakeup data access" begin
    expanded = measurements_for_file(RUO2_PROJECT, _WAKEUP_FIXTURE)

    pn_item   = first(filter(m -> m.measurement_kind === :wakeup_pn,   expanded))
    pund_item = first(filter(m -> m.measurement_kind === :wakeup_pund, expanded))
    data = MeasurementBrowser.read_measurement_data(RUO2_PROJECT, [pn_item, pund_item])

    @test length(data) == 2
    @test all(df -> nrow(df) > 0, data)
    @test names(data[1]) == names(data[2])
end
