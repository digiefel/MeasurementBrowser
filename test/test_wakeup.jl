using Test
using DataFrames: nrow
using MeasurementBrowser

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

    source = index_source_file(_WAKEUP_FIXTURE)
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

    @test all(m -> haskey(m.parameters, :wakeup_count), expanded)
    @test all(m -> haskey(m.parameters, :wakeup_f), expanded)
    @test all(m -> haskey(m.parameters, :wakeup_V), expanded)
    @test all(m -> haskey(m.parameters, :fatigue_count), expanded)
    @test all(m -> haskey(m.parameters, :fatigue_f), expanded)
    @test all(m -> haskey(m.parameters, :fatigue_V), expanded)
    @test Set(get(m.parameters, :wakeup_count, NaN) for m in expanded) == Set([1000.0])
    @test Set(get(m.parameters, :wakeup_f, NaN) for m in expanded) == Set([1000.0])
    @test all(m -> get(m.parameters, :fatigue_count, NaN) == 0, expanded)
    @test all(m -> isnan(get(m.parameters, :fatigue_f, NaN)), expanded)
    @test all(m -> isnan(get(m.parameters, :fatigue_V, NaN)), expanded)

    wakeup_voltages = Set(get(m.parameters, :wakeup_V, NaN) for m in expanded)
    @test wakeup_voltages == Set([1.5, 2.5, 3.5])
    @test Set(
        (
            m.parameters[:wakeup_V],
            m.stats[:V_base],
            m.stats[:V_min],
            m.stats[:V_max],
            m.stats[:V_amp],
        )
        for m in expanded
    ) == Set([
        (1.5, 0.0, -1.5, 1.5, 1.5),
        (2.5, 0.0, -2.5, 2.5, 2.5),
        (3.5, 0.0, -3.5, 3.5, 3.5),
    ])

    # Every expanded item shares the same source file
    @test all(m -> m.filepath == _WAKEUP_FIXTURE, expanded)
end

@testset "pund_wakeup count accumulation across files" begin
    @test isfile(_WAKEUP_ACCUM_FIRST)
    @test isfile(_WAKEUP_ACCUM_SECOND)

    source = scan_source(joinpath(@__DIR__, "fixtures", "RuO2"); project=RUO2_PROJECT)
    measurements = source.hierarchy.all_measurements

    first = only(filter(
        m -> m.filepath == _WAKEUP_ACCUM_FIRST &&
             m.measurement_kind === :wakeup_pund &&
             m.parameters[:wakeup_V] == 3.5,
        measurements,
    ))
    second = only(filter(
        m -> m.filepath == _WAKEUP_ACCUM_SECOND &&
             m.measurement_kind === :wakeup_pund &&
             m.parameters[:wakeup_V] == 3.5,
        measurements,
    ))

    @test first.parameters[:wakeup_count] == 1000.0
    @test second.parameters[:wakeup_count] == 2000.0
end

@testset "pund_wakeup expanded IDs and parameters" begin
    expanded = measurements_for_file(RUO2_PROJECT, _WAKEUP_FIXTURE)

    @test all(m -> startswith(m.unique_id, _WAKEUP_FIXTURE * "#wakeup_V="), expanded)
    @test length(unique(m.unique_id for m in expanded)) == length(expanded)
    @test all(m -> m.device_info.location == ["RuO2test_A11", "XI", "FeCap", "A9"], expanded)
end

@testset "pund_wakeup staged plot pipeline" begin
    expanded = measurements_for_file(RUO2_PROJECT, _WAKEUP_FIXTURE)

    pn_item   = first(filter(m -> m.measurement_kind === :wakeup_pn,   expanded))
    pund_item = first(filter(m -> m.measurement_kind === :wakeup_pund, expanded))
    pn_params   = merge(pn_item.device_info.parameters,   pn_item.parameters)
    pund_params = merge(pund_item.device_info.parameters, pund_item.parameters)

    loaded_pn = MeasurementBrowser.load_plot_for_file(
        RUO2_PROJECT, _WAKEUP_FIXTURE, :wakeup_pn; device_params=pn_params,
    )
    @test loaded_pn !== nothing
    @test nrow(loaded_pn.df) > 0

    loaded_pund = MeasurementBrowser.load_plot_for_file(
        RUO2_PROJECT, _WAKEUP_FIXTURE, :wakeup_pund; device_params=pund_params,
    )
    @test loaded_pund !== nothing
    @test nrow(loaded_pund.df) > 0

    analyzed = MeasurementBrowser.analyze_plot_for_file(
        RUO2_PROJECT, :wakeup_pund, loaded_pund,
    )
    @test analyzed !== nothing

    fig = MeasurementBrowser.draw_plot_for_file(RUO2_PROJECT, :wakeup_pund, analyzed)
    @test fig !== nothing
end
