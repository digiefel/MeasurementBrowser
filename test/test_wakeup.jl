using Test
using MeasurementBrowser

const _WAKEUP_FIXTURE = joinpath(
    @__DIR__,
    "fixtures", "RuO2",
    "RuO2test_A11_XI_FeCap_A9_20260509_231259_PUND_Wakeup.csv",
)

@testset "pund_wakeup detection and single MeasurementInfo" begin
    @test isfile(_WAKEUP_FIXTURE)
    @test detect_kind(RUO2_PROJECT, basename(_WAKEUP_FIXTURE)) == :pund_wakeup

    meas = MeasurementInfo(_WAKEUP_FIXTURE, RUO2_PROJECT)
    @test meas.measurement_kind == :pund_wakeup
    @test meas.device_info.location == ["RuO2test_A11", "XI", "FeCap", "A9"]
end

@testset "pund_wakeup expansion" begin
    expanded = measurements_for_file(RUO2_PROJECT, _WAKEUP_FIXTURE)

    # 3 amplitudes × 2 segment kinds (pn + pund) = 6 items
    @test length(expanded) == 6

    pn_items   = filter(m -> m.measurement_kind === :wakeup_pn,   expanded)
    pund_items = filter(m -> m.measurement_kind === :wakeup_pund, expanded)
    @test length(pn_items)   == 3
    @test length(pund_items) == 3

    amplitudes = Set(m.parameters[:amplitude_V] for m in expanded)
    @test amplitudes == Set([1.5, 2.5, 3.5])

    # Every expanded item shares the same source file
    @test all(m -> m.filepath == _WAKEUP_FIXTURE, expanded)
end

@testset "pund_wakeup display labels" begin
    expanded = measurements_for_file(RUO2_PROJECT, _WAKEUP_FIXTURE)

    item_1p5_pn   = only(filter(m -> m.measurement_kind === :wakeup_pn   && m.parameters[:amplitude_V] == 1.5, expanded))
    item_1p5_pund = only(filter(m -> m.measurement_kind === :wakeup_pund && m.parameters[:amplitude_V] == 1.5, expanded))

    @test kind_label(RUO2_PROJECT, :wakeup_pn)   == "Wakeup PN"
    @test kind_label(RUO2_PROJECT, :wakeup_pund) == "Wakeup PUND"

    label_pn   = display_label(RUO2_PROJECT, item_1p5_pn)
    label_pund = display_label(RUO2_PROJECT, item_1p5_pund)
    @test occursin("Wakeup", label_pn)
    @test occursin("1.5", label_pn)
    @test occursin("PN", label_pn)
    @test occursin("Wakeup", label_pund)
    @test occursin("1.5", label_pund)
    @test occursin("PUND", label_pund)
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
