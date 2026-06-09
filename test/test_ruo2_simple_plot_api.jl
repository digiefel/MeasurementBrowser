using MeasurementBrowser
using DataFrames: nrow
using Test

@testset "RuO2 simple plot API" begin
    fixtures_dir = joinpath(@__DIR__, "fixtures", "RuO2")
    iv_path = joinpath(fixtures_dir, "RuO2test_A11_XI_FeCapBD_A1A2_20260509_184021_IVSweep.csv")
    tlm_path = joinpath(fixtures_dir, "TLM_4P [RuO2test_A9_VI_TLML100W2(12) ; 2025-10-01 16_21_45].csv")

    iv = only(measurements_for_file(RUO2_PROJECT, iv_path))
    tlm = only(measurements_for_file(RUO2_PROJECT, tlm_path))

    @test MeasurementBrowser._plot_job_data(RUO2_PROJECT, RuO2IVSweepPlot, [iv]) === nothing
    @test MeasurementBrowser._plot_job_data(RUO2_PROJECT, RuO2TLM4PointPlot, [tlm]) === nothing

    iv_data = only(process_measurement_data(RUO2_PROJECT, [iv]))
    @test names(iv_data) == ["v", "i_abs", "i_positive"]
    @test nrow(iv_data) > 0

    tlm_data = only(process_measurement_data(RUO2_PROJECT, [tlm]))
    @test all(name in names(tlm_data) for name in [
        "current_uA",
        "voltage_mV",
        "resistance_kohm",
        "valid_resistance",
        "fit_resistance_ohm",
        "fit_offset_V",
        "rho_sheet",
    ])
    @test nrow(tlm_data) > 0

    iv_fig = setup_plot(RUO2_PROJECT, RuO2IVSweepPlot, [iv, iv])
    @test plot_data!(RUO2_PROJECT, RuO2IVSweepPlot, [iv, iv], iv_fig) === nothing

    tlm_fig = setup_plot(RUO2_PROJECT, RuO2TLM4PointPlot, [tlm, tlm])
    @test plot_data!(RUO2_PROJECT, RuO2TLM4PointPlot, [tlm, tlm], tlm_fig) === nothing
end
