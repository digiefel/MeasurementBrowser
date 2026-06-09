using MeasurementBrowser
using DataFrames: nrow
using Test

@testset "RuO2 simple plot API" begin
    fixtures_dir = joinpath(@__DIR__, "fixtures", "RuO2")
    pund_path = joinpath(fixtures_dir, "3V FE PUND [RuO2test_A9_VI_D1(2) ; 2025-10-01 17_12_33].csv")
    iv_path = joinpath(fixtures_dir, "RuO2test_A11_XI_FeCapBD_A1A2_20260509_184021_IVSweep.csv")
    tlm_path = joinpath(fixtures_dir, "TLM_4P [RuO2test_A9_VI_TLML100W2(12) ; 2025-10-01 16_21_45].csv")

    pund = only(measurements_for_file(RUO2_PROJECT, pund_path))
    iv = only(measurements_for_file(RUO2_PROJECT, iv_path))
    tlm = only(measurements_for_file(RUO2_PROJECT, tlm_path))
    tlm_l100 = MeasurementInfo(
        tlm;
        unique_id="$(tlm.unique_id):L100",
        clean_title="$(tlm.clean_title) L100",
        device_info=DeviceInfo(
            copy(tlm.device_info.location),
            merge(tlm.device_info.parameters, Dict{Symbol,Any}(:length_um => 100.0, :width_um => 2.0)),
        ),
    )
    tlm_l50 = MeasurementInfo(
        tlm;
        unique_id="$(tlm.unique_id):L50",
        clean_title="$(tlm.clean_title) L50",
        device_info=DeviceInfo(
            copy(tlm.device_info.location),
            merge(tlm.device_info.parameters, Dict{Symbol,Any}(:length_um => 50.0, :width_um => 2.0)),
        ),
    )
    tlm_analysis_measurements = [tlm_l100, tlm_l50]

    @test MeasurementBrowser._plot_job_data(RUO2_PROJECT, RuO2PUNDPlot, [pund]) === nothing
    @test MeasurementBrowser._plot_job_data(RUO2_PROJECT, RuO2IVSweepPlot, [iv]) === nothing
    @test MeasurementBrowser._plot_job_data(RUO2_PROJECT, RuO2TLM4PointPlot, [tlm]) === nothing
    @test MeasurementBrowser._plot_job_data(RUO2_PROJECT, RuO2TLMAnalysisPlot, tlm_analysis_measurements) === nothing

    pund_data = only(process_measurement_data(RUO2_PROJECT, [pund]))
    @test all(name in names(pund_data) for name in [
        "time_us",
        "current_uA",
        "i_fe_uA",
        "q_fe_pC",
        "q_centered_pC",
        "pulse_group_size",
    ])
    @test nrow(pund_data) > 0

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

    pund_fig = setup_plot(RUO2_PROJECT, RuO2PUNDPlot, [pund, pund])
    @test plot_data!(RUO2_PROJECT, RuO2PUNDPlot, [pund, pund], pund_fig) === nothing

    iv_fig = setup_plot(RUO2_PROJECT, RuO2IVSweepPlot, [iv, iv])
    @test plot_data!(RUO2_PROJECT, RuO2IVSweepPlot, [iv, iv], iv_fig) === nothing

    tlm_fig = setup_plot(RUO2_PROJECT, RuO2TLM4PointPlot, [tlm, tlm])
    @test plot_data!(RUO2_PROJECT, RuO2TLM4PointPlot, [tlm, tlm], tlm_fig) === nothing

    tlm_analysis_fig = setup_plot(RUO2_PROJECT, RuO2TLMAnalysisPlot, tlm_analysis_measurements)
    @test_logs (:info, r"TLM combined analysis completed") (:info, r"Sheet/contact analysis") begin
        @test plot_data!(RUO2_PROJECT, RuO2TLMAnalysisPlot, tlm_analysis_measurements, tlm_analysis_fig) === nothing
    end
end
