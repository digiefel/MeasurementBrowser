using MeasurementBrowser
using DataFrames: nrow
using Test

const RUO2_PROJECT = MeasurementBrowser.RUO2_PROJECT
const RuO2PUNDPlot = MeasurementBrowser.RuO2PUNDPlot
const RuO2IVSweepPlot = MeasurementBrowser.RuO2IVSweepPlot
const RuO2TLM4PointPlot = MeasurementBrowser.RuO2TLM4PointPlot
const RuO2TLMAnalysisPlot = MeasurementBrowser.RuO2TLMAnalysisPlot
const RuO2TLMTemperaturePlot = MeasurementBrowser.RuO2TLMTemperaturePlot
const RuO2PUNDFatiguePlot = MeasurementBrowser.RuO2PUNDFatiguePlot

@testset "RuO2 simple plot API" begin
    fixtures_dir = joinpath(@__DIR__, "fixtures", "RuO2")
    pund_path = joinpath(fixtures_dir, "3V FE PUND [RuO2test_A9_VI_D1(2) ; 2025-10-01 17_12_33].csv")
    iv_path = joinpath(fixtures_dir, "RuO2test_A11_XI_FeCapBD_A1A2_20260509_184021_IVSweep.csv")
    tlm_device_iv_path =
        joinpath(fixtures_dir, "RuO2test_A11_XI_TLM_L100W4_20260609_183048_14K_IVSweep.csv")
    tlm_path = joinpath(fixtures_dir, "TLM_4P [RuO2test_A9_VI_TLML100W2(12) ; 2025-10-01 16_21_45].csv")
    fatigue_path = joinpath(fixtures_dir, "RuO2test_A12_XI_FeCap_D6_20260511_222203_PUND_Fatigue_V2.csv")

    pund = only(measurements_for_file(RUO2_PROJECT, pund_path))
    iv = only(measurements_for_file(RUO2_PROJECT, iv_path))
    tlm_device_iv = only(measurements_for_file(RUO2_PROJECT, tlm_device_iv_path))
    tlm = only(measurements_for_file(RUO2_PROJECT, tlm_path))
    fatigue = measurements_for_file(RUO2_PROJECT, fatigue_path)[1:3]

    function tlm_measurement(length_um, temperature_K=nothing)
        parameters = deepcopy(tlm.parameters)
        suffix = "L$(Int(length_um))"
        if temperature_K !== nothing
            parameters[:temperature_K] = Float64(temperature_K)
            parameters[:site] = "test_site"
            parameters[:oxygen_percent] = 10.0
            parameters[:oxygen_flow_sccm] = 5.0
            suffix *= "_$(Int(temperature_K))K"
        end
        return MeasurementInfo(
            tlm;
            unique_id="$(tlm.unique_id):$suffix",
            clean_title="$(tlm.clean_title) $suffix",
            parameters,
            device_info=DeviceInfo(
                copy(tlm.device_info.location),
                merge(tlm.device_info.parameters, Dict{Symbol,Any}(:length_um => Float64(length_um), :width_um => 2.0)),
            ),
        )
    end

    tlm_l100 = tlm_measurement(100.0)
    tlm_l50 = tlm_measurement(50.0)
    tlm_analysis_measurements = [tlm_l100, tlm_l50]
    tlm_temperature_measurements = [
        tlm_measurement(100.0, 298.0),
        tlm_measurement(50.0, 298.0),
        tlm_measurement(100.0, 350.0),
        tlm_measurement(50.0, 350.0),
    ]
    workspace = MeasurementBrowser.Workspace.Workspace(RUO2_PROJECT, fixtures_dir)

    pund_data = only(process_measurement_data(workspace, [pund]))
    @test all(name in names(pund_data) for name in [
        "time_us",
        "current_uA",
        "i_fe_uA",
        "q_fe_pC",
        "q_centered_pC",
        "pulse_group_size",
    ])
    @test nrow(pund_data) > 0

    iv_data = only(process_measurement_data(workspace, [iv]))
    @test names(iv_data) == ["v", "i_abs", "i_positive"]
    @test nrow(iv_data) > 0

    tlm_device_iv_data = only(process_measurement_data(workspace, [tlm_device_iv]))
    @test tlm_device_iv.measurement_kind === :iv
    @test maximum(abs, tlm_device_iv_data.v) ≈ 0.001
    @test maximum(tlm_device_iv_data.i_abs) < 1e-6

    tlm_data = only(process_measurement_data(workspace, [tlm]))
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

    pund_fig = setup_plot(workspace, RuO2PUNDPlot, [pund, pund])
    @test plot_data!(workspace, RuO2PUNDPlot, [pund, pund], pund_fig) === nothing

    iv_fig = setup_plot(workspace, RuO2IVSweepPlot, [iv, iv])
    @test plot_data!(workspace, RuO2IVSweepPlot, [iv, iv], iv_fig) === nothing
    tlm_device_iv_fig = setup_plot(workspace, RuO2IVSweepPlot, [tlm_device_iv])
    @test plot_data!(
        workspace,
        RuO2IVSweepPlot,
        [tlm_device_iv],
        tlm_device_iv_fig,
    ) === nothing

    tlm_fig = setup_plot(workspace, RuO2TLM4PointPlot, [tlm, tlm])
    @test plot_data!(workspace, RuO2TLM4PointPlot, [tlm, tlm], tlm_fig) === nothing

    tlm_analysis_fig = setup_plot(workspace, RuO2TLMAnalysisPlot, tlm_analysis_measurements)
    @test_logs (:info, r"TLM combined analysis completed") (:info, r"Sheet/contact analysis") begin
        @test plot_data!(workspace, RuO2TLMAnalysisPlot, tlm_analysis_measurements, tlm_analysis_fig) === nothing
    end

    tlm_temperature_fig = setup_plot(workspace, RuO2TLMTemperaturePlot, tlm_temperature_measurements)
    @test_logs (:info, r"TLM combined analysis completed") (:info, r"Sheet/contact analysis") (:info, r"TLM combined analysis completed") (:info, r"Sheet/contact analysis") begin
        @test plot_data!(workspace, RuO2TLMTemperaturePlot, tlm_temperature_measurements, tlm_temperature_fig) === nothing
    end

    fatigue_fig = setup_plot(workspace, RuO2PUNDFatiguePlot, fatigue)
    @test plot_data!(workspace, RuO2PUNDFatiguePlot, fatigue, fatigue_fig) === nothing
end
