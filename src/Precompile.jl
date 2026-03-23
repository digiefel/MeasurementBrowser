@setup_workload begin
    fixture_dir = normpath(joinpath(@__DIR__, "..", "test"))
    pund_path = joinpath(fixture_dir, "3V FE PUND [RuO2test_A9_VI_D1(2) ; 2025-10-01 17_12_33].csv")
    wakeup_path = joinpath(fixture_dir, "Wakeup 3V [RuO2test_A9_VI_D1(2) ; 2025-10-01 17_10_48].csv")
    tlm_path = joinpath(fixture_dir, "TLM_4P [RuO2test_A9_VI_TLML100W2(12) ; 2025-10-01 16_21_45].csv")

    @compile_workload begin
        if isfile(pund_path)
            pund_meas = MeasurementInfo(pund_path, RUO2_PROJECT)
            pund_params = merge(pund_meas.device_info.parameters, pund_meas.parameters)
            loaded = load_plot_for_file(RUO2_PROJECT, pund_path, :pund; device_params=pund_params)
            analyzed = analyze_plot_for_file(RUO2_PROJECT, :pund, loaded; device_params=pund_params)
            draw_plot_for_file(RUO2_PROJECT, :pund, analyzed; device_params=pund_params)
        end

        if isfile(wakeup_path)
            loaded = load_plot_for_file(RUO2_PROJECT, wakeup_path, :wakeup)
            analyzed = analyze_plot_for_file(RUO2_PROJECT, :wakeup, loaded)
            draw_plot_for_file(RUO2_PROJECT, :wakeup, analyzed)
        end

        if isfile(tlm_path)
            tlm_meas = MeasurementInfo(tlm_path, RUO2_PROJECT)
            tlm_params = merge(tlm_meas.device_info.parameters, tlm_meas.parameters)
            loaded = load_plot_for_file(RUO2_PROJECT, tlm_path, :tlm4p; device_params=tlm_params)
            analyzed = analyze_plot_for_file(RUO2_PROJECT, :tlm4p, loaded; device_params=tlm_params)
            draw_plot_for_file(RUO2_PROJECT, :tlm4p, analyzed; device_params=tlm_params)
        end
    end
end
