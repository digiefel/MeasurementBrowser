@setup_workload begin
    fixture_dir = normpath(joinpath(@__DIR__, "..", "test", "fixtures", "RuO2"))
    pund_path = joinpath(fixture_dir, "3V FE PUND [RuO2test_A9_VI_D1(2) ; 2025-10-01 17_12_33].csv")
    tlm_path  = joinpath(fixture_dir, "TLM_4P [RuO2test_A9_VI_TLML100W2(12) ; 2025-10-01 16_21_45].csv")
    scan_fixture_dir = mktempdir()
    cp(pund_path, joinpath(scan_fixture_dir, basename(pund_path)); force=true)
    cp(tlm_path,  joinpath(scan_fixture_dir, basename(tlm_path));  force=true)

    @compile_workload begin
        pund_meas = nothing
        tlm_meas  = nothing

        if isfile(pund_path)
            pund_meas = only(measurements_for_file(RUO2_PROJECT, pund_path))
            pund_params = merge(pund_meas.device_info.parameters, pund_meas.parameters)
            pund_df = only(data_of_measurements(RUO2_PROJECT, [pund_meas]))
            loaded = _ruo2_plot_data(pund_meas, pund_df)
            analyzed = _analyze_ruo2_file_plot(
                RUO2_PROJECT,
                RuO2PUNDPlot,
                loaded;
                device_params=pund_params,
            )
            _draw_ruo2_file_plot(RUO2_PROJECT, RuO2PUNDPlot, analyzed; device_params=pund_params)
        end

        if isfile(tlm_path)
            tlm_meas = only(measurements_for_file(RUO2_PROJECT, tlm_path))
            tlm_params = merge(tlm_meas.device_info.parameters, tlm_meas.parameters)
            tlm_df = only(data_of_measurements(RUO2_PROJECT, [tlm_meas]))
            loaded = _ruo2_plot_data(tlm_meas, tlm_df)
            analyzed = _analyze_ruo2_file_plot(
                RUO2_PROJECT,
                RuO2TLM4PointPlot,
                loaded;
                device_params=tlm_params,
            )
            _draw_ruo2_file_plot(
                RUO2_PROJECT,
                RuO2TLM4PointPlot,
                analyzed;
                device_params=tlm_params,
            )
        end

        if pund_meas !== nothing && tlm_meas !== nothing
            infer_measurement_group(
                "precompile_group",
                MeasurementInfo[pund_meas],
                MeasurementInfo[pund_meas, tlm_meas],
            )
        end

        scan_source(scan_fixture_dir; project=RUO2_PROJECT)
        scan_source(scan_fixture_dir; project=RUO2_PROJECT, count_first=true)
    end
end
