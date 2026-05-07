using MeasurementBrowser
using DataLoader
using Test

function _write_csv(path::String, body::String)
    mkpath(dirname(path))
    open(path, "w") do io
        write(io, body)
    end
end

@testset "RuO2 CVSweep support" begin
    @test detect_kind(RUO2_PROJECT, "RuO2test_A8_VI_FeCap_B1_20260320_115801_298K_CVSweep.csv") == :cvsweep
    @test kind_label(RUO2_PROJECT, :cvsweep) == "C-V Sweep"

    mktempdir() do dir
        g_path = joinpath(
            dir,
            "RuO2test_A8",
            "VI",
            "FeCap",
            "B1",
            "RuO2test_A8_VI_FeCap_B1_20260320_115801_298K_CVSweep.csv",
        )
        _write_csv(
            g_path,
            """
Frequency_Hz,Bias_V,Cp (F),G (S),Time_sec,Status_Cp,Status_G,Status_Combined
1000.0,1.0,5.0e-12,1.0e-9,0.1,0,0,0
1000.0,-1.0,6.0e-12,1.5e-9,0.2,0,0,0
10000.0,1.0,4.0e-12,2.0e-9,0.3,0,0,0
10000.0,-1.0,3.0e-12,2.5e-9,0.4,0,0,0
""",
        )

        rp_path = joinpath(
            dir,
            "RuO2test_A8",
            "VI",
            "30um",
            "C3",
            "RuO2test_A8_VI_30um_C3_20260320_015049_298K_CVSweep.csv",
        )
        _write_csv(
            rp_path,
            """
Frequency_Hz,Bias_V,Cp (F),Rp (Ohm),Time_sec,Status_Cp,Status_Rp,Status_Combined
1000.0,2.0,4.9e-11,1.2e8,0.1,0,0,0
1000.0,-2.0,5.1e-11,1.1e8,0.2,0,0,0
10000.0,2.0,4.7e-11,1.5e7,0.3,0,0,0
10000.0,-2.0,4.8e-11,1.6e7,0.4,0,0,0
""",
        )

        no_status_path = joinpath(
            dir,
            "RuO2test_A9",
            "VI",
            "FeCap",
            "E4",
            "RuO2test_A9_VI_FeCap_E4_20260505_172807_CVSweep.csv",
        )
        _write_csv(
            no_status_path,
            """
Frequency_Hz,Bias_V,Cp (F),G (S),Time_sec
1000.0,1.0,5.0e-12,1.0e-9,0.1
1000.0,-1.0,6.0e-12,1.5e-9,0.2
""",
        )

        headerless_path = joinpath(
            dir,
            "RuO2test_A8",
            "VI",
            "25um",
            "E1",
            "RuO2test_A8_VI_25um_E1_20260316_175153_298K_CVSweep.csv",
        )
        _write_csv(
            headerless_path,
            """
Bias_V,Cp (F),Rp (Ohm),Time_sec,Status_Cp,Status_Rp,Status_Combined
-3.0,4.93e-11,9.2e6,0.1,0,0,0
3.0,4.67e-11,1.3e7,0.2,0,0,0
""",
        )

        g_items = MeasurementBrowser.interpret_file(RUO2_PROJECT, index_csv_file(g_path))
        @test length(g_items) == 1
        @test only(g_items).kind == :cvsweep
        @test only(g_items).device_path == ["RuO2test_A8", "VI", "FeCap", "B1"]

        @test isempty(MeasurementBrowser.interpret_file(RUO2_PROJECT, index_csv_file(headerless_path)))

        g_loaded = read_cv_sweep(basename(g_path), dirname(g_path))
        @test g_loaded.secondary_kind == :conductance
        @test names(g_loaded.df) == [
            "frequency_Hz",
            "bias_V",
            "cp_F",
            "secondary_value",
            "time_s",
            "status_cp",
            "status_secondary",
            "status_combined",
        ]
        @test sort(unique(g_loaded.df.frequency_Hz)) == [1000.0, 10000.0]

        no_status_loaded = read_cv_sweep(basename(no_status_path), dirname(no_status_path))
        @test all(no_status_loaded.df.status_cp .== 0)
        @test all(no_status_loaded.df.status_secondary .== 0)
        @test all(no_status_loaded.df.status_combined .== 0)

        g_plot_loaded = MeasurementBrowser.load_plot_for_file(RUO2_PROJECT, g_path, :cvsweep)
        g_plot_analyzed = MeasurementBrowser.analyze_plot_for_file(RUO2_PROJECT, :cvsweep, g_plot_loaded)
        @test g_plot_analyzed.secondary_kind == :conductance
        @test g_plot_analyzed.secondary_label == "G (nS)"
        @test g_plot_analyzed.frequencies_Hz == [1000.0, 10000.0]
        @test MeasurementBrowser.draw_plot_for_file(RUO2_PROJECT, :cvsweep, g_plot_analyzed) !== nothing

        rp_plot_loaded = MeasurementBrowser.load_plot_for_file(RUO2_PROJECT, rp_path, :cvsweep)
        rp_plot_analyzed = MeasurementBrowser.analyze_plot_for_file(RUO2_PROJECT, :cvsweep, rp_plot_loaded)
        @test rp_plot_analyzed.secondary_kind == :parallel_resistance
        @test rp_plot_analyzed.secondary_label == "Rp (MΩ)"
        @test MeasurementBrowser.draw_plot_for_file(RUO2_PROJECT, :cvsweep, rp_plot_analyzed) !== nothing
    end
end
