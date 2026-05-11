using MeasurementBrowser
using DataLoader
using DataFrames: nrow
using Test

const _CV_FIXTURES = (
    z_comment = "RuO2test_A13_XI_FeCap_B1_20260511_004120_CVSweep.csv",
    g_plain   = "RuO2test_A9_VI_FeCap_A4_20260320_181744_298K_CVSweep.csv",
)

function _write_csv(path::String, body::String)
    mkpath(dirname(path))
    open(path, "w") do io
        write(io, body)
    end
end

@testset "RuO2 CVSweep support" begin
    @test detect_kind(RUO2_PROJECT, "RuO2test_A8_VI_FeCap_B1_20260320_115801_298K_CVSweep.csv") == :cvsweep
    @test kind_label(RUO2_PROJECT, :cvsweep) == "C-V Sweep"

    fixtures_dir = joinpath(@__DIR__, "fixtures", "RuO2")

    # Fixture with '#'-prefixed comment header and a direct Z (Ohm) column.
    z_path = joinpath(fixtures_dir, _CV_FIXTURES.z_comment)
    @test cv_sweep_has_schema(z_path)
    z_df = read_cv_sweep(basename(z_path), dirname(z_path))
    @test names(z_df) == ["frequency_Hz", "bias_V", "Cp_F", "Z_Ohm", "time_s", "status_cp", "status_combined"]
    @test nrow(z_df) > 0
    @test all(z_df.Z_Ohm .> 0)
    @test all(z_df.Cp_F .> 0)

    # Fixture without comments, has G (S) — Z must be computed from G + Cp.
    g_path = joinpath(fixtures_dir, _CV_FIXTURES.g_plain)
    @test cv_sweep_has_schema(g_path)
    g_df = read_cv_sweep(basename(g_path), dirname(g_path))
    @test names(g_df) == ["frequency_Hz", "bias_V", "Cp_F", "Z_Ohm", "time_s", "status_cp", "status_combined"]
    @test nrow(g_df) > 0
    @test all(isfinite, g_df.Z_Ohm)
    @test all(g_df.Z_Ohm .> 0)

    mktempdir() do dir
        # Synthetic file with Rp (Ohm) — Z computed from Rp + Cp.
        rp_path = joinpath(
            dir, "RuO2test_A8", "VI", "30um", "C3",
            "RuO2test_A8_VI_30um_C3_20260320_015049_298K_CVSweep.csv",
        )
        _write_csv(rp_path, """
Frequency_Hz,Bias_V,Cp (F),Rp (Ohm),Time_sec,Status_Cp,Status_Combined
1000.0,2.0,4.9e-11,1.2e8,0.1,0,0
1000.0,-2.0,5.1e-11,1.1e8,0.2,0,0
10000.0,2.0,4.7e-11,1.5e7,0.3,0,0
10000.0,-2.0,4.8e-11,1.6e7,0.4,0,0
""")
        @test cv_sweep_has_schema(rp_path)
        rp_df = read_cv_sweep(basename(rp_path), dirname(rp_path))
        @test all(isfinite, rp_df.Z_Ohm)
        @test all(rp_df.Z_Ohm .> 0)

        # File missing Cp/Z/G/Rp must be rejected.
        bad_path = joinpath(dir, "bad", "RuO2test_A8_VI_25um_E1_20260316_175153_298K_CVSweep.csv")
        _write_csv(bad_path, """
Bias_V,Cp (F),Rp (Ohm),Time_sec
-3.0,4.93e-11,9.2e6,0.1
3.0,4.67e-11,1.3e7,0.2
""")
        @test !cv_sweep_has_schema(bad_path)
    end

    # Full pipeline check: interpret, analyze, draw for the comment-header fixture.
    items = measurements_for_file(MeasurementBrowser.RUO2_PROJECT, z_path)
    @test length(items) == 1
    @test only(items).measurement_kind == :cvsweep

    loaded = MeasurementBrowser.load_plot_for_file(RUO2_PROJECT, z_path, :cvsweep)
    analyzed = MeasurementBrowser.analyze_plot_for_file(RUO2_PROJECT, :cvsweep, loaded)
    @test hasproperty(analyzed, :df)
    @test hasproperty(analyzed, :frequencies_Hz)
    @test !isempty(analyzed.frequencies_Hz)
    @test MeasurementBrowser.draw_plot_for_file(RUO2_PROJECT, :cvsweep, analyzed) !== nothing
end
