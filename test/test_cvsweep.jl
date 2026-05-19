using MeasurementBrowser
using DataLoader
using DataFrames: nrow
using Test

const _CV_FIXTURES = (
    z_comment = "RuO2test_A13_XI_FeCap_B1_20260511_004120_CVSweep.csv",
    g_plain   = "RuO2test_A9_VI_FeCap_A4_20260320_181744_298K_CVSweep.csv",
)

@testset "RuO2 CVSweep support" begin
    @test detect_kind(RUO2_PROJECT, "RuO2test_A8_VI_FeCap_B1_20260320_115801_298K_CVSweep.csv") == :cvsweep

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
