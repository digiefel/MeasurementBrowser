using MeasurementBrowser
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
    @test MeasurementBrowser.cv_sweep_has_schema(z_path)
    z_df = MeasurementBrowser.read_cv_sweep(basename(z_path), dirname(z_path))
    @test names(z_df) == ["frequency_Hz", "bias_V", "Cp_F", "Z_Ohm", "time_s", "status_cp", "status_combined"]
    @test nrow(z_df) > 0
    @test all(z_df.Z_Ohm .> 0)
    @test all(z_df.Cp_F .> 0)

    # Fixture without comments, has G (S) — Z must be computed from G + Cp.
    g_path = joinpath(fixtures_dir, _CV_FIXTURES.g_plain)
    @test MeasurementBrowser.cv_sweep_has_schema(g_path)
    g_df = MeasurementBrowser.read_cv_sweep(basename(g_path), dirname(g_path))
    @test names(g_df) == ["frequency_Hz", "bias_V", "Cp_F", "Z_Ohm", "time_s", "status_cp", "status_combined"]
    @test nrow(g_df) > 0
    @test all(isfinite, g_df.Z_Ohm)
    @test all(g_df.Z_Ohm .> 0)

    # Full pipeline check: interpret, load, and plot for the comment-header fixture.
    items = measurements_for_file(MeasurementBrowser.RUO2_PROJECT, z_path)
    @test length(items) == 1
    @test only(items).measurement_kind == :cvsweep

    data = MeasurementBrowser.read_measurement_data(RUO2_PROJECT, items)
    @test length(data) == 1
    @test nrow(only(data)) == nrow(z_df)

    fig = MeasurementBrowser.setup_plot(RUO2_PROJECT, RuO2CVSweepPlot, items)
    @test MeasurementBrowser.plot_data!(RUO2_PROJECT, RuO2CVSweepPlot, items, fig) === nothing
end
