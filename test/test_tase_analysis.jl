using MeasurementBrowser
using Test
using DataFrames: nrow

const TASE_PROJECT = MeasurementBrowser.TASE_PROJECT
const TASEFourTerminalIVPlot = MeasurementBrowser.TASEFourTerminalIVPlot

@testset "TASE Analysis" begin
    fixture1 = joinpath(@__DIR__, "fixtures", "TASE", "TASESNS1c1f_A_2TSNJunction_11_20260224_111623_298K_FourTerminalIV.csv")
    fixture2 = joinpath(@__DIR__, "fixtures", "TASE", "TASESNS1c1f_A_2TSNJunction_31_20260224_111700_298K_FourTerminalIV.csv")
    source1 = MeasurementBrowser.index_source_file(fixture1)
    measurements = [
        only(measurements_for_file(TASE_PROJECT, fixture1)),
        only(measurements_for_file(TASE_PROJECT, fixture2)),
    ]

    @test parse_device_info(TASE_PROJECT, source1).location == measurements[1].device_info.location

    @testset "plot data api" begin
        workspace = MeasurementBrowser.Workspace.Workspace(TASE_PROJECT, dirname(fixture1))
        data = read_measurement_data(workspace, measurements)
        @test length(data) == 2
        @test all(nrow(df) == 3 for df in data)

        fig = setup_plot(workspace, TASEFourTerminalIVPlot, measurements)
        @test plot_data!(workspace, TASEFourTerminalIVPlot, measurements, fig) === nothing
    end
end
