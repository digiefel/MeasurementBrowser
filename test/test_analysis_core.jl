using Test
using DataFrames

using MeasurementBrowser

@testset "analysis core" begin
    result = AnalysisResult(
        :demo,
        "Demo",
        :device,
        DataFrame(x=[1, 2], y=[3.0, 4.0]),
        [(kind=:table, label="Table")],
        Dict{Symbol,Any}(:source => :test),
    )

    @test result.key == :demo
    @test result.label == "Demo"
    @test result.row_kind == :device
    @test nrow(result.table) == 2
    @test result.view_presets == [(kind=:table, label="Table")]
    @test result.meta[:source] == :test

    @test MeasurementBrowser.available_analyses(RUO2_PROJECT, MeasurementInfo[]) == NamedTuple[]
    @test MeasurementBrowser.run_analysis(RUO2_PROJECT, :demo, MeasurementInfo[]) === nothing
    @test MeasurementBrowser.draw_analysis_view(result, (kind=:table,)) === nothing
end
