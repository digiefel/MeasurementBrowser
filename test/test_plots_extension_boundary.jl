using DataBrowserPlots
using DataBrowserCore: merge_item_tables
using DataFrames: DataFrame
using Test

@testset "table plot typed extraction" begin
    # Float32 (and Rational) columns must plot via typed values, not display strings.
    df = DataFrame(x=Float32[1.5, 2.5, 3.5], y=[1 // 2, 3 // 4, 5 // 4])
    table, _ = merge_item_tables(Tuple{Any,Any}[("i", df)])
    x, y = DataBrowserPlots._table_plot_vectors(table, 1, 2)
    @test x == [1.5, 2.5, 3.5]
    @test y == [0.5, 0.75, 1.25]

    # Non-numeric and missing cells are skipped, not fatal
    df_mixed = DataFrame(x=[1.0, 2.0], label=["a", "b"])
    mixed, _ = merge_item_tables(Tuple{Any,Any}[("i", df_mixed)])
    x2, y2 = DataBrowserPlots._table_plot_vectors(mixed, 1, 2)
    @test isempty(x2) && isempty(y2)
end
