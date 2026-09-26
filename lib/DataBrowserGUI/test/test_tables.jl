using Test, DataBrowserGUI

@testset "table union preserves values and item provenance" begin
    table, warnings = merge_item_tables(Tuple{Any,Any}[
        ("A", (x=Float32[1, 2], y=[3, 4])), ("B", (y=[5], z=[6]))])
    @test isempty(warnings)
    x, y, z = [findfirst(==(name), table.columns) for name in ["x", "y", "z"]]
    @test table.getvalue(1, x) === 1.0f0
    @test table.getvalue(3, x) === missing
    @test table.getvalue(3, y) == 5
    @test table.getvalue(3, z) == 6
    @test table.item_labels[table.row_item] == ["A", "A", "B"]
end

@testset "row tables and mixed payloads" begin
    table, warnings = merge_item_tables([
        ("rows", [(x=1, y=2), (x=3, y=4)]),
        ("opaque", nothing),
        ("columns", (y=[5], z=[6])),
    ])
    @test length(warnings) == 1
    @test table.rows == 3
    @test table.item_labels[table.row_item] == ["rows", "rows", "columns"]
    x, y, z = [findfirst(==(name), table.columns) for name in ["x", "y", "z"]]
    @test [table.getvalue(row, y) for row in 1:table.rows] == [2, 4, 5]
    @test table.getvalue(2, x) == 3
    @test table.getvalue(1, z) === missing
end
