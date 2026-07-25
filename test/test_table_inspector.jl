using DataBrowserGUI
using DataBrowser: inspect_table
using DataFrames: DataFrame
using Test

const Browser = DataBrowserGUI.Browser
using DataBrowserCore: merge_item_tables

@testset "table inspector" begin
    fixture_root = joinpath(@__DIR__, "fixtures", "tables")
    preview = inspect_table(joinpath(fixture_root, "with_preamble.dat"))
    @test preview.delimiter == ','
    @test preview.header_row == 3
    @test preview.columns == ["time_s", "current_A", "voltage_V"]
    @test preview.row_count == 2
    @test isempty(preview.warnings)

    headerless = inspect_table(joinpath(fixture_root, "no_header.dat"); max_rows=1)
    @test headerless.delimiter == '\t'
    @test headerless.header_row === nothing
    @test headerless.row_count == 2
    @test headerless.preview_rows == 1
    @test !isempty(headerless.warnings)

    # --- merge_item_tables: single item ---
    @testset "merge_item_tables single item" begin
        df = DataFrame(x=[1.0, 2.0, 3.0], y=[4.0, 5.0, 6.0])
        pairs = [("item_a", df)]
        table, warnings = merge_item_tables(Tuple{Any,Any}[(lbl, d) for (lbl, d) in pairs])
        @test isempty(warnings)
        @test table.columns == ["x", "y"]
        @test table.rows == 3
        @test length(table.row_item) == 3
        @test all(table.row_item .== 1)
        @test table.item_labels == ["item_a"]
        # Check cell values
        @test table.getcell(1, 1) == "1.0"
        @test table.getcell(1, 2) == "4.0"
        @test table.getcell(3, 2) == "6.0"
    end

    # --- merge_item_tables: multiple items, column union ---
    @testset "merge_item_tables multi-item column union" begin
        df1 = DataFrame(x=[1, 2], y=[3, 4])
        df2 = DataFrame(y=[5, 6], z=[7, 8])  # y shared, z new, x missing
        pairs = Tuple{Any,Any}[("A", df1), ("B", df2)]
        table, warnings = merge_item_tables(pairs)
        @test isempty(warnings)
        @test Set(table.columns) == Set(["x", "y", "z"])
        @test table.rows == 4  # 2 from df1 + 2 from df2
        @test table.item_labels == ["A", "B"]
        # Provenance: rows 1-2 are from item 1, rows 3-4 from item 2
        @test table.row_item[1:2] == [1, 1]
        @test table.row_item[3:4] == [2, 2]
        # df2 row 1 (table row 3): z column = 7, x column = "" (missing)
        x_col = findfirst(==("x"), table.columns)
        z_col = findfirst(==("z"), table.columns)
        @test table.getcell(3, x_col) == ""   # df2 has no x
        @test table.getcell(3, z_col) == "7"
    end

    # --- merge_item_tables: empty input ---
    @testset "merge_item_tables empty" begin
        table, warnings = merge_item_tables(Tuple{Any,Any}[])
        @test isempty(warnings)
        @test table.rows == 0
        @test isempty(table.columns)
        @test isempty(table.item_labels)
    end

    # --- merge_item_tables: non-tabular items are skipped, not fatal ---
    @testset "merge_item_tables non-tabular skip" begin
        df = DataFrame(x=[1, 2], y=[3, 4])
        pairs = Tuple{Any,Any}[
            ("good", df),
            ("bad", Dict(:not => "a table")),
            ("matrix", [1 2; 3 4]),
        ]
        table, warnings = merge_item_tables(pairs)
        @test table.rows == 2
        @test table.item_labels == ["good"]
        @test length(warnings) == 2
        @test any(occursin("'bad'", w) for w in warnings)
        @test any(occursin("'matrix'", w) for w in warnings)
    end

    # --- InspectorTable: typed cell access ---
    @testset "InspectorTable typed cell access" begin
        df = DataFrame(a=Float32[1.5, 2.5], b=[1 // 2, 3 // 4])
        table, _ = merge_item_tables(Tuple{Any,Any}[("i", df)])
        @test table.getvalue(1, 1) === 1.5f0
        @test table.getvalue(2, 2) === 3 // 4
        # Display text stays a string; typed access must not go through it
        @test table.getcell(1, 1) == "1.5f0"

        # A column absent from one merged item reads as missing there
        df2 = DataFrame(c=[9])
        merged, _ = merge_item_tables(Tuple{Any,Any}[("i", df), ("j", df2)])
        c_col = findfirst(==("c"), merged.columns)
        @test merged.getvalue(1, c_col) === missing
        @test merged.getvalue(3, c_col) === 9
    end

    # --- _update_multi_selection!: basic operations ---
    @testset "multi_selection Int rows" begin
        all_rows = [1, 2, 3, 4, 5]

        # Single click replaces
        sel = Int[]
        Browser._update_multi_selection!(sel, 3, all_rows, false, false)
        @test sel == [3]

        # Second single click replaces
        Browser._update_multi_selection!(sel, 5, all_rows, false, false)
        @test sel == [5]

        # Ctrl+click toggles on
        Browser._update_multi_selection!(sel, 3, all_rows, false, true)
        @test 3 in sel && 5 in sel

        # Ctrl+click toggles off
        Browser._update_multi_selection!(sel, 3, all_rows, false, true)
        @test !(3 in sel)
        @test 5 in sel

        # Shift+click from anchor 5 to 2: range 2..5
        sel = [5]
        Browser._update_multi_selection!(sel, 2, all_rows, true, false)
        @test Set(sel) == Set([2, 3, 4, 5])

    end

    @testset "cell range clipboard text" begin
        state = Browser.DataGridState(
            cell_anchor=(3, 2),
            cell_focus=(1, 3),
        )
        values = [
            "r1c1" "r1c2" "with\ttab";
            "r2c1" "r2c2" "with \"quote\"";
            "r3c1" "r3c2" "r3c3";
        ]
        cell(row, column)::String = values[row, column]

        @test Browser._selected_cell_bounds(state, 3, 3) == (1, 3, 2, 3)
        @test Browser.selected_cell_text(state, 3, 3, cell) ==
            "r1c2\t\"with\ttab\"\n" *
            "r2c2\t\"with \"\"quote\"\"\"\n" *
            "r3c2\tr3c3"

        state.cell_anchor = (20, -2)
        state.cell_focus = (20, -2)
        @test Browser._selected_cell_bounds(state, 3, 3) == (3, 3, 1, 1)
        @test Browser.selected_cell_text(state, 3, 3, cell) == "r3c1"
    end

end
