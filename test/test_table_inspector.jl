using MeasurementBrowser
using MeasurementBrowser: inspect_table
using DataFrames: DataFrame
using Test

const Browser = MeasurementBrowser.Browser
const TI = MeasurementBrowser.TableInspector

@testset "table inspector" begin

    # --- legacy buffer helpers ---
    mktempdir() do dir
        buffer = fill(UInt8(0), 8)
        Browser._set_buffer_string!(buffer, "abc")
        @test Browser._buffer_string(buffer) == "abc"
        Browser._set_buffer_string!(buffer, "abcdefghij")
        @test Browser._buffer_string(buffer) == "abcdefg"
        @test buffer[end] == UInt8(0)

        comma_path = joinpath(dir, "with_preamble.csv")
        write(
            comma_path,
            "# instrument note\n" *
            "sweep id: 3\n" *
            "time_s,current_A,voltage_V\n" *
            "0.0,1e-9,0.0\n" *
            "1.0,2e-9,0.1\n",
        )

        preview = inspect_table(comma_path)
        @test preview.delimiter == ','
        @test preview.header_row == 3
        @test preview.data_start_row == 4
        @test preview.columns == ["time_s", "current_A", "voltage_V"]
        @test preview.row_count == 2
        @test preview.preview_rows == 2

        # All-rows mode: no cap, warnings empty
        preview_all = inspect_table(comma_path; max_rows=typemax(Int))
        @test preview_all.row_count == 2
        @test preview_all.preview_rows == 2
        @test isempty(preview_all.warnings)

        tab_path = joinpath(dir, "no_header.tsv")
        write(tab_path, "1\t2\t3\n4\t5\t6\n")
        preview = inspect_table(tab_path; max_rows=1)
        @test preview.delimiter == '\t'
        @test preview.header_row === nothing
        @test preview.data_start_row == 1
        @test preview.row_count == 2
        @test preview.preview_rows == 1
        @test !isempty(preview.warnings)
    end

    # --- merge_item_tables: single item ---
    @testset "merge_item_tables single item" begin
        df = DataFrame(x=[1.0, 2.0, 3.0], y=[4.0, 5.0, 6.0])
        pairs = [("item_a", df)]
        table, warnings = TI.merge_item_tables(
            Tuple{Any,Any}[(lbl, d) for (lbl, d) in pairs],
        )
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
        table, warnings = TI.merge_item_tables(pairs)
        @test isempty(warnings)
        # Column union: x, y, z (x first from df1, then y, then z from df2)
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

    # --- merge_item_tables: non-DataFrame skipped ---
    @testset "merge_item_tables non-DataFrame warning" begin
        df = DataFrame(a=[1, 2])
        pairs = Tuple{Any,Any}[("good", df), ("bad", "not a dataframe")]
        table, warnings = TI.merge_item_tables(pairs)
        @test length(warnings) == 1
        @test occursin("non-tabular", warnings[1])
        @test table.rows == 2
        @test table.item_labels == ["good"]
    end

    # --- merge_item_tables: empty input ---
    @testset "merge_item_tables empty" begin
        table, warnings = TI.merge_item_tables(Tuple{Any,Any}[])
        @test table.rows == 0
        @test isempty(table.columns)
        @test isempty(table.item_labels)
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

        # Ctrl+A pattern (select all)
        sel = Int[]
        for r in all_rows
            Browser._update_multi_selection!(sel, r, all_rows, false, true)
        end
        @test Set(sel) == Set(all_rows)
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

    # --- file-mode grid model builder ---
    @testset "file grid model from TablePreview" begin
        mktempdir() do dir
            path = joinpath(dir, "sample.csv")
            write(path, "a,b,c\n1,2,3\n4,5,6\n7,8,9\n")
            # Load all rows
            preview = inspect_table(path; max_rows=typemax(Int))
            @test preview.row_count == 3
            @test preview.preview_rows == 3
            @test isempty(preview.warnings)

            columns, n_rows, cell = Browser._file_grid_model(preview)
            @test columns == ["a", "b", "c"]
            @test n_rows == 3
            # First row, first column
            @test cell(1, 1) == "1"
            @test cell(1, 2) == "2"
            @test cell(3, 3) == "9"
        end
    end

    # --- PersistedProjectView round-trip (no table_column_widths) ---
    @testset "PersistedProjectView round-trip" begin
        view = Browser.PersistedProjectView(project="test")
        toml_data = Browser._project_view_to_toml(view)
        @test haskey(toml_data, "project")
        @test toml_data["project"] == "test"
        # table_column_widths field is gone — ensure no key pollution
        @test !haskey(toml_data, "table_column_widths")

        restored = Browser._project_view_from_toml(Browser.PersistedProjectView, toml_data)
        @test restored.project == "test"
    end
end
