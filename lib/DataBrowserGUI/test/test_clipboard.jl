using Test, DataBrowserGUI
const Browser = DataBrowserGUI.Browser

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

    @test Browser.selected_cell_text(state, 3, 3, cell) ==
        "r1c2\t\"with\ttab\"\n" *
        "r2c2\t\"with \"\"quote\"\"\"\n" *
        "r3c2\tr3c3"

    state.cell_anchor = (20, -2)
    state.cell_focus = (20, -2)
    @test Browser.selected_cell_text(state, 3, 3, cell) == "r3c1"
end
