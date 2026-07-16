import CImGui as ig
import CImGui.CSyntax: @c
using GLMakie: Axis, Figure, lines!, scatter!

import DataBrowserGUI
const Browser = DataBrowserGUI.Browser

using DataBrowserCore: InspectorTable, merge_item_tables
import DataBrowserCore.Workspace
using DataBrowserAPI: item_data
using .MakieImguiIntegration: MakieFigure

function _column_combo!(
    label::String,
    current::Int,
    columns::Vector{String},
    width::Real,
)::Int
    isempty(columns) && return current
    current = clamp(current, 1, length(columns))
    ig.SetNextItemWidth(width)
    if ig.BeginCombo(label, columns[current])
        for (index, column) in enumerate(columns)
            if ig.Selectable(column, index == current)
                current = index
            end
        end
        ig.EndCombo()
    end
    return current
end

function _table_plot_vectors(
    table::InspectorTable,
    x_col::Int,
    y_col::Int,
)::Tuple{Vector{Float64},Vector{Float64}}
    x_values = Float64[]
    y_values = Float64[]
    for row in 1:table.rows
        x_value = table.getvalue(row, x_col)
        y_value = table.getvalue(row, y_col)
        x_value isa Real && y_value isa Real || continue
        x_float = Float64(x_value)
        y_float = Float64(y_value)
        isfinite(x_float) && isfinite(y_float) || continue
        push!(x_values, x_float)
        push!(y_values, y_float)
    end
    return x_values, y_values
end

"""Clear the cached selection table so the next frame rebuilds it."""
function _clear_table_plot_table!(table_plot::TablePlotState)::Nothing
    table_plot.table = nothing
    table_plot.table_key = nothing
    table_plot.table_error = ""
    return nothing
end

"""
Return the merged table for the current selection, rebuilding only when the selection changes.

Errors are cached under the same key (shown via `table_error`) so a failing build does not
retry — and re-log — every frame.
"""
function _sync_table_plot_table!(
    state::Browser.BrowserState,
    table_plot::TablePlotState,
)::Union{Nothing,InspectorTable}
    workspace = state.workspace
    if !(workspace isa Workspace.Workspace)
        _clear_table_plot_table!(table_plot)
        return nothing
    end

    _, selected_records, _ = Browser._project_visible_selection(state)
    if isempty(selected_records)
        _clear_table_plot_table!(table_plot)
        return nothing
    end

    key = Tuple(sort!([r.id for r in selected_records]))
    table_plot.table_key == key && return table_plot.table

    table_plot.table_key = key
    table_plot.table_error = ""
    table_plot.table = try
        materialized = Workspace.materialize_items(workspace, selected_records)
        pairs = Tuple{Any,Any}[]
        for i in 1:length(selected_records)
            record = get(workspace.index.items, selected_records[i].id, nothing)
            label = record !== nothing ? record.label : string(materialized[i])
            data = try
                item_data(materialized[i])
            catch err
                bt = catch_backtrace()
                @error "Table plot: failed to load item data" label exception=(err, bt)
                continue
            end
            push!(pairs, (label, data))
        end
        table, _ = merge_item_tables(pairs)
        table
    catch err
        bt = catch_backtrace()
        @error "Table plot: failed to build table" exception=(err, bt)
        table_plot.table_error = first(split(sprint(showerror, err), '\n'; limit=2))
        table_plot.figure = nothing
        table_plot.plot_key = nothing
        nothing
    end
    return table_plot.table
end

function _ensure_table_plot!(
    state::Browser.BrowserState,
    table_plot::TablePlotState,
    table::InspectorTable,
)::Nothing
    table.rows == 0 && return nothing
    length(table.columns) < 2 && return nothing

    x_index = clamp(table_plot.x_column, 1, length(table.columns))
    y_index = clamp(table_plot.y_column, 1, length(table.columns))
    plot_key = (table_plot.table_key, x_index, y_index)
    table_plot.plot_key == plot_key && table_plot.figure !== nothing && return nothing

    x_col_name = table.columns[x_index]
    y_col_name = table.columns[y_index]
    x_values, y_values = _table_plot_vectors(table, x_index, y_index)

    if isempty(x_values)
        table_plot.figure = nothing
        table_plot.plot_key = nothing
        table_plot.plot_error = "No finite numeric points for the selected columns."
        return nothing
    end

    fig = Figure(size=(600, 400))
    ax = Axis(fig[1, 1], xlabel=x_col_name, ylabel=y_col_name, title="Table Plot")
    lines!(ax, x_values, y_values; linewidth=1)
    scatter!(ax, x_values, y_values; markersize=4)
    table_plot.figure = fig
    table_plot.plot_key = plot_key
    table_plot.plot_error = ""
    return nothing
end

function render_table_plot_menu!(state::Browser.BrowserState, table_plot::TablePlotState)::Nothing
    if ig.BeginMenu("Plot")
        if ig.MenuItem("Table Plot", C_NULL, table_plot.visible)
            table_plot.visible = !table_plot.visible
        end
        ig.EndMenu()
    end
    return nothing
end

function render_table_plot_window!(
    state::Browser.BrowserState,
    table_plot::TablePlotState,
)::Nothing
    table_plot.visible || return nothing

    open_ref = Ref(true)
    ig.SetNextWindowSize((700, 500), ig.ImGuiCond_FirstUseEver)
    if ig.Begin("Table Plot", open_ref)
        table = _sync_table_plot_table!(state, table_plot)

        if !isempty(table_plot.table_error)
            ig.TextColored((1.0, 0.35, 0.35, 1.0), table_plot.table_error)
        elseif table isa InspectorTable && table.rows > 0 && length(table.columns) >= 2
            width = max(80.0, (ig.GetContentRegionAvail().x - 12.0) / 2.0)
            table_plot.x_column = _column_combo!(
                "X##table_plot_x",
                table_plot.x_column,
                table.columns,
                width,
            )
            ig.SameLine()
            table_plot.y_column = _column_combo!(
                "Y##table_plot_y",
                table_plot.y_column,
                table.columns,
                width,
            )
            _ensure_table_plot!(state, table_plot, table)
            if table_plot.figure !== nothing
                MakieFigure(
                    "table_plot_window",
                    table_plot.figure;
                    auto_resize_x=true,
                    auto_resize_y=true,
                )
            elseif !isempty(table_plot.plot_error)
                ig.TextDisabled(table_plot.plot_error)
            end
        else
            ig.TextDisabled("Select items with tabular data to plot columns.")
        end
    end
    open_ref[] || (table_plot.visible = false)
    ig.End()
    return nothing
end
