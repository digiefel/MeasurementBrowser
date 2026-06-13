using DataFrames: nrow
using GLMakie: Axis, Figure, lines!, scatter!
import CImGui as ig
using NativeFileDialog: pick_file

using ..TableInspector:
    TablePreview,
    inspect_table
using .MakieImguiIntegration: MakieFigure

"""Read one table into the inspector state and reset view-local plot choices."""
function _inspect_table_path!(
    state::BrowserState,
    path::AbstractString,
)::Nothing
    inspector = state.table_inspector
    _set_buffer_string!(inspector.path_buffer, String(path))
    try
        preview = inspect_table(path)
        _set_buffer_string!(inspector.path_buffer, preview.path)
        inspector.preview = preview
        inspector.error = ""
        inspector.x_column = 1
        inspector.y_column = min(2, length(preview.columns))
        inspector.figure = nothing
        inspector.plot_key = nothing
        inspector.plot_error = ""
    catch err
        bt = catch_backtrace()
        inspector.preview = nothing
        inspector.error = first(split(sprint(showerror, err), '\n'; limit=2))
        inspector.figure = nothing
        inspector.plot_key = nothing
        @error("Table inspection failed\nPath: $(String(path))", exception=(err, bt))
    end
    inspector.visible = true
    return nothing
end

"""Update the inspector from browser selection when Live is enabled."""
function _sync_live_table_source!(state::BrowserState)::Nothing
    inspector = state.table_inspector
    inspector.live || return nothing
    selected_path = _selected_table_source(state)
    selected_path === nothing && return nothing
    current_path = strip(_buffer_string(inspector.path_buffer))
    !isempty(current_path) && normpath(current_path) == normpath(selected_path) &&
        return nothing
    _inspect_table_path!(state, selected_path)
    return nothing
end

"""Return the first selected measurement source file, if a source selection exists."""
function _selected_table_source(state::BrowserState)::Union{Nothing,String}
    _, selected_measurements, _ = _project_visible_selection(state)
    isempty(selected_measurements) && return nothing
    return first(selected_measurements).filepath
end

"""Render the menu commands that open the table inspector."""
function _render_table_inspector_menu!(state::BrowserState)::Nothing
    if ig.BeginMenu("Inspect")
        if ig.MenuItem("Table Inspector", C_NULL, state.table_inspector.visible)
            state.table_inspector.visible = !state.table_inspector.visible
        end
        ig.EndMenu()
    end
    return nothing
end

"""Return display text for one table cell."""
function _table_cell_text(value::Any)::String
    text = sprint(show, value)
    return length(text) > 90 ? first(text, 87) * "..." : text
end

"""Render the preview dataframe as a scrollable ImGui table."""
function _render_table_preview(preview::TablePreview)::Nothing
    table = preview.table
    columns = preview.columns
    if isempty(columns) || nrow(table) == 0
        ig.TextDisabled("No rows to preview")
        return nothing
    end

    flags = ig.ImGuiTableFlags_RowBg |
            ig.ImGuiTableFlags_Borders |
            ig.ImGuiTableFlags_Resizable |
            ig.ImGuiTableFlags_ScrollX |
            ig.ImGuiTableFlags_ScrollY
    if ig.BeginTable("table_inspector_preview", length(columns), flags, (0, 360))
        for column in columns
            ig.TableSetupColumn(column, ig.ImGuiTableColumnFlags_WidthFixed, 140)
        end
        ig.TableHeadersRow()
        for row in 1:nrow(table)
            ig.TableNextRow()
            for column in columns
                ig.TableNextColumn()
                ig.TextUnformatted(_table_cell_text(table[row, column]))
            end
        end
        ig.EndTable()
    end
    return nothing
end

"""Choose one dataframe column with a compact combo box."""
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

"""Convert two preview columns into finite float vectors for plotting."""
function _table_plot_vectors(
    preview::TablePreview,
    x_column::String,
    y_column::String,
)::Tuple{Vector{Float64},Vector{Float64}}
    x_values = Float64[]
    y_values = Float64[]
    for (x, y) in zip(preview.table[!, x_column], preview.table[!, y_column])
        x_float = try
            x isa Real ? Float64(x) : parse(Float64, string(x))
        catch
            continue
        end
        y_float = try
            y isa Real ? Float64(y) : parse(Float64, string(y))
        catch
            continue
        end
        isfinite(x_float) && isfinite(y_float) || continue
        push!(x_values, x_float)
        push!(y_values, y_float)
    end
    return x_values, y_values
end

"""Build or reuse the quick plot for the current X/Y column choices."""
function _ensure_table_plot!(inspector::TableInspectorState)::Nothing
    preview = inspector.preview
    preview isa TablePreview || return nothing
    columns = preview.columns
    length(columns) >= 2 || return nothing
    x_index = clamp(inspector.x_column, 1, length(columns))
    y_index = clamp(inspector.y_column, 1, length(columns))
    plot_key = (preview.path, x_index, y_index, preview.preview_rows)
    inspector.plot_key == plot_key && inspector.figure !== nothing && return nothing

    x_column = columns[x_index]
    y_column = columns[y_index]
    x_values, y_values = _table_plot_vectors(preview, x_column, y_column)
    if isempty(x_values)
        inspector.figure = nothing
        inspector.plot_key = nothing
        inspector.plot_error = "No finite numeric points for the selected columns."
        return nothing
    end

    fig = Figure(size=(650, 360))
    ax = Axis(fig[1, 1], xlabel=x_column, ylabel=y_column, title="Table Preview")
    lines!(ax, x_values, y_values; linewidth=1)
    scatter!(ax, x_values, y_values; markersize=4)
    inspector.figure = fig
    inspector.plot_key = plot_key
    inspector.plot_error = ""
    return nothing
end

"""Render quick plot controls and the plot generated from the preview dataframe."""
function _render_table_quick_plot!(inspector::TableInspectorState)::Nothing
    preview = inspector.preview
    preview isa TablePreview || return nothing
    columns = preview.columns
    length(columns) < 2 && return nothing

    width = max(96.0, (ig.GetContentRegionAvail().x - 16.0) / 2.0)
    inspector.x_column = _column_combo!(
        "X##table_inspector_x",
        inspector.x_column,
        columns,
        width,
    )
    ig.SameLine()
    inspector.y_column = _column_combo!(
        "Y##table_inspector_y",
        inspector.y_column,
        columns,
        width,
    )
    _ensure_table_plot!(inspector)
    if inspector.figure !== nothing
        MakieFigure(
            "table_inspector_plot",
            inspector.figure;
            auto_resize_x=true,
            auto_resize_y=false,
        )
    elseif !isempty(inspector.plot_error)
        ig.TextDisabled(inspector.plot_error)
    end
    return nothing
end

"""Render the generic table inspector window."""
function render_table_inspector_window(state::BrowserState)::Nothing
    inspector = state.table_inspector
    inspector.visible || return nothing

    open_ref = Ref(true)
    ig.SetNextWindowSize((900, 720), ig.ImGuiCond_FirstUseEver)
    if ig.Begin("Table Inspector", open_ref, ig.ImGuiWindowFlags_NoDocking)
        _sync_live_table_source!(state)
        ig.InputText(
            "##table_inspector_path",
            inspector.path_buffer,
            length(inspector.path_buffer),
        )
        ig.SameLine()
        live = inspector.live
        if @c ig.Checkbox("Live##table_inspector_live", &live)
            inspector.live = live
            _sync_live_table_source!(state)
        end
        ig.SameLine()
        if ig.Button("Open...")
            path = pick_file()
            if !isnothing(path) && !isempty(path)
                inspector.live = false
                _inspect_table_path!(state, path)
            end
        end
        ig.SameLine()
        selected_path = _selected_table_source(state)
        selected_path === nothing && ig.BeginDisabled()
        if ig.Button("Open Selection")
            _inspect_table_path!(state, selected_path)
        end
        selected_path === nothing && ig.EndDisabled()
        if selected_path === nothing && ig.BeginItemTooltip()
            ig.TextUnformatted("Select a measurement first.")
            ig.EndTooltip()
        end
        ig.SameLine()
        if ig.Button("Reload")
            path = strip(_buffer_string(inspector.path_buffer))
            isempty(path) || _inspect_table_path!(state, path)
        end

        if !isempty(inspector.error)
            ig.TextColored((1.0, 0.35, 0.35, 1.0), inspector.error)
        end

        preview = inspector.preview
        if preview isa TablePreview
            delimiter = preview.delimiter == '\t' ? "\\t" : string(preview.delimiter)
            header = preview.header_row === nothing ? "none" : string(preview.header_row)
            ig.TextDisabled(
                "Delimiter: $delimiter   Header row: $header   Data starts: " *
                "$(preview.data_start_row)   Rows: $(preview.row_count)",
            )
            for warning in preview.warnings
                ig.TextDisabled(warning)
            end
            ig.Separator()
            _render_table_preview(preview)
            ig.Separator()
            _render_table_quick_plot!(inspector)
        else
            ig.TextDisabled("Open a delimited text file or select a measurement source.")
        end
    end
    open_ref[] || (inspector.visible = false)
    ig.End()
    return nothing
end
