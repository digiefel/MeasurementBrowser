using DataFrames: nrow, DataFrame, names
using GLMakie: Axis, Figure, lines!, scatter!
import CImGui as ig
using NativeFileDialog: pick_file

using ...DataBrowserSources

using ..TableInspector:
    InspectorTable,
    merge_item_tables
using .MakieImguiIntegration: MakieFigure

import ..Workspace
using DataBrowserAPI: item_data

# ---------------------------------------------------------------------------
# Buffer helpers (kept for the raw-preview mode)
# ---------------------------------------------------------------------------

"""Return the null-terminated text currently stored in an ImGui byte buffer."""
function _buffer_string(buffer::Vector{UInt8})::String
    terminator = findfirst(==(UInt8(0)), buffer)
    last_index = terminator === nothing ? length(buffer) : terminator - 1
    return String(buffer[1:last_index])
end

"""Write null-terminated text into an ImGui byte buffer."""
function _set_buffer_string!(buffer::Vector{UInt8}, text::AbstractString)::Nothing
    isempty(buffer) && return nothing
    fill!(buffer, UInt8(0))
    bytes = codeunits(String(text))
    byte_count = min(length(bytes), length(buffer) - 1)
    byte_count > 0 && copyto!(buffer, 1, bytes, 1, byte_count)
    return nothing
end

# ---------------------------------------------------------------------------
# Item-data source helpers
# ---------------------------------------------------------------------------

const _ITEM_TINT_COLORS = UInt32[
    0x20_4C_80_FF,  # blue-ish
    0x20_80_4C_FF,  # green-ish
    0x80_4C_20_FF,  # orange-ish
    0x80_20_4C_FF,  # purple-ish
    0x20_80_80_FF,  # teal-ish
    0x80_80_20_FF,  # yellow-ish
]

"""Return an ImGui-packed RGBA color for provenance tinting of item index `i`."""
function _provenance_tint(i::Int)::UInt32
    return _ITEM_TINT_COLORS[mod1(i, length(_ITEM_TINT_COLORS))]
end

"""
Resolve selected items and refresh the inspector table when the selection changes.

Materialization is synchronous and may be slow for very large items; for now we accept that
cost (same as the plot panel). A key based on item ids prevents redundant reloads.
"""
function _sync_item_data_inspector!(state::BrowserState)::Nothing
    inspector = state.table_inspector
    workspace = state.workspace
    if !(workspace isa Workspace.Workspace)
        inspector.inspector_table = nothing
        inspector.inspector_warnings = String[]
        inspector.inspector_key = nothing
        return nothing
    end

    _, selected_records, _ = _project_visible_selection(state)
    isempty(selected_records) && return nothing

    # Compute a cache key so we don't rebuild on every frame
    key = tuple(sort([r.id for r in selected_records])..., inspector.show_provenance_column)
    inspector.inspector_key == key && return nothing

    # Materialize (may load from cache or origin)
    materialized = try
        Workspace.materialize_items(workspace, selected_records)
    catch err
        bt = catch_backtrace()
        @error "Table inspector: failed to materialize items" exception=(err, bt)
        inspector.inspector_table = nothing
        inspector.inspector_warnings = ["Error loading items: $(first(split(sprint(showerror, err), '\n'; limit=2)))"]
        inspector.inspector_key = key
        return nothing
    end

    # Build label from record (more useful than the item object's show)
    pairs = [
        (materialized[i], get(workspace.index.items, selected_records[i].id, nothing))
        for i in 1:length(selected_records)
    ]
    labeled_pairs = Tuple{Any,Any}[]
    warnings = String[]
    for (mat_item, record) in pairs
        d = try
            data_raw = item_data(mat_item)
            data_raw isa DataFrame ? data_raw : nothing
        catch
            nothing
        end
        if d === nothing
            lbl = record !== nothing ? record.item_label : "item"
            push!(warnings, "Item '$lbl' has non-tabular data; skipped.")
            continue
        end
        lbl = record !== nothing ? record.item_label : string(mat_item)
        push!(labeled_pairs, (lbl, d))
    end

    # Build InspectorTable from labeled pairs (label is the "item" here)
    col_set = Set{String}()
    columns = String[]
    dfs = DataFrame[]
    labels = String[]
    for (lbl, df) in labeled_pairs
        cols = String.(names(df))
        for c in cols
            if c ∉ col_set
                push!(col_set, c)
                push!(columns, c)
            end
        end
        push!(dfs, df)
        push!(labels, lbl)
    end

    if isempty(dfs)
        inspector.inspector_table = InspectorTable(columns, 0, Int[], labels, (_, _) -> "")
        inspector.inspector_warnings = warnings
        inspector.inspector_key = key
        inspector.grid.selected_rows = Int[]
        return nothing
    end

    # Add provenance column if requested and multi-item
    show_prov = inspector.show_provenance_column && length(dfs) > 1
    if show_prov
        pushfirst!(columns, "_item_")
    end

    row_item = Int[]
    row_offsets = Int[]
    for (i, df) in enumerate(dfs)
        for r in 1:nrow(df)
            push!(row_item, i)
            push!(row_offsets, r)
        end
    end

    total_rows = length(row_item)
    col_indices = [
        Dict(c => j for (j, c) in enumerate(String.(names(df))))
        for df in dfs
    ]

    # Offset for provenance column
    prov_offset = show_prov ? 1 : 0

    function getcell(row::Int, col::Int)::String
        item_i = row_item[row]
        row_in_item = row_offsets[row]
        # Provenance column is col=1 when show_prov
        if show_prov && col == 1
            return labels[item_i]
        end
        actual_col = col - prov_offset
        actual_col < 1 && return ""
        df = dfs[item_i]
        col_name = columns[col]
        ci = get(col_indices[item_i], col_name, nothing)
        ci === nothing && return ""
        text = sprint(show, df[row_in_item, ci])
        return length(text) > 90 ? first(text, 87) * "..." : text
    end

    inspector.inspector_table = InspectorTable(columns, total_rows, row_item, labels, getcell)
    inspector.inspector_warnings = warnings
    inspector.inspector_key = key
    inspector.grid.selected_rows = Int[]
    inspector.figure = nothing
    inspector.plot_key = nothing

    # Track current kind to form a stable per-kind DataGrid table id (used by imgui.ini)
    kinds = unique([r.kind for r in selected_records])
    inspector.current_kind = length(kinds) == 1 ? first(kinds) : nothing

    return nothing
end

# ---------------------------------------------------------------------------
# Row tint helper
# ---------------------------------------------------------------------------

"""Return an ImGui packed color for provenance tinting (semi-transparent row overlay)."""
function _row_tint_for_item(row_item_index::Int)::UInt32
    colors = [
        ig.IM_COL32(100, 150, 255, 40),   # blue-ish
        ig.IM_COL32(100, 220, 150, 40),   # green-ish
        ig.IM_COL32(255, 160, 80,  40),   # orange-ish
        ig.IM_COL32(200, 100, 220, 40),   # purple-ish
        ig.IM_COL32(80,  210, 210, 40),   # teal-ish
        ig.IM_COL32(210, 200, 80,  40),   # yellow-ish
    ]
    return colors[mod1(row_item_index, length(colors))]
end

# ---------------------------------------------------------------------------
# Quick-plot helpers (driven from InspectorTable)
# ---------------------------------------------------------------------------

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

"""Extract finite float64 pairs from the inspector table for the given column indices."""
function _inspector_plot_vectors(
    table::InspectorTable,
    x_col::Int,
    y_col::Int,
    row_subset::Union{Nothing,Vector{Int}},
)::Tuple{Vector{Float64},Vector{Float64}}
    x_values = Float64[]
    y_values = Float64[]
    rows = row_subset !== nothing ? row_subset : (1:table.rows)
    for row in rows
        x_str = table.getcell(row, x_col)
        y_str = table.getcell(row, y_col)
        x_float = tryparse(Float64, x_str)
        y_float = tryparse(Float64, y_str)
        x_float === nothing && continue
        y_float === nothing && continue
        isfinite(x_float) && isfinite(y_float) || continue
        push!(x_values, x_float)
        push!(y_values, y_float)
    end
    return x_values, y_values
end

"""Build or reuse the quick plot from the inspector table."""
function _ensure_inspector_plot!(
    inspector::TableInspectorState,
)::Nothing
    table = inspector.inspector_table
    table isa InspectorTable || return nothing
    table.rows == 0 && return nothing
    columns = table.columns
    length(columns) < 2 && return nothing

    x_index = clamp(inspector.x_column, 1, length(columns))
    y_index = clamp(inspector.y_column, 1, length(columns))

    selected_rows = inspector.grid.selected_rows
    row_subset = (inspector.plot_selected_only && !isempty(selected_rows)) ?
        copy(selected_rows) : nothing
    plot_key = (inspector.inspector_key, x_index, y_index, inspector.plot_selected_only,
                row_subset !== nothing ? sort(row_subset) : nothing)
    inspector.plot_key == plot_key && inspector.figure !== nothing && return nothing

    x_col_name = columns[x_index]
    y_col_name = columns[y_index]
    x_values, y_values = _inspector_plot_vectors(table, x_index, y_index, row_subset)

    if isempty(x_values)
        inspector.figure = nothing
        inspector.plot_key = nothing
        inspector.plot_error = "No finite numeric points for the selected columns."
        return nothing
    end

    fig = Figure(size=(600, 400))
    ax = Axis(fig[1, 1], xlabel=x_col_name, ylabel=y_col_name, title="Inspector Plot")
    lines!(ax, x_values, y_values; linewidth=1)
    scatter!(ax, x_values, y_values; markersize=4)
    inspector.figure = fig
    inspector.plot_key = plot_key
    inspector.plot_error = ""
    return nothing
end

"""Render X/Y column chooser and the quick plot."""
function _render_inspector_plot!(inspector::TableInspectorState)::Nothing
    table = inspector.inspector_table
    table isa InspectorTable || return nothing
    columns = table.columns
    length(columns) < 2 && return nothing

    width = max(80.0, (ig.GetContentRegionAvail().x - 12.0) / 2.0)
    inspector.x_column = _column_combo!(
        "X##ti_x",
        inspector.x_column,
        columns,
        width,
    )
    ig.SameLine()
    inspector.y_column = _column_combo!(
        "Y##ti_y",
        inspector.y_column,
        columns,
        width,
    )

    selected_rows = inspector.grid.selected_rows
    can_plot_selected = !isempty(selected_rows)
    plot_sel = inspector.plot_selected_only
    if @c ig.Checkbox("Plot selected only##ti_pso", &plot_sel)
        inspector.plot_selected_only = plot_sel
        inspector.figure = nothing
        inspector.plot_key = nothing
    end
    if inspector.plot_selected_only && !can_plot_selected
        ig.SameLine()
        ig.TextDisabled("(select rows)")
    end

    _ensure_inspector_plot!(inspector)
    if inspector.figure !== nothing
        MakieFigure(
            "table_inspector_plot",
            inspector.figure;
            auto_resize_x=true,
            auto_resize_y=true,
        )
    elseif !isempty(inspector.plot_error)
        ig.TextDisabled(inspector.plot_error)
    end
    return nothing
end

# ---------------------------------------------------------------------------
# Raw file mode — second DataGrid mode
# ---------------------------------------------------------------------------

# Row count above which we show a soft warning in the file mode
const _FILE_MODE_LARGE_ROW_WARN = 100_000

"""Read one table into the inspector state (raw file mode, all rows)."""
function _inspect_table_path!(
    state::BrowserState,
    path::AbstractString,
)::Nothing
    inspector = state.table_inspector
    _set_buffer_string!(inspector.path_buffer, String(path))
    try
        # Full file load (default); explicit max_rows caps rows for tests or callers that need it
        preview = inspect_table(path)
        _set_buffer_string!(inspector.path_buffer, preview.path)
        inspector.preview = preview
        inspector.file_grid = DataGridState()
        inspector.error = ""
    catch err
        bt = catch_backtrace()
        inspector.preview = nothing
        inspector.error = first(split(sprint(showerror, err), '\n'; limit=2))
        @error("Table inspection failed\nPath: $(String(path))", exception=(err, bt))
    end
    inspector.visible = true
    return nothing
end

"""
Build a DataGrid model (columns + row_count + cell callback) from a `TabularFileSource`.

Returns `(columns::Vector{String}, n_rows::Int, cell::Function)`.
"""
function _file_grid_model(
    preview::TabularFileSource,
)::Tuple{Vector{String},Int,Function}
    table = preview.table
    columns = preview.columns
    n_rows = nrow(table)
    col_indices = Dict(c => j for (j, c) in enumerate(columns))

    function cell(row::Int, col::Int)::String
        col_name = columns[col]
        ci = get(col_indices, col_name, nothing)
        ci === nothing && return ""
        text = sprint(show, table[row, ci])
        return length(text) > 90 ? first(text, 87) * "..." : text
    end

    return columns, n_rows, cell
end

"""Render the raw-file mode via DataGrid (virtualized, all rows, no provenance)."""
function _render_file_mode!(inspector::TableInspectorState)::Nothing
    preview = inspector.preview
    preview isa TabularFileSource || return nothing

    delimiter = preview.delimiter == '\t' ? "\\t" : string(preview.delimiter)
    header    = preview.header_row === nothing ? "none" : string(preview.header_row)
    ig.TextDisabled(
        "Delimiter: $delimiter   Header row: $header   Data starts: " *
        "$(preview.data_start_row)   Rows: $(preview.row_count)",
    )

    if preview.row_count > _FILE_MODE_LARGE_ROW_WARN
        ig.TextColored(
            (1.0f0, 0.8f0, 0.2f0, 1.0f0),
            "Large file: $(preview.row_count) rows",
        )
    end

    ig.Separator()

    columns, n_rows, cell = _file_grid_model(preview)
    if isempty(columns) || n_rows == 0
        ig.TextDisabled("No rows to display")
        return nothing
    end

    render_data_grid!(
        "file",
        inspector.file_grid;
        n_rows,
        columns,
        cell,
    )
    return nothing
end

# ---------------------------------------------------------------------------
# Menu integration
# ---------------------------------------------------------------------------

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

# ---------------------------------------------------------------------------
# Main window
# ---------------------------------------------------------------------------

"""
Render the table inspector window.

Primary mode: shows the selected items' `.data` (as DataFrames), merged by column union
with per-row provenance tinting when multiple items are selected.

Secondary raw-file mode: accessible via the path bar / Open... button; shows an arbitrary
delimited text file without provenance.
"""
function render_table_inspector_window(state::BrowserState)::Nothing
    inspector = state.table_inspector
    inspector.visible || return nothing

    open_ref = Ref(true)
    ig.SetNextWindowSize((1100, 720), ig.ImGuiCond_FirstUseEver)
    if ig.Begin("Table Inspector", open_ref, ig.ImGuiWindowFlags_NoDocking)

        # --- sync item-data on every frame (fast no-op when key hasn't changed) ---
        _sync_item_data_inspector!(state)

        table = inspector.inspector_table
        has_item_data = table isa InspectorTable && table.rows > 0

        # --- toolbar: provenance toggle + raw-file controls ---
        if has_item_data && length(table.item_labels) > 1
            show_prov = inspector.show_provenance_column
            if @c ig.Checkbox("Provenance column##ti_prov", &show_prov)
                inspector.show_provenance_column = show_prov
                inspector.inspector_key = nothing  # force rebuild
            end
            ig.SameLine()
        end

        ig.PushItemWidth(300)
        ig.InputText(
            "##ti_path",
            inspector.path_buffer,
            length(inspector.path_buffer),
        )
        ig.PopItemWidth()
        ig.SameLine()
        live = inspector.live
        if @c ig.Checkbox("Live##ti_live", &live)
            inspector.live = live
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
        if ig.Button("Reload")
            path = strip(_buffer_string(inspector.path_buffer))
            if !isempty(path)
                _inspect_table_path!(state, path)
            end
        end

        if !isempty(inspector.error)
            ig.TextColored((1.0, 0.35, 0.35, 1.0), inspector.error)
        end

        # --- item-data view (primary) ---
        if has_item_data
            _render_item_data_view!(state, table)

        # --- raw file mode (secondary) ---
        elseif inspector.preview isa TabularFileSource
            _render_file_mode!(inspector)

        else
            # Show warnings and hint
            for w in inspector.inspector_warnings
                ig.TextDisabled(w)
            end
            if table isa InspectorTable && table.rows == 0 && !isempty(table.item_labels)
                ig.TextDisabled("No tabular data found in the selected items.")
            else
                ig.TextDisabled("Select items to view their data, or open a file with 'Open...'.")
            end
        end
    end
    open_ref[] || (inspector.visible = false)
    ig.End()
    return nothing
end

"""Render the primary item-data view: DataGrid on the left, plot on the right."""
function _render_item_data_view!(state::BrowserState, table::InspectorTable)::Nothing
    inspector = state.table_inspector

    # Warnings (non-tabular items skipped, etc.)
    for w in inspector.inspector_warnings
        ig.TextDisabled(w)
    end

    # Row count / large-data warning
    if table.rows >= 1_000_000
        ig.TextColored(
            (1.0f0, 0.8f0, 0.2f0, 1.0f0),
            "Large dataset: $(table.rows) rows",
        )
    end

    # Multi-item: show provenance legend
    multi_item = length(table.item_labels) > 1
    if multi_item
        avail_x = ig.GetContentRegionAvail().x
        label_width = max(80.0f0, avail_x / Float32(length(table.item_labels) + 1))
        for (i, lbl) in enumerate(table.item_labels)
            tint_u32 = _row_tint_for_item(i)
            r = Float32((tint_u32 >> 0)  & 0xFF) / 255.0f0
            g = Float32((tint_u32 >> 8)  & 0xFF) / 255.0f0
            b = Float32((tint_u32 >> 16) & 0xFF) / 255.0f0
            ig.TextColored((r, g, b, 1.0f0), lbl)
            i < length(table.item_labels) && ig.SameLine()
        end
    end

    avail = ig.GetContentRegionAvail()
    left_w = avail.x * 0.55f0

    # Left: data grid
    ig.BeginChild("##ti_table", (left_w, 0.0f0), false)

    row_tint = if multi_item
        (row::Int) -> _row_tint_for_item(table.row_item[row])
    else
        (_) -> nothing
    end

    # Stable per-kind table id so imgui.ini keys column widths per item kind.
    # "mixed" when multiple kinds are selected; kinds are always Symbols so String() is safe.
    grid_id = inspector.current_kind !== nothing ? string(inspector.current_kind) : "mixed"

    render_data_grid!(
        grid_id,
        inspector.grid;
        n_rows=table.rows,
        columns=table.columns,
        cell=table.getcell,
        row_tint,
        on_selection_change=(_) -> begin
            # invalidate plot cache when selection changes with plot_selected_only
            inspector.plot_key = nothing
        end,
    )

    ig.EndChild()
    ig.SameLine()

    # Right: quick plot
    ig.BeginChild("##ti_plot", (0.0f0, 0.0f0), false)
    _render_inspector_plot!(inspector)
    ig.EndChild()

    return nothing
end
