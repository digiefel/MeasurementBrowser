using Tables
import CImGui as ig
import CImGui.CSyntax: @c
using NativeFileDialog: pick_file

using DataBrowserSources

using DataBrowserCore: InspectorTable, merge_item_tables

import DataBrowserCore.Workspace
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

    # Build labels from records (more useful than the item object's show); skip items whose
    # data fails to load so one bad item never hides its siblings.
    warnings = String[]
    labeled_pairs = Tuple{Any,Any}[]
    for i in 1:length(selected_records)
        mat_item = materialized[i]
        record = get(workspace.index.items, selected_records[i].id, nothing)
        label = record !== nothing ? record.item_label : string(mat_item)
        data = try
            item_data(mat_item)
        catch err
            bt = catch_backtrace()
            @error "Table inspector: failed to load item data" label exception=(err, bt)
            push!(warnings, "Item '$label': failed to load data; skipped.")
            continue
        end
        push!(labeled_pairs, (label, data))
    end

    table, merge_warnings = try
        merge_item_tables(labeled_pairs)
    catch err
        bt = catch_backtrace()
        @error "Table inspector: failed to build table" exception=(err, bt)
        inspector.inspector_table = nothing
        inspector.inspector_warnings = [
            "Error building table: $(first(split(sprint(showerror, err), '\n'; limit=2)))",
        ]
        inspector.inspector_key = key
        return nothing
    end
    append!(warnings, merge_warnings)

    show_prov = inspector.show_provenance_column && length(table.item_labels) > 1
    if show_prov
        columns = vcat(["_item_"], table.columns)
        inner = table
        function getcell(row::Int, col::Int)::String
            col == 1 && return inner.item_labels[inner.row_item[row]]
            return inner.getcell(row, col - 1)
        end
        function getvalue(row::Int, col::Int)::Any
            col == 1 && return inner.item_labels[inner.row_item[row]]
            return inner.getvalue(row, col - 1)
        end
        table = InspectorTable(
            columns, inner.rows, inner.row_item, inner.item_labels, getcell, getvalue)
    end

    inspector.inspector_table = table
    inspector.inspector_warnings = warnings
    inspector.inspector_key = key
    inspector.grid.selected_rows = Int[]

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
    n_rows = Tables.rowcount(table)
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

"""Render the primary item-data view as a full-width DataGrid."""
function _render_item_data_view!(state::BrowserState, table::InspectorTable)::Nothing
    inspector = state.table_inspector

    for w in inspector.inspector_warnings
        ig.TextDisabled(w)
    end

    if table.rows >= 1_000_000
        ig.TextColored(
            (1.0f0, 0.8f0, 0.2f0, 1.0f0),
            "Large dataset: $(table.rows) rows",
        )
    end

    multi_item = length(table.item_labels) > 1
    if multi_item
        for (i, lbl) in enumerate(table.item_labels)
            tint_u32 = _row_tint_for_item(i)
            r = Float32((tint_u32 >> 0)  & 0xFF) / 255.0f0
            g = Float32((tint_u32 >> 8)  & 0xFF) / 255.0f0
            b = Float32((tint_u32 >> 16) & 0xFF) / 255.0f0
            ig.TextColored((r, g, b, 1.0f0), lbl)
            i < length(table.item_labels) && ig.SameLine()
        end
    end

    row_tint = if multi_item
        (row::Int) -> _row_tint_for_item(table.row_item[row])
    else
        (_) -> nothing
    end

    grid_id = inspector.current_kind !== nothing ? string(inspector.current_kind) : "mixed"

    render_data_grid!(
        grid_id,
        inspector.grid;
        n_rows=table.rows,
        columns=table.columns,
        cell=table.getcell,
        row_tint,
    )

    return nothing
end
