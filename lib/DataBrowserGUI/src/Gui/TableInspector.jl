import CImGui as ig
import CImGui.CSyntax: @c

using DataBrowserCore: InspectorTable, merge_item_tables

import DataBrowserCore.Workspace
using DataBrowserAPI: item_data, label

# ---------------------------------------------------------------------------
# Item-data source helpers
# ---------------------------------------------------------------------------

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
    key = (workspace.scan.epoch, sort([r.id for r in selected_records]), inspector.show_provenance_column)
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
        label = record !== nothing ? record.label : string(mat_item)
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
    types = unique([r.type for r in selected_records])
    inspector.current_kind = length(types) == 1 ? label(only(types)) : nothing

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

Shows the selected items' data (Tables.jl-compatible payloads), merged by column union with
per-row provenance tinting when multiple items are selected.
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

        if has_item_data && length(table.item_labels) > 1
            show_prov = inspector.show_provenance_column
            if @c ig.Checkbox("Provenance column##ti_prov", &show_prov)
                inspector.show_provenance_column = show_prov
                inspector.inspector_key = nothing  # force rebuild
            end
        end

        if has_item_data
            _render_item_data_view!(state, table)
        else
            for w in inspector.inspector_warnings
                ig.TextDisabled(w)
            end
            if table isa InspectorTable && table.rows == 0 && !isempty(table.item_labels)
                ig.TextDisabled("No tabular data found in the selected items.")
            else
                ig.TextDisabled("Select items to view their data.")
            end
        end
    end
    open_ref[] || (inspector.visible = false)
    ig.End()
    return nothing
end

"""Render the item-data view as a full-width DataGrid."""
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
