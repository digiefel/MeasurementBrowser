import CImGui as ig
import CImGui.CSyntax: @c

# ---------------------------------------------------------------------------
# Shared multi-select helper
# ---------------------------------------------------------------------------

"""
Apply range, toggle, or replacement selection to a list of items.

Generic over T so it can be used with item ids (String), row indices (Int), etc.
`shift_held` extends the range from the last selected item; `ctrl_held` toggles the item.
"""
function _update_multi_selection!(
    selected_items::Vector{T},
    item::T,
    all_items::Vector{T},
    shift_held::Bool,
    ctrl_held::Bool,
)::Nothing where {T}
    if shift_held && !isempty(selected_items)
        last_item = selected_items[end]
        start_idx = findfirst(x -> x == last_item, all_items)
        end_idx = findfirst(x -> x == item, all_items)
        if start_idx !== nothing && end_idx !== nothing
            start_idx > end_idx && ((start_idx, end_idx) = (end_idx, start_idx))
            append!(selected_items, all_items[start_idx:end_idx])
            unique!(selected_items)
        end
    elseif ctrl_held
        if item in selected_items
            filter!(x -> x != item, selected_items)
        else
            push!(selected_items, item)
        end
    else
        empty!(selected_items)
        push!(selected_items, item)
    end
    return nothing
end

# ---------------------------------------------------------------------------
# Render API
# ---------------------------------------------------------------------------

# NOTE: DataGridState is defined in Browser/State.jl (included before this file)
# so the struct is available here without redefinition.

"""
Render a virtualized, multi-select table grid.

Arguments:
- `id`: stable ImGui identifier string (no `##` prefix needed; the function adds it)
- `state`: mutable `DataGridState` owned by the consumer; written back each frame
- `n_rows`: total number of data rows
- `columns`: ordered column header names
- `cell`: `(row::Int, col::Int) -> String` callback for cell text
- `row_tint`: `(row::Int) -> Union{Nothing,UInt32}` for provenance colouring (ImGui packed color)
- `on_selection_change`: called with the new `Vector{Int}` whenever the selection changes
- `height`: child height in pixels; 0.0f0 = fill available

Column widths are persisted across restarts by ImGui via the ini file (keyed by `id`).
On first appearance, `ImGuiTableFlags_SizingFixedFit` auto-sizes each column to its content.

Multi-select: click = replace, Shift+click = extend range, Ctrl/Cmd+click = toggle,
Ctrl/Cmd+A = select all, ↑/↓ = move (Shift+↑/↓ = extend), Escape = clear.
"""
function render_data_grid!(
    id::String,
    state::DataGridState;
    n_rows::Int,
    columns::Vector{String},
    cell::Function,
    row_tint::Function = _ -> nothing,
    on_selection_change::Function = identity,
    height::Float32 = 0.0f0,
)::Nothing
    isempty(columns) && return nothing
    n_cols = length(columns)

    flags = ig.ImGuiTableFlags_RowBg          |
            ig.ImGuiTableFlags_Borders        |
            ig.ImGuiTableFlags_Resizable      |
            ig.ImGuiTableFlags_ScrollX        |
            ig.ImGuiTableFlags_ScrollY        |
            ig.ImGuiTableFlags_SizingFixedFit

    if !ig.BeginTable("##datagrid_$id", n_cols, flags, (0.0f0, height))
        return nothing
    end

    # Sticky header (freeze first row)
    ig.TableSetupScrollFreeze(0, 1)

    # Column headers — width 0 with SizingFixedFit = auto-fit on first appearance;
    # imgui.ini restores saved widths on subsequent opens.
    for col in columns
        ig.TableSetupColumn(col, ig.ImGuiTableColumnFlags_WidthFixed, 0.0f0)
    end
    ig.TableHeadersRow()

    # Keyboard navigation (only when the table window is hovered/focused)
    io = ig.GetIO()
    window_focused = ig.IsWindowFocused(ig.ImGuiFocusedFlags_RootAndChildWindows)
    window_hovered = ig.IsWindowHovered(ig.ImGuiHoveredFlags_RootAndChildWindows)
    active = window_focused || window_hovered

    selection_changed = false
    selected = state.selected_rows
    all_rows = 1:n_rows |> collect   # Vector{Int} for the selection helper

    if active && n_rows > 0
        ctrl_held = unsafe_load(io.KeyCtrl)
        shift_held = unsafe_load(io.KeyShift)

        # Ctrl/Cmd+A: select all
        if ig.IsKeyPressed(ig.ImGuiKey_A) && ctrl_held
            if length(selected) != n_rows
                state.selected_rows = collect(1:n_rows)
                selection_changed = true
            end
        end

        # Escape: clear selection
        if ig.IsKeyPressed(ig.ImGuiKey_Escape)
            if !isempty(selected)
                empty!(state.selected_rows)
                selection_changed = true
            end
        end

        # Arrow navigation
        anchor = isempty(selected) ? 0 : selected[end]
        if ig.IsKeyPressed(ig.ImGuiKey_DownArrow) && anchor < n_rows
            new_anchor = anchor + 1
            _update_multi_selection!(state.selected_rows, new_anchor, all_rows, shift_held, ctrl_held && !shift_held)
            state.scroll_to_row = new_anchor
            selection_changed = true
        end
        if ig.IsKeyPressed(ig.ImGuiKey_UpArrow) && anchor > 1
            new_anchor = anchor - 1
            _update_multi_selection!(state.selected_rows, new_anchor, all_rows, shift_held, ctrl_held && !shift_held)
            state.scroll_to_row = new_anchor
            selection_changed = true
        end
    end

    state.focused = window_focused

    # Programmatic scroll to row
    if state.scroll_to_row !== nothing
        target_row = state.scroll_to_row
        if 1 <= target_row <= n_rows
            ig.SetScrollY(Float32((target_row - 1) * ig.GetTextLineHeightWithSpacing()))
        end
        state.scroll_to_row = nothing
    end

    # Virtualized rows via ImGuiListClipper
    selected_set = Set(selected)
    clipper = ig.lib.ImGuiListClipper()
    try
        ig.Begin(clipper, n_rows, ig.GetTextLineHeightWithSpacing())
        while ig.Step(clipper)
            start_idx = Int(unsafe_load(clipper.DisplayStart)) + 1
            end_idx   = Int(unsafe_load(clipper.DisplayEnd))
            for row in start_idx:end_idx
                ig.PushID(row)
                ig.TableNextRow()

                # Row tint for provenance
                tint = row_tint(row)
                if tint !== nothing
                    ig.TableSetBgColor(ig.ImGuiTableBgTarget_RowBg0, tint)
                end

                ig.TableSetColumnIndex(0)
                is_selected = row in selected_set

                # Full-row selectable in first column
                if ig.Selectable(
                    "##row_$row",
                    is_selected,
                    ig.ImGuiSelectableFlags_SpanAllColumns,
                    (0.0f0, ig.GetTextLineHeightWithSpacing()),
                )
                    ctrl_held = unsafe_load(io.KeyCtrl)
                    shift_held = unsafe_load(io.KeyShift)
                    _update_multi_selection!(state.selected_rows, row, all_rows, shift_held, ctrl_held)
                    selected_set = Set(state.selected_rows)
                    selection_changed = true
                end

                # Overlay: first column text on same line as the selectable
                ig.SameLine()
                ig.TextUnformatted(cell(row, 1))

                # Remaining columns
                for col in 2:n_cols
                    ig.TableSetColumnIndex(col - 1)
                    ig.TextUnformatted(cell(row, col))
                end

                ig.PopID()
            end
        end
    finally
        ig.Destroy(clipper)
    end

    ig.EndTable()

    if selection_changed
        on_selection_change(copy(state.selected_rows))
    end

    return nothing
end
