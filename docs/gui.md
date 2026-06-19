# GUI Internals

> The CImGui frontend. This doc describes GUI state and interaction contracts, not a map of files.

## Entry point and frame loop

`open_browser(workspace)` creates one `BrowserState` and runs the render loop for an already-opened
workspace. The folder-open UI still uses the high-level callback project path internally:
`open_workspace(project, root_path)`. Docking layout is configured once at startup: left side for
navigation and information, right side for plot-oriented work.

## State boundary

Render functions receive one `BrowserState`. Its `workspace` field is the browser's single reference
to the open workspace. `PlotState` owns plot windows and choices, `TableInspectorState` owns the
generic file-inspection window, and `PerformanceState` owns samples used only for diagnostics.

The workspace owns the source, item index, selected collection and item identities, loaded cache,
loaded item memory, and the state of source-scan and cache work. Those values must not be copied
into browser state.

Source-scan and cache work each have independent state, progress, errors, cancellation, and task
ownership. This prevents one operation from overwriting the status shown for the other.

## Panels

| Panel | Role |
|---|---|
| Hierarchy tree | Multi-select tree, primary navigation. |
| Plot Area | Main plot window with plot-kind chooser, Live toggle, Detach, Export, and Help. |
| Information | Selected collection and item details. |
| Table Inspector | Materialized item-data viewer with per-row provenance, multi-select, and quick X/Y plot. |

### Selection flow (tree → plot)

1. A tree or future spatial-browser interaction writes collection paths and item keys to
   `workspace.selection`.
2. Each frame resolves those stable ids against `workspace.index` and applies tag visibility.
3. Plot windows with `Live` enabled use that visible selection and redraw when it changes.

A new panel changes selection only through the workspace. It does not maintain its own selected
objects.

## Plot rendering

Plots render directly when their inputs change. `PlotViewState` stores stable item keys,
plot type, figure, errors, and Live setting. It resolves those ids from the workspace index when
rendering, so browser state never keeps a second copy of item records. The main and detached
plots use the same type and the same rendering path.

Each plot window owns its plot kind and Live setting. The main plot starts with Live enabled, so it follows the browser selection. Detached plot windows start with Live disabled, so they keep the items they were created with unless the user enables Live in that window.

While the app runs, plot choices are stored as `item kind => plot type`. Project-local
persistence writes the type names to `measurementbrowser.toml`. When the main Live plot sees one
item kind, it restores that kind's last plot choice.

## Context menus (right-click)

Two existing popups, both attached to selected browser items:

- Tree nodes: "Mark Bad" / "Unmark Bad".
- Item rows: "Open Plot in New Window", mark/unmark bad.

Both apply to multi-selection via `_set_collections_bad!` / `_set_items_bad!`.

## DataGrid component

`src/Gui/DataGrid.jl` is a reusable, model-driven virtualized table widget. Any future panel can
use it by providing callbacks rather than passing DataFrames or item objects directly.

**State**: `DataGridState` (defined in `Browser/State.jl`) is owned by the consumer and passed by
reference each frame. It tracks `selected_rows::Vector{Int}`, a `scroll_to_row` request, and a
`focused` flag written back each frame.

**Render API**:
```julia
render_data_grid!(id, state;
    n_rows, columns, cell,         # cell(row, col) -> String
    row_tint = _ -> nothing,       # row -> Union{Nothing,UInt32} packed RGBA
    initial_widths = nothing,      # Vector{Float32} to restore persisted column widths
    on_selection_change = identity,
    height = 0.0f0)                # 0 = fill available
```

The grid uses `ImGuiListClipper` for O(visible) rendering, sticky header via
`TableSetupScrollFreeze(0, 1)`, `ImGuiTableFlags_Resizable` for drag-to-resize columns, and the
shared `_update_multi_selection!` helper for click/Shift/Ctrl/arrow/Escape/Ctrl+A selection.

**Column width persistence**: `TableGetColumnWidth` is not available in this CImGui build.
Provide `initial_widths` to restore saved widths on open; ImGui's resizable table state persists
within a session. To persist across sessions, record the widths externally before closing.

**Multi-select shared helper**: `_update_multi_selection!(selected, item, all_items, shift, ctrl)`
is defined once in `DataGrid.jl` and used by both `DataGrid` (row indices as `Int`) and the item
panel in `TreePanel.jl` (item ids as `String`).

## Table Inspector

The Table Inspector shows the **selected logical items' own data** (`item_data`) as a merged,
virtualized table. It replaces the former raw file re-parse approach.

**Primary mode (item data)**:
- On each frame, resolves the browser selection to `ItemRecord`s via `_project_visible_selection`
  and materializes items with `Workspace.materialize_items`.
- Merges multiple items' `DataFrame`s by column union; columns present in only some items render
  blank for the others. A stable key based on item ids prevents redundant rebuilds.
- Multi-item selections: subtle per-item alternating row background tint for provenance; an optional
  leading "\_item\_" column (toggled by the "Provenance column" checkbox) shows each row's source label.
- All rows are rendered (no cap); virtualization keeps rendering O(visible).
- Header is sticky (frozen first row).
- Multi-row selection with click / Shift+click / Ctrl+click / Ctrl+A / ↑↓ / Escape.
- Plot on the right (55/45 split): X/Y column chooser; "Plot selected only" checkbox restricts the
  plot to the grid selection.

**Per-kind column width persistence**: When a single-kind selection is open, the initial column widths
are restored from `BrowserState.table_column_widths_by_kind[kind]`. These are serialized to / from
`PersistedProjectView.table_column_widths` (comma-joined `Float32` widths keyed by `String(kind)`).
Mixed-kind selections skip width persistence and use auto-fit.

**Secondary raw-file mode**: `Inspect → Table Inspector` still exposes the path bar, `Live`
checkbox, `Open...`, and `Reload` controls for inspecting arbitrary delimited files. This mode
uses the legacy `inspect_table(path)` preview with the old non-virtualized table renderer (no
provenance). The `merge_item_tables` / `InspectorTable` model is item-data only.

`inspect_table(path)` remains exported for external use. It returns a `TablePreview` with the
file's detected delimiter, header row, first data row, approximate row count, and bounded preview
DataFrame. It does not call project readers or write cache entries.

## Annotations in the GUI

Tree nodes and item rows colour their text from the dominant effective tag for that item.
Untagged items render in the default colour. The `bad` tag's catalog entry is added on demand
whenever the user marks something bad. Visibility filtering uses the same effective-tag lookup over
the item key and its ancestors.

This is currently the browser's only annotation interface. Coordinates, spatial layout, and
notes are not yet shown or edited by the GUI.

## MakieIntegration constraint

Makie integration keeps **one Makie screen per `title_id`** in a global context. Multiple
independent canvases require either a unique `title_id` per canvas or a non-Makie renderer for the
secondary surface.
