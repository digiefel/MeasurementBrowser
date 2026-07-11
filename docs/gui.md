# GUI Internals

> The CImGui frontend. This doc describes GUI state and interaction contracts, not a map of files.

## Entry point and frame loop

`open_browser(workspace)` creates one `BrowserState` and runs the render loop for an already-opened
workspace. The first visible frames are a small preparation surface while the browser pays one-time
GLMakie/ImGui figure warmup runs through the registered `PlotsExtension`; the normal docked layout
is shown only after extensions report ready. The folder-open UI still uses the high-level callback project path internally:
`open_workspace(project, root_path)`. Docking layout is configured once at startup: left side for
navigation and information, right side for plot-oriented work.

## State boundary

Render functions receive one `BrowserState`. Its `workspace` field is the browser's single reference
to the open workspace. `TableInspectorState` owns the generic table-inspection window,
`PerformanceState` owns frame/UI samples used only for diagnostics, and optional GUI extensions
(such as `DataBrowserPlots.PlotsExtension`) own Makie plot windows, table plotting, and plot
persistence under `extensions.PlotsExtension` in `databrowser.toml`. Scan timings live with the
project performing the callbacks. Internal structured traces
and optional CPU samples belong to the workspace that captured them.

The workspace owns the source, item index, selected collection and item identities, DuckDB cache,
work graph, and source-scan/cache state. Those values must not be copied into browser state.

The GUI reads background-work state through one snapshot, `workspace.status` (a `WorkspaceStatus`):
its level (color), short label, one detail line, busy flag, optional progress fraction, cache/source
counts, and the live error list. `poll_workspace!` recomputes that snapshot only when something
changes or work is active, so idle frames are allocation-free, and the GUI never reaches into the
scan/analysis/cache internals directly. This is the stable contract: the cache's internals can change
without touching the UI as long as it still produces the same status. Source changes enter through
the startup scan or the source watcher; completed work folds progress, failures, and cache state into
the snapshot.

The cache popup shows one progress bar plus four fixed count cells: Sources, entries, analyzed, and
failures. Tooltips are multi-row and read only from `WorkspaceStatus`. Sources use the source-owned
noun, so `DirectorySource` says "source files". Entry, processed, analyzed, and failure counts come
from persisted cache tables. Item analysis can be cached before a processed payload is stored when
background processing is off, so the analyzed count may temporarily exceed processed; this is a
current engine wrinkle, not a GUI invariant.

Directory sources are watched continuously.

## Panels

| Panel | Role |
|---|---|
| Hierarchy tree | Multi-select tree, primary navigation. |
| Plot Area | Main plot window with plot-kind chooser, Live toggle, Detach, Export, and Help. |
| Information | Selected collection and item details. |
| Table Inspector | Materialized item-data viewer with per-row provenance and multi-select. |
| Table Plot (DBPlots) | Independent X/Y plot over the workspace selection's merged table; toggled from the Plot menu. |
| Performance | Frame/memory diagnostics, scan phase/source timings, plot timings, and opt-in internal profiling. |

## Scan profiling

The Performance window's always-on scan profile keeps one bounded row per source item. It shows
interpretation time, expanded item count, and the scheduler thread that performed the source work.
Processing and statistics belong to the separate processing/summarizing activities and are not
reported as source-read time.

Internal engine controls are absent unless the workspace was opened with
`profile_internal=true`. Its Internal tab starts and stops a bounded structured capture, groups span
latencies, separates writer wait from service, shows process counters and recent filtered events, and
exports Perfetto JSON. General live histories stay in the normal Live tab as CImGui sparklines
fed from bounded ring buffers (CSV export, no Makie). Starting during active work drains it and
starts one clean rebuild. With
`profile_cpu=true`, the same manual capture also retains a reduced Julia source-line hotspot report;
raw sampling buffers are cleared when capture stops. See [profiling.md](profiling.md).

## Development tools

The Debug menu exposes Dear ImGui's Metrics/Debugger, Debug Log, ID Stack Tool, debugger-gated Item
Picker, runtime debug options, style/font editor, user guide, and build information beside the
DataBrowser Performance window. An ImPlot submenu exposes its metrics, style/configuration tools,
user guide, and demo; its context is created lazily and destroyed before the Dear ImGui context. These
tools are available for diagnostics that need axes, zooming, or dense interaction without treating
ImPlot as the default for simple tables and bars.

The Dear ImGui and ImPlot demos are secondary reference entries after their diagnostic tools.
Tool-window state is transient and does not change `databrowser.toml`.

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

When a plot view no longer owns a visible figure — closed detached window, empty selection, invalid
plot kind, or draw failure — the browser destroys its embedded Makie screen instead of only dropping
the `Figure` reference. This keeps old scene graphs and OpenGL render resources from accumulating
while browsing.

One render materializes its selection once. The same processed item objects are passed to plot setup
and drawing; the resulting Makie figure owns the plotted values for that plot state. There is no
package-level object LRU and no separate debug-plot path.

Each plot window owns its plot kind and Live setting. The main plot starts with Live enabled, so it follows the browser selection. Detached plot windows start with Live disabled, so they keep the items they were created with unless the user enables Live in that window.

While the app runs, plot choices are stored as `item kind => plot type`. Project-local
persistence writes the type names to `databrowser.toml`. When the main Live plot sees one
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
reference each frame. Row mode tracks `selected_rows::Vector{Int}`. Cell mode tracks the anchor and
focus corners of one rectangular selection. Both modes share the scroll request and focus state.

**Render API**:
```julia
render_data_grid!(id, state;
    n_rows, columns, cell,         # cell(row, col) -> String
    cell_link = (_, _) -> nothing, # optional local path or URL for one cell
    row_tint = _ -> nothing,       # row -> Union{Nothing,UInt32} packed RGBA
    on_selection_change = identity,
    selection_mode = :rows,        # :rows or :cells
    height = 0.0f0)                # 0 = fill available
```

The grid uses `ImGuiListClipper` for O(visible) rendering, sticky header via
`TableSetupScrollFreeze(0, 1)`, `ImGuiTableFlags_Resizable` for drag-to-resize columns, and the
shared `_update_multi_selection!` helper for row selection. Cell mode supports click-drag and
Shift/arrow range extension, select-all, and spreadsheet-compatible TSV clipboard copy. The Table
Inspector deliberately stays in row mode; diagnostic tables use cell mode.

**Column width persistence**: Column widths are persisted across restarts via ImGui's ini file
(see below). `ImGuiTableFlags_SizingFixedFit` auto-sizes columns on first appearance; on
subsequent opens the ini file restores the saved widths. The `id` argument forms part of the
ini table key, so different callers with different ids get independent saved widths.

**Multi-select shared helper**: `_update_multi_selection!(selected, item, all_items, shift, ctrl)`
is defined once in `DataGrid.jl` and used by both `DataGrid` (row indices as `Int`) and the item
panel in `TreePanel.jl` (item ids as `String`).

## ImGui ini file

ImGui's ini file persists window sizes, positions, and — most importantly — `[Table]` column
widths. The browser re-enables ini persistence and points it at a stable per-machine path:

```
<DEPOT_PATH[1]>/databrowser/imgui.ini
```

(`homedir()` is used when `DEPOT_PATH` is empty.) The directory is created at startup if it does
not exist. The ini pointer is held in a module-level `const _IMGUI_INI_BYTES` in `Layout.jl` so
the byte buffer lives for the whole process (dear imgui stores the pointer, not a copy).

**DockBuilder is still the authority for docking layout.** `_setup_docking_layout!` uses a
one-shot Ref so it runs exactly once per process; the ini file's `[Docking]` section is
overwritten by that one-shot setup and does not take effect.

**Side-effect note**: enabling the ini also causes ImGui to save and restore window positions and
sizes (anything not controlled by DockBuilder). In practice this is desirable, but it means
window positions persist per machine rather than being reset each launch.

## Table Inspector

The Table Inspector shows data in two modes, both rendered through the same `DataGrid` component.

**Primary mode (item data)**:
- On each frame, resolves the browser selection to `ItemRecord`s via `_project_visible_selection`
  and materializes items with `Workspace.materialize_items`.
- Merges multiple items' Tables.jl-compatible data by column union; columns present in only some
  items render blank for the others. Non-tabular items are skipped with a per-item warning rather
  than failing the whole view. A stable key based on item ids prevents redundant rebuilds.
- Multi-item selections: subtle per-item alternating row background tint for provenance; an optional
  leading "\_item\_" column (toggled by the "Provenance column" checkbox) shows each row's source label.
- All rows are rendered (no cap); virtualization keeps rendering O(visible).
- Header is sticky (frozen first row).
- Multi-row selection with click / Shift+click / Ctrl+click / Ctrl+A / ↑↓ / Escape.
- The inspector only shows tables. Plotting columns lives in the separate Table Plot window
  (DBPlots), an independent visualizer over the same workspace selection.

**Column width persistence (item mode)**: The DataGrid is called with a per-kind `id` —
`string(kind)` for a single-kind selection, `"mixed"` when multiple kinds are selected. ImGui
keys `[Table]` entries in the ini file by that id, so column widths are stored and restored
per kind globally across restarts. There is no toml layer for column widths.

**Secondary raw-file mode**: `Inspect → Table Inspector` exposes the path bar, `Live` checkbox,
`Open...`, and `Reload` controls for inspecting arbitrary delimited files. The full file is
loaded (no row cap; a soft warning appears above 100 000 rows) and rendered through `DataGrid`
with id `"file"`. This mode has no provenance tinting. The grid model
is built by `_file_grid_model(preview)` which returns `(columns, n_rows, cell)` from the parsed
`TablePreview`; its column widths are persisted separately under the `"file"` ini key.

`inspect_table(path)` remains exported for external use. It returns a `TablePreview` with the
file's detected delimiter, header row, first data row, approximate row count, and a DataFrame of
all rows (or a bounded subset when called with an explicit `max_rows`). It does not call project
readers or write cache entries.

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

ImGui never samples the Makie framebuffer's color attachment directly. Each embedded figure owns
a plain 2D texture; after every Makie render the color attachment is copied into it and flushed
(`_sync_display_texture!`), and that copy is what reaches the draw list. Sampling an FBO-attached
texture from another GL context — a plot window dragged outside the main window becomes its own
platform viewport with multi-viewport enabled — makes the macOS driver log "GLD_TEXTURE_INDEX_2D
is unloadable" and substitute a zero texture; the plain copy is legal from every shared context,
so detached windows render on any monitor. The hidden warmup window additionally stays pinned to
the main viewport (`SetNextWindowViewport`) so startup never creates a platform window.
