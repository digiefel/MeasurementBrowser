# GUI Internals

> The CImGui frontend. This doc describes GUI state and interaction contracts, not a map of files.

## Entry point and frame loop

`create_window_and_run_loop` initializes an empty `ui_state::Dict{Symbol,Any}` and runs the render
loop. Docking layout is configured once at startup: left side for navigation and information, right
side for plot-oriented work.

## Current State Shape

The GUI currently passes a single `ui_state::Dict{Symbol,Any}` by reference through render
functions. That is an implementation fact, not a long-term architecture decision. App state and
persistence are active design questions, especially as the browser moves toward GUI/API parity,
saved workflows, and non-blocking REPL use.

Until that is redesigned, treat `ui_state` as transient GUI state. Avoid adding durable project,
cache, workflow, or figure state there unless there is no clearer owner yet. Initialization is split
by domain:

- `_init_scan_state!` (scan/progress)
- `_init_cache_state!` (HDF5 cache)
- `_init_tag_state!` (Tags)
- `_init_figure_script_state!` (figure-script export)
- `_init_plot_state!` (plot windows)

### Key `ui_state` entries you'll touch most

| Key | Type | Purpose |
|---|---|---|
| `:selected_device_paths` | `Vector{String}` | Selected paths (canonical, source-of-truth). |
| `:selected_devices` | `Vector{HierarchyNode}` | Resolved nodes for the above. |
| `:selected_measurements` | `Vector{MeasurementInfo}` | Filtered down by `_apply_visible_selection!`. |
| `:plot_figure` | Makie `Figure` or `nothing` | Current single-plot figure. |
| `:main_plot_kind` | `Type{<:PlotKind}` or missing | Plot kind selected in the main plot toolbar. |
| `:main_plot_live` | `Bool` | Whether the main plot follows the current browser selection. Defaults to `true`. |
| `:main_plot_measurements` | `Vector{MeasurementInfo}` | Measurements used by the main plot when Live is disabled. |
| `:plot_kind_preferences` | `Dict{String,Dict{String,String}}` | Remembered plot kind by project and measurement kind. |
| `:open_plot_windows` | `Vector{Dict{Symbol,Any}}` | Detached plot windows. Each stores its own plot kind, measurements, Live flag, figure, and error. |
| `:tag_state` | `Annotations.Tags.TagState` | Catalog + per-key assignments (loaded at project init). |

## Panels

| Panel | Role |
|---|---|
| Hierarchy tree | Multi-select tree, primary navigation. |
| Plot Area | Main plot window with plot-kind chooser, Live toggle, Detach, Export, and Help. |
| Information | Device modal and figure scripts. |

### Selection flow (tree → plot)

1. A tree or future spatial-browser interaction writes canonical paths to `:selected_device_paths`.
2. The visible-selection pass resolves paths, applies tag-driven visibility, and updates
   `:selected_measurements`.
3. Plot windows with `Live` enabled pick up the change next frame and redraw if their plot key changed.

**Any new panel that drives selection must write to `:selected_device_paths` and call `_apply_visible_selection!`** — that's how plots stay in sync.

## Plot rendering

Plots render directly when their plot key changes. The main figure is stored in
`ui_state[:plot_figure]`; detached windows store figures in their own entry.

Each plot window owns its plot kind and Live setting. The main plot starts with Live enabled, so it follows the browser selection. Detached plot windows start with Live disabled, so they keep the measurements they were created with unless the user enables Live in that window.

Plot kind choices are persisted in the app preferences as project name → measurement kind → plot kind type name. When the main Live plot sees a single measurement kind, it uses that remembered plot kind or clears the chooser if none has been chosen for that kind yet.

## Context menus (right-click)

Two existing popups, both attached to selected browser items:

- Tree nodes: "Mark Bad" / "Unmark Bad".
- Measurement rows: "Open Plot in New Window", mark/unmark bad.

Both apply to multi-selection via `_set_devices_bad!` / `_set_measurements_bad!`.

## Tag-driven styling

Tree nodes and measurement rows colour their text from the dominant effective tag for that item.
Untagged items render in the default colour. The `bad` tag's catalog entry is added on demand
whenever the user marks something bad. Visibility filtering uses the same effective-tag lookup over
the item key and its ancestors.

## MakieIntegration constraint

Makie integration keeps **one Makie screen per `title_id`** in a global context. Multiple
independent canvases require either a unique `title_id` per canvas or a non-Makie renderer for the
secondary surface.
