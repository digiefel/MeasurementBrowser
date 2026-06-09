# GUI Internals

> The CImGui frontend. The code lives under [src/Gui/](../src/Gui/), included from the slim [src/Gui.jl](../src/Gui.jl) shell. This doc tells you where to land for any given task.

## File map

| File | Owns |
|---|---|
| [src/Gui.jl](../src/Gui.jl) | Imports (GLFW, GLMakie, CImGui, ModernGL, TOML, NativeFileDialog) and `MakieIntegration`; `include`s the files below into the parent `MeasurementBrowser` module. |
| [src/Gui/State.jl](../src/Gui/State.jl) | `ui_state` plumbing: prefs (`_prefs_path`, `_load_prefs`, `_save_prefs`, recents handling, `_open_project_path!`), per-domain init helpers (`_init_scan_state!`, `_init_cache_state!`, `_init_tag_state!`, `_init_plot_state!`, `_init_figure_script_state!`), buffer helpers (`_text_buffer`, `_buffer_string`, `_set_buffer_string!`), figure-script state machine (group ops, job request/launch/poll, fact index), scan/cache state machine (`_begin_scan!`, snapshot appliers, `_launch_source_scan_job!`, `_launch_project_reload_job!`, `_launch_cache_update_job!`, `_poll_scan_events!`, `_poll_source_scan_events!`), perf utilities (`_time!`, `_memory_snapshot`, `_print_perf_summary`), and `_helpmarker`. |
| [src/Gui/BadAndStyling.jl](../src/Gui/BadAndStyling.jl) | Tag-driven visibility and styling: `_load_tag_state_for_root!`, `_tag_state_ready`, `_tag_state_or_error`, `_show_bad_effective`, `_ancestor_keys_for_path`, `_has_bad_tag`, `_device_is_visible`, `_measurement_is_visible`, `_apply_visible_selection!`, `_set_devices_bad!`, `_set_measurements_bad!`, `_ensure_bad_catalog_entry!`, `_selection_targets`, `_render_tag_state_error!`, `_push_tag_text_style!`. |
| [src/Gui/PlotPanel.jl](../src/Gui/PlotPanel.jl) | Plot job queue (`_plot_job`, `_launch_plot_job!`, `_queue_plot_job!`, `_finish_plot_job!`, `_poll_plot_events!`, `_plot_target_loading`), plot toolbar helpers, main plot panel (`render_plot_window`), detached plot windows (`render_additional_plot_windows`), plot indicator/error helpers, runtime warmup. |
| [src/Gui/TreePanel.jl](../src/Gui/TreePanel.jl) | Hierarchy tree (`_render_hierarchy_tree_panel`), measurement panel (`_render_measurements_panel`), `render_selection_window`, multi-select utility (`_update_multi_selection!`), filter helpers (`_visible_measurements`, `_measurement_matches_filter`, …), `_render_selection_status!`, `_render_selection_toolbar!`, leaf-collection helpers. |
| [src/Gui/InfoModal.jl](../src/Gui/InfoModal.jl) | Information panel (`render_info_window`), figure-script window (`render_figure_script_window`, `_write_figure_script_from_ui!`, `_figure_script_output_path`, `_render_figure_script_group_tooltip`), missing-metadata modal (`render_device_info_modal`). |
| [src/Gui/Layout.jl](../src/Gui/Layout.jl) | Frame loop (`create_window_and_run_loop`, `start_browser`), docking (`_setup_docking_layout!`), top menu bar (`render_menu_bar`), project window (`render_project_window`), perf window (`render_perf_window`), cache toolbar / progress / controls (`_cache_activity_model`, `_cache_progress_models`, `_source_progress_models`, `_render_cache_toolbar_button!`, `_render_cache_controls!`, `_render_progress_indicator!`). |

## Entry point and frame loop

`create_window_and_run_loop` ([Layout.jl](../src/Gui/Layout.jl)) initializes an empty `ui_state::Dict{Symbol,Any}` and runs the render loop. Docking layout is configured once at startup in `_setup_docking_layout!` ([Layout.jl](../src/Gui/Layout.jl)) — a 2-column split: left (Hierarchy, Info), right (Plot Area).

## State: a single Dict

There is **no AppState struct**. Everything lives in `ui_state`, threaded by reference into every render function. Initialization is split across:

- `_init_scan_state!` (scan/progress)
- `_init_cache_state!` (HDF5 cache)
- `_init_tag_state!` (Tags)
- `_init_figure_script_state!` (figure-script export)
- `_init_plot_state!` (plot job queue)

All defined near the top of [State.jl](../src/Gui/State.jl).

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
| `:open_plot_windows` | `Vector{Dict{Symbol,Any}}` | Detached plot windows. Each stores its own plot kind, measurements, Live flag, figure, and error. |
| `:tag_state` | `Annotations.Tags.TagState` | Catalog + per-key assignments (loaded at project init). |

## Panels

| Panel | Render fn | Role |
|---|---|---|
| Hierarchy tree | `_render_hierarchy_tree_panel` ([TreePanel.jl](../src/Gui/TreePanel.jl)) | Multi-select tree, primary navigation. |
| Plot Area | embedded `MakieFigure` ([PlotPanel.jl](../src/Gui/PlotPanel.jl)) | Main plot window with plot-kind chooser, Live toggle, Detach, Export, and Help. |
| Information | `render_info_window`, `render_device_info_modal`, `render_figure_script_window` ([InfoModal.jl](../src/Gui/InfoModal.jl)) | Device modal and figure scripts. |

### Selection flow (tree → plot)

1. Click on tree node → `ig.IsItemClicked()` ([TreePanel.jl:178](../src/Gui/TreePanel.jl)).
2. Modifier check via `ig.GetIO()` ([TreePanel.jl:180](../src/Gui/TreePanel.jl)).
3. `_update_multi_selection!` ([TreePanel.jl:285](../src/Gui/TreePanel.jl)) writes to `:selected_device_paths`.
4. `_apply_visible_selection!` ([BadAndStyling.jl:125](../src/Gui/BadAndStyling.jl)) filters measurements and updates `:selected_measurements`.
5. Plot windows with `Live` enabled pick up the change next frame and queue a `PlotJob`.

**Any new panel that drives selection must write to `:selected_device_paths` and call `_apply_visible_selection!`** — that's how plots stay in sync.

## Async plot rendering

Plots render off the UI thread via the `PlotJob` queue ([PlotJobs.jl](../src/PlotJobs.jl) + [PlotPanel.jl](../src/Gui/PlotPanel.jl) `_queue_plot_job!`). Results are polled each frame in `_poll_plot_events!` ([PlotPanel.jl](../src/Gui/PlotPanel.jl)). The main figure is stored in `ui_state[:plot_figure]`; detached windows store figures in their own entry.

Each plot window owns its plot kind and Live setting. The main plot starts with Live enabled, so it follows the browser selection. Detached plot windows start with Live disabled, so they keep the measurements they were created with unless the user enables Live in that window.

## Context menus (right-click)

Two existing popups, both via `ig.BeginPopupContextItem`:

- Tree nodes ([TreePanel.jl:196](../src/Gui/TreePanel.jl)): "Mark Bad" / "Unmark Bad".
- Measurement rows ([TreePanel.jl:485](../src/Gui/TreePanel.jl)): "Open Plot in New Window", mark/unmark bad.

Both apply to multi-selection via `_set_devices_bad!` / `_set_measurements_bad!`.

## Tag-driven styling

Tree nodes and measurement rows colour their text via `_push_tag_text_style!` ([BadAndStyling.jl:228](../src/Gui/BadAndStyling.jl)), which resolves the dominant tag for a given key + ancestor list using `Annotations.Tags.dominant_color`. Untagged items render in the default colour. The `bad` tag's catalog entry (red, priority 100) is added on demand by `_ensure_bad_catalog_entry!` whenever the user marks something bad. Visibility filtering happens in `_device_is_visible` ([BadAndStyling.jl:86](../src/Gui/BadAndStyling.jl)) and `_measurement_is_visible` ([BadAndStyling.jl:91](../src/Gui/BadAndStyling.jl)); both compose `_has_bad_tag` over the key and its ancestors via `_ancestor_keys_for_path`.

## MakieIntegration constraint

[MakieIntegration.jl](../src/MakieIntegration.jl) (`MakieFigure` ~line 162) keeps **one Makie screen per `title_id`** in a global `makie_context` Dict. Multiple independent canvases (e.g., a second plot panel) require either a unique `title_id` per canvas or a non-Makie renderer for the secondary surface.

## Where to look

| Concern | Approx. location |
|---|---|
| Frame loop | [Layout.jl:817](../src/Gui/Layout.jl) `create_window_and_run_loop` |
| Docking layout | [Layout.jl:786](../src/Gui/Layout.jl) `_setup_docking_layout!` |
| State init | top of [State.jl](../src/Gui/State.jl) |
| Tree panel | [TreePanel.jl:36](../src/Gui/TreePanel.jl) |
| Plot embedding | [PlotPanel.jl](../src/Gui/PlotPanel.jl) |
| Plot job queue | [PlotPanel.jl](../src/Gui/PlotPanel.jl) |
| Tag styling + visibility | [BadAndStyling.jl:86, 91, 125, 228](../src/Gui/BadAndStyling.jl), [TreePanel.jl:196, 485](../src/Gui/TreePanel.jl) |
| Project + tag-state init | [BadAndStyling.jl:1](../src/Gui/BadAndStyling.jl) `_load_tag_state_for_root!` |
