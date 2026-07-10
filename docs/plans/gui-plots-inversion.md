# GUI/Plots Inversion — Execution Plan

**Status (2026-07-11): executed.** All phases landed on `dbplots`, followed by a review pass
fixing eight findings (inspector crash, sparkline render, typed cell access, table caching,
warmup ordering, shell/extension boundary hooks, docs). Kept for reference; current behavior
is documented in [../gui.md](../gui.md).

The concrete, file-level plan for the migration in
[gui-extension-architecture.md](gui-extension-architecture.md): flip the dependency edge so
`DataBrowserPlots -> DataBrowserGUI`, leaving a GUI shell that loads without GLMakie. Written to be
executed phase by phase; **the full test suite must pass and the app must be committable after every
phase**. Run tests with `julia --project --threads=4 -e 'using Pkg; Pkg.test()'`.

Out of scope (do not touch): the Cache edge flip (`DataBrowserCache` still depends on Core — separate
PR, see roadmap), the full diagnostics overhaul in [diagnostics.md](diagnostics.md) (Perfetto removal,
diagnostic runs — only the Performance *window's* Makie usage is in scope here), any generic plotter
work, and the table inspector's raw file-preview mode (`preview::TabularFileSource` in
`TableInspectorState`) — its planned removal into proper tabular data sources is separate work; leave
it in GUI as is.

## Current state (verified 2026-07-10)

- `lib/DataBrowserGUI/Project.toml` depends on `DataBrowserPlots`, `GLMakie`, `DataFrames`.
- `lib/DataBrowserPlots` is small: `DataBrowserPlots.jl` (the `setup_plot`/`plot_data!`
  materialize-and-dispatch bridge) and `table_model.jl` (`InspectorTable`, `merge_item_tables`).
- GLMakie appears in GUI at: `Browser/State.jl` (`PlotViewState.figure`, `PlotState.warmup_figure`,
  `LivePlotsState` figures/Observables), `Browser/MakieIntegration.jl` (the whole
  `MakieImguiIntegration` module), `Gui/PlotPanel.jl`, `Gui/PerformanceWindow.jl` (two live Makie
  figures + `Makie.save` export), `Gui/TableInspector.jl` (quick X/Y plot).
- The render loop (`Gui/Layout.jl`, `_run_browser`) hardwires `render_plot_window`,
  `render_additional_plot_windows`, and a warmup gate on `state.plots.runtime_warmed`.
- `src/Precompile.jl` exercises `Browser._make_timings_figure` / `_make_build_figure` /
  `_sample_build!` and the plot pipeline.
- Tests touching moved code: `test_registry_plot.jl`, `test_type_api.jl`, `test_tase_analysis.jl`,
  `test_table_inspector.jl`, `test_event_driven_workspace.jl`, `test_performance_window.jl`,
  `test_project_view_state.jl`, `test_glfw_scan_startup_stress.jl`.

## Phase 0 — Table model moves to Core

The typed table model is Core's "base support for tables" (vision §5/§6); the inspector *widget*
stays in GUI; its Makie quick-plot moves to Plots in Phase 3.

1. Move `lib/DataBrowserPlots/src/table_model.jl` to `lib/DataBrowserCore/src/TableModel.jl`. Keep
   the code as is; carry its imports (`DataFrames`, item/label helpers) — Core already depends on
   DataFrames and DataBrowserAPI. Include it from `DataBrowserCore.jl` and export `InspectorTable`,
   `merge_item_tables` from `DataBrowserCore`.
2. Update consumers: `lib/DataBrowserGUI/src/Browser/State.jl:9` and
   `lib/DataBrowserGUI/src/Gui/TableInspector.jl:8` import from `DataBrowserCore` instead of
   `DataBrowserPlots`. Remove the include and re-export from `DataBrowserPlots.jl`.
3. Make the table path Tables.jl-generic instead of DataFrames-specific. `InspectorTable`
   construction accepts any Tables.jl table (`Tables.columns`, column names via
   `Tables.columnnames`), and the payload checks in `Gui/TableInspector.jl` (~lines 100, 116) use
   the Tables.jl interface instead of `isa DataFrame`. Do not add wrapper helpers — multiple
   dispatch and the Tables.jl API are the contract. Add the lightweight `Tables` dep to GUI now;
   `DataFrames` leaves GUI's deps in Phase 4.
4. Fix test imports (`test_table_inspector.jl` and any other test importing these names from
   `DataBrowserPlots`). Run the suite.

## Phase 1 — Performance window loses Makie

Per [diagnostics.md](diagnostics.md): CImGui-native sparklines, no GLMakie in GUI diagnostics.

1. In `Browser/State.jl`, rewrite `LivePlotsState`: delete `timings_figure`, `build_figure`, and all
   `Observable` fields. Keep plain `Vector{Float32}` ring buffers for every series (the build series
   already have `_buf` vectors; give the four timing series x/y ring buffers too). Keep `capacity`
   and the finite-difference sampling fields.
2. In `Gui/PerformanceWindow.jl`: remove the GLMakie imports and the `MakieImguiIntegration` use.
   Delete `_make_timings_figure` and `_make_build_figure`. Rewrite `_update_live_timings!` and
   `_sample_build!` to fill the plain ring buffers. Render each series with `ig.PlotLines` (one
   sparkline per series, current value as overlay text) inside the existing tab layout. Replace the
   `Makie.save` PNG export with a CSV export of the ring buffers (same button), or drop the button —
   CSV preferred, matching diagnostics.md.
3. Update `src/Precompile.jl:92-94` to exercise the new update/render helpers instead of the figure
   constructors. Update `test_performance_window.jl`.
4. `Browser/State.jl` still imports GLMakie for `PlotViewState`/`warmup_figure` — that goes in
   Phase 3. Run the suite.

## Phase 2 — GUI extension registry (additive; nothing moves yet)

1. New file `lib/DataBrowserGUI/src/Browser/Extensions.jl`, included early in `Browser.jl`. The
   registry is dispatch-based: an extension is a mutable struct subtyping an abstract type, and the
   shell calls generic functions that have no-op fallbacks, so an extension implements only what it
   needs:

   ```julia
   abstract type GuiExtension end

   extension_id(ext::GuiExtension)::String = String(nameof(typeof(ext)))
   init!(ext::GuiExtension, state::BrowserState) = nothing
   menu!(ext::GuiExtension, state::BrowserState) = nothing
   draw!(ext::GuiExtension, state::BrowserState) = nothing
   reset!(ext::GuiExtension, state::BrowserState) = nothing       # workspace switch/close
   shutdown!(ext::GuiExtension, state::BrowserState) = nothing
   is_ready(ext::GuiExtension, state::BrowserState)::Bool = true  # warmup gating
   save_view(ext::GuiExtension, state::BrowserState)::Dict{String,Any} = Dict{String,Any}()
   load_view!(ext::GuiExtension, state::BrowserState, view::Dict{String,Any}) = nothing

   register_gui_extension!(::Type{<:GuiExtension})  # appends the type in load order; same type replaces
   ```

2. `BrowserState` gains `extensions::Vector{GuiExtension}`; `open_browser` instantiates each
   registered type with its zero-argument constructor. The instance is per-browser and mutable —
   it *is* the extension's state. There is no separate keyed state slot; extension packages put
   their state in their own struct fields.
3. Wire the shell (`Gui/Layout.jl`): instantiate and `init!` in `open_browser` after `BrowserState`
   construction; call each `menu!` at the end of `render_menu_bar`; call each `draw!` after the
   core panels in the frame loop; call `shutdown!` in the `on_exit` handler; call `reset!` where
   `Operations.jl:68` currently does `state.plots = PlotState()`.
4. Persistence: add `extensions::Dict{String,Dict{String,Any}} = Dict()` to `PersistedProjectView`,
   keyed by `extension_id`; on save, collect each extension's `save_view`; on load, dispatch
   `load_view!`. Unknown ids are kept and rewritten as-is (an extension may not be loaded this
   session).
5. Warmup gate: replace the `!state.plots.runtime_warmed` gate at `Gui/Layout.jl:587` with
   "startup preparation until `all(is_ready(ext, state) for ext in state.extensions)`", keeping the
   existing `_render_startup_preparation!` frames. Until Phase 3 the plot warmup still runs in the
   shell, so keep the old call temporarily under the new gate.
6. Add a test: register a dummy extension type, assert instantiation and callback order, instance
   state surviving across frames, and `save_view`/`load_view!` roundtrip through
   `PersistedProjectView`. Run the suite.

## Phase 3 — Makie code moves to DataBrowserPlots

The big move. `DataBrowserPlots` temporarily depends on both GUI and its old direct deps; GUI's
`Project.toml` is cleaned in Phase 4. Move, don't rewrite — code should land verbatim plus import
fixes unless a step says otherwise.

1. Move `lib/DataBrowserGUI/src/Browser/MakieIntegration.jl` (module `MakieImguiIntegration`) to
   `lib/DataBrowserPlots/src/MakieIntegration.jl`.
2. Move `PlotViewState` and `PlotState` from `Browser/State.jl` to a new
   `lib/DataBrowserPlots/src/PlotState.jl`, held by `mutable struct PlotsExtension <: GuiExtension`
   (either as fields or as one `plots::PlotState` field — keep `PlotState` as a type, it is
   referenced widely). Delete the `plots::PlotState` field from `BrowserState`; plot code receives
   the extension instance through the `GuiExtension` methods.
3. Move `lib/DataBrowserGUI/src/Gui/PlotPanel.jl` to `lib/DataBrowserPlots/src/PlotPanel.jl`
   (`render_plot_window`, `render_additional_plot_windows`, export, `profile_next_plot` console
   breakdown, `ensure_plot_runtime_warmed!`). It needs `CImGui`, `NativeFileDialog`, `Printf` in
   Plots' deps. `profile_next_plot::Bool` stays on `BrowserState` only if other panels set it —
   otherwise move it into `PlotState`.
4. Move the plot parts of `Browser/Persistence.jl` (`PersistedPlotView`, `_persisted_plot_view`,
   `plot_kinds` save/restore) into Plots, expressed as its `save_view`/`load_view` callbacks writing
   plain strings under `extensions.PlotsExtension` in `databrowser.toml`. No migration and no
   compatibility: old plot tables in existing `databrowser.toml` files are ignored, and a file
   without the new table starts fresh.
5. Decouple table plotting from the table inspector completely. Delete the quick-plot state from
   `TableInspectorState` (`x_column`, `y_column`, `figure`, `plot_key`, `plot_error`,
   `plot_selected_only`) and the Makie render block (~`TableInspector.jl:337`); the inspector only
   shows tables. `DataBrowserPlots` provides its own table-plot window, drawn through its `draw!`
   like its other windows, as an independent visualizer over the *workspace* selection: it
   materializes the selected items and builds its own `InspectorTable` from them, with the same
   live/pinned behavior as other plot windows. No accessors into the inspector, no shared GUI
   state — when both windows are live they show the same items because they follow the same
   workspace selection. The two share only the Core table model, and the plotter reads typed
   values, never display strings.
6. Remove the plot menu items and the two plot-window render calls from `Gui/Layout.jl`; Plots
   provides them via its `menu!`/`draw!` methods. Plots' `is_ready` wraps
   `ensure_plot_runtime_warmed!` (returns `runtime_warmed`); delete the shell's temporary warmup
   call from Phase 2 step 5. Its `reset!` reinitializes the instance's plot state.
7. `DataBrowserPlots.__init__` calls `register_gui_extension!(PlotsExtension)`. Keep the existing
   `setup_plot`/`plot_data!` bridge methods as they are.
8. Fix `Browser.jl` (drop `using DataBrowserPlots` and moved includes), `src/Precompile.jl`
   (moved names now live in `DataBrowserPlots`), and the affected tests
   (`test_project_view_state.jl`, `test_glfw_scan_startup_stress.jl`, plot tests). Run the suite.

## Phase 4 — Flip the Project.toml edges

1. `lib/DataBrowserGUI/Project.toml`: remove `DataBrowserPlots`, `GLMakie`, and `DataFrames` from
   `[deps]`/`[compat]`; keep the lightweight `Tables` dep added in Phase 0. If anything still fails
   to compile without the removed deps, it was missed in an earlier phase — move it, don't re-add
   the dep.
2. `lib/DataBrowserPlots/Project.toml`: add `DataBrowserGUI`, `CImGui`, `NativeFileDialog`
   (+ `Printf`). Root `Project.toml` `[sources]` already has both packages.
3. `src/DataBrowser.jl` keeps `using DataBrowserGUI` and `using DataBrowserPlots`; Julia loads GUI
   first because Plots depends on it, and Plots' `__init__` registers the extension.
4. Run the suite.

## Phase 5 — Acceptance

- `GLMakie` is absent from `lib/DataBrowserGUI`'s resolved Manifest. The Phase 4 dep removal is the
  enforcement; the Manifest check is the whole verification.
- Full test suite green.
- Launch the app on a `test/fixtures` project through the umbrella: main plot follows selection,
  detached plot windows open/close, plot export works, the table inspector opens and the table-plot
  window plots columns of the selected items, the Performance window renders sparklines with no
  Makie, and the warmup screen still shows at startup.
- Behavioral boundary tests (add if missing): with only `DataBrowserGUI` loaded, the table inspector
  opens and no plot menu exists; with `DataBrowserPlots` loaded, the plot menu and the table-plot
  window appear.
- Update docs in the final commit: `docs/gui.md` (new ownership), `CLAUDE.md` table rows for
  `lib/DataBrowserPlots/` / `lib/DataBrowserGUI/`, and mark the migration steps done in
  [gui-extension-architecture.md](gui-extension-architecture.md).

## Decisions already made (do not relitigate)

- `InspectorTable` lives in **Core** (base table support), the inspector widget in GUI, table
  plotting in Plots as its own window.
- The registry is dispatch-based: `abstract type GuiExtension` with no-op fallback methods;
  per-browser mutable instances own their own state. No keyed state slots, no plugin discovery.
- The tabular contract in GUI is the Tables.jl interface — no DataFrames types, no wrapper helpers.
- The table inspector and the table plotter are fully decoupled: both are independent visualizers
  over the workspace selection, sharing only the Core table model. Neither reads the other's state.
- Persisted plot view moves to Plots under an `extensions.<id>` TOML table. No migration, no
  compatibility: old plot tables in `databrowser.toml` are ignored.
- The Performance window keeps its tables and gains `ig.PlotLines` sparklines; PNG export becomes
  CSV export. Perfetto/trace removal is a later, separate change.
- No `ext/` package extensions in this pass; `DataBrowserPlots` is a hard dep of the umbrella only.
