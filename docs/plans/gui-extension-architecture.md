# GUI Extension Architecture

## Purpose

`DataBrowserGUI` is the host shell and base visualizer layer. `DataBrowserPlots`
is the default Makie extension package: it adds Makie embedding, high-level
plotting, enhanced default visualizers, and plot-related windows to that shell.

This replaces the earlier split where `DataBrowserGUI` depended on
`DataBrowserPlots`. That direction made plotting part of the shell and left no
real extension boundary to reuse for analysis packages.

## Target Package Graph

```text
DataBrowserAPI
DataBrowserProfiling
DataBrowserAnnotations -> API
DataBrowserSources     -> API
DataBrowserCache       -> API, Profiling, DataFrames/DuckDB
DataBrowserCore        -> API, Sources, Cache, Annotations, Profiling
DataBrowserGUI         -> API, Core, Annotations, Sources, CImGui/GLFW
DataBrowserPlots       -> API, Core, GUI, GLMakie
DataBrowser            -> GUI, default extension set
```

The important edge is `Plots -> GUI`, not `GUI -> Plots`.

`DataBrowserGUI` should load and run without GLMakie. It may expose empty plot
regions, menu slots, or extension hooks, but it must not import Makie code.
`DataBrowserPlots` depends on GUI and registers the Makie plotting experience.
The umbrella loads `DataBrowserPlots` by default for now because Makie plotting
is central to the current product, but the package graph should allow a base
GUI without it.

## Ownership

`DataBrowserGUI` owns:

- browser lifecycle and `BrowserState`;
- docking, layout, menus, status bars, and modal plumbing;
- shared CImGui widgets such as `DataGrid`;
- lightweight visualizers for Core-supported data shapes, such as table
  inspectors and later simple array/image inspectors that need no heavy plotting
  stack or data-library dependency;
- low-level user-defined GUI visualizer hooks;
- the window/menu/panel registry used by first-party and external GUI packages
  (extension instances own their own state).

`DataBrowserPlots` owns:

- `MakieIntegration.jl` and the GLMakie-in-ImGui bridge;
- plot panels, detached plot windows, plot warmup, export, and plot profiling UI;
- the Makie visualizer API currently exposed through `register_plot!`;
- default Makie visualizers: generic plotter, enhanced table plots, images,
  signals, arrays, metadata, nested data, and already-created Makie figures;
- a table-plot window: an independent Makie visualizer over the workspace
  selection, decoupled from the base table inspector;
- plot-specific persisted state, registered through the GUI state hooks.

`DataBrowserCore` owns materialization and workspace data access. Plot packages
ask Core for materialized items through public Core/API functions; they do not
make GUI load cache files or source files directly.

## GUI Registry

The smallest useful GUI extension surface is dispatch on an abstract type owned
by `DataBrowserGUI`. An extension is a mutable struct; the shell calls generic
functions with no-op fallbacks, so an extension implements only what it needs:

```julia
abstract type GuiExtension end

init!(ext::GuiExtension, state::BrowserState) = nothing
menu!(ext::GuiExtension, state::BrowserState) = nothing
draw!(ext::GuiExtension, state::BrowserState) = nothing
reset!(ext::GuiExtension, state::BrowserState) = nothing
shutdown!(ext::GuiExtension, state::BrowserState) = nothing
is_ready(ext::GuiExtension, state::BrowserState) = true
save_view(ext::GuiExtension, state::BrowserState) = Dict{String,Any}()
load_view!(ext::GuiExtension, state::BrowserState, view) = nothing

register_gui_extension!(::Type{<:GuiExtension})
```

The registry stores extension *types* in load order; each browser instantiates
one instance per type at startup. The per-browser mutable instance is the
extension's state — no separate keyed state slot. Each frame the shell draws
its own core panels, then calls `menu!` and `draw!` on every instance; the
remaining lifecycle methods cover startup, warmup gating, workspace switches,
shutdown, and the per-project saved view (persisted under a namespaced
`extensions.<id>` table).

The registry should stay boring: no dependency injection container, no plugin
discovery yet, no persistence machinery beyond the namespaced view table.

## Visualizer Layers

Core provides support for common data shapes without owning a GUI:

- tables (implemented now);
- vectors, arrays, and image-like values later, when the data contracts are
  clear. Shapes needing heavy dependencies (audio codecs, say) belong in
  extension packages, not Core.

GUI provides simple visualizers for those shapes using CImGui and native widgets:

- table inspection through `DataGrid`;
- lightweight matrix/image inspection when those Core shapes exist;
- low-level custom visualizer callbacks for packages that want to draw directly
  in the GUI shell.

Plots provides Makie-based visualizers on top:

- the current registered plot API;
- generic high-level plot builders;
- richer defaults for tables, arrays, images, signals, and domain packages.

This means the table inspector stays in `DataBrowserGUI`; its connected Makie
plotter belongs in `DataBrowserPlots`.

## DataBrowserPlots Load Path

`DataBrowserPlots.__init__()` registers its GUI callbacks if `DataBrowserGUI` is
loaded because `DataBrowserPlots` directly depends on `DataBrowserGUI`.

The umbrella package imports both:

```julia
using DataBrowserGUI
using DataBrowserPlots
```

That gives the normal app the plotting extension by default. A headless or
minimal GUI environment can load `DataBrowserGUI` without loading
`DataBrowserPlots`. Later, the umbrella can make the default extension set
configurable without changing the package boundary.

## Migration Plan

The file-level execution plan, phased to keep tests green, is
[gui-plots-inversion.md](gui-plots-inversion.md).

1. Correct the package docs and `Project.toml` direction:
   `DataBrowserGUI` stops depending on `DataBrowserPlots`; `DataBrowserPlots`
   depends on `DataBrowserGUI`.
2. Remove GLMakie from GUI diagnostics by following
   [diagnostics.md](diagnostics.md): replace Makie performance plots with CImGui
   tables, bars, and sparklines.
3. Add the tiny GUI registry to `DataBrowserGUI`.
4. Move `MakieIntegration.jl`, plot state, `PlotPanel.jl`, and plot warmup from
   GUI into `DataBrowserPlots`.
5. Move `InspectorTable` to Core and make its construction Tables.jl-generic.
   The GUI table inspector reads the Tables.jl interface rather than depending
   on DataFrames; no wrapper helpers — dispatch is the contract.
6. Keep TableInspector and shared tabular inspection in `DataBrowserGUI`; drop
   its embedded quick-plot. `DataBrowserPlots` provides a separate table-plot
   window: an independent visualizer over the workspace selection that builds
   its own table from the selected items and plots typed values, never display
   strings. The two windows share only the Core table model; when both are
   live they agree because they follow the same workspace selection.
7. Make the umbrella load `DataBrowserGUI` and then the default GUI extension
   set, starting with `DataBrowserPlots`.
8. Only after the first-party path works, decide whether any optional glue
   belongs in Julia `ext/`.

## Checks

- `using DataBrowserGUI` does not load GLMakie.
- `using DataBrowserPlots` loads GLMakie and registers plot windows with GUI.
- `using DataBrowser` gives the full default app with plots.
- `DataBrowserCore` tests run headless with no GUI or GLMakie imports.
- One built-in plotting window is registered through the same registry an
  external package would use.
- The base table inspector opens without DBPlots; the table-plot window exists
  only when DBPlots is loaded.
