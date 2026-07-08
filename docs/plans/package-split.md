# Package Split — Execution Plan

The concrete steps behind [roadmap.md](roadmap.md) Stage 1: how the single `MeasurementBrowser`
package becomes the `DataBrowser` family in [../vision.md](../vision.md) §6. This doc owns the
file-by-file mapping and the ordering; the roadmap owns *why* and *when*, the vision owns the target
boundaries.

## Method: extract one named package at a time, bottom-up

The split is a sequence of extractions, not one big move. At each step we lift the code for one
function out of the main package into `lib/DataBrowser<Name>`, give it its final name and its own
`Project.toml`, and wire it back with a `[sources]` path dependency so the main package keeps
resolving it without a registry entry. After each step the whole suite still builds and passes.

Two properties make this work:

- **Each step is independently reviewable.** The unit of review is a real, finally-named package with
  a real boundary — "is this what `DataBrowserCache` should contain, and are its deps right?" — not a
  reshuffle of `MeasurementBrowser` internals under names that no longer mean anything.
- **The rename falls out for free.** We never do a separate "rename" pass. Everything is extracted
  into `DataBrowser*` packages; what remains in the original package shrinks until it is just the thin
  umbrella, and the final step renames that residual `MeasurementBrowser` → `DataBrowser`.

Extraction goes bottom-up (leaves first) so a package only ever depends on packages already extracted.

## Where the code is today

The API-shaped code already lives largely apart from the engine, which is what makes step 2 tractable:

| Current file / dir | Holds | Lands in |
|---|---|---|
| `src/Profiling.jl` | instrumentation, trace spans | `DataBrowserProfiling` |
| `src/Projects.jl` | abstract extension types, value types (`Project`, recipe types), **generic-function declarations** — imports nothing today | `DataBrowserAPI` |
| `src/Project.jl` | `define_project` / `register_*` **plus** the project's engine-interface *methods*, serialization, and `using DataFrames` | split: construction → API, the rest → Core |
| `src/Annotations/` (nested pkg) | tags / notes / coords / layout, keyed on `String` paths | `DataBrowserAnnotations` |
| `src/DataSources/DirectorySource.jl` | the concrete directory source | `DataBrowserSources` |
| `src/TableInspector.jl` — `TablePreview` / `inspect_table` | raw delimited-**file** preview (detect delimiter/header) | `DataBrowserSources` |
| `src/TableInspector.jl` — `InspectorTable` / `merge_item_tables` | merge item `.data` frames for the table **viewer** | `DataBrowserPlots` |
| `src/Gui/TableInspector.jl` | the ImGui table panel + its GLMakie plot | `DataBrowserGUI` |
| `src/Cache.jl`, `src/Cache/` | DuckDB persistence, cache buffer, build metrics | `DataBrowserCache` |
| `src/ItemIndex.jl`, `src/ItemIndex/` | `DataItem`, `ItemRecord`, `Hierarchy`, scanning, cancellation | `DataBrowserCore` |
| `src/Workspace.jl`, `src/Workspace/` | workspace, operations, processing, data access, status | `DataBrowserCore` (+ `WorkGraph` extracted) |
| `src/Visualization.jl` | plot-kind declarations **and** GLMakie rendering | split: declarations → API, rendering → `DataBrowserPlots` |
| `src/Browser.jl`, `src/Browser/`, `src/Gui/` | the CImGui shell, panels, browser | `DataBrowserGUI` |
| `src/MeasurementBrowser.jl` | top-level module, includes, exports | residual → `DataBrowser` umbrella |

## Dependency leaks the split has to sever

These are the couplings that keep the current package one indivisible unit; the extraction order is
chosen to cut them:

- **`DataFrame` in the "API" region.** `Project.jl` does `using DataFrames`, and `ItemIndex.jl` types
  signatures on `AbstractDataFrame`. `DataBrowserAPI` must be payload-agnostic, so none of that goes
  into API: the abstract `AbstractDataItem` contract (already in `Projects.jl`, dependency-free) is
  API; the concrete `DataItem` / `ItemRecord` and any `DataFrame`-typed signature stay in Core.
- **`GLMakie: Figure` in `Visualization.jl`.** Split it: the plot-kind *declarations* (`PlotKind`,
  `plot_kinds`, `setup_plot` / `plot_data!` as bare generics) go to API; the GLMakie rendering goes to
  `DataBrowserPlots`. After this, Core no longer transitively needs GLMakie.
- **`MB_BENCH_ENGINE_ONLY`.** The `ENGINE_ONLY_BENCHMARK_LOAD` branch in `MeasurementBrowser.jl`
  exists only to fake a headless engine inside the monolith. Once `DataBrowserCore` is a real headless
  package, `bench/` depends on Core directly and this branch is deleted.

## Target layout

```
Project.toml                 # residual → DataBrowser umbrella (→ GUI)
src/DataBrowser.jl
lib/
  DataBrowserProfiling/
  DataBrowserAPI/
  DataBrowserAnnotations/
  DataBrowserSources/
  DataBrowserCache/
  DataBrowserCore/
  DataBrowserPlots/
  DataBrowserGUI/
```

The umbrella's `Project.toml` lists the `lib/*` packages under `[sources]` as path deps; each `lib/*`
package lists its own `[sources]` entries for the DataBrowser packages it depends on. Only packages we
actually publish become registry entries; the rest resolve purely through paths.

## The steps

Each step: create `lib/DataBrowser<Name>` with a `Project.toml`, move the code, add `[sources]` wiring,
run `julia --project --threads=4 -e 'using Pkg; Pkg.test()'`, commit. The step is done only when the
suite is green.

1. **`DataBrowserProfiling`** — move `Profiling.jl`. Deps: `JSON`, `Profile`, `Statistics`. No
   DataBrowser deps. Smallest extraction; shakes out the `[sources]` mechanics before anything
   load-bearing.
2. **`DataBrowserAPI`** — move `Projects.jl` wholesale (the types + generic declarations), then move
   *only* the construction half of `Project.jl` (`define_project`, `register_item!`,
   `register_collection_analysis!`, `register_plot!` — the parts that just build a `Project` value)
   and the plot-kind declarations carved out of `Visualization.jl`. Leave every method that drives
   callbacks, and everything touching `DataFrame`, behind for Core. Deps: `Dates` only. *Boundary to
   enforce:* no `DataFrame`, no engine type, no heavy dep enters API.
3. **`DataBrowserAnnotations`** — move `src/Annotations/` to `lib/`, add the dep on `DataBrowserAPI`,
   and retype tags/notes/coords/layout off `String` paths onto the API's item-identity type. (Doing
   the retype here is what makes it a genuine `→ API` package rather than the current string-keyed
   leaf.)
4. **`DataBrowserSources`** — move `DataSources/DirectorySource.jl` and the raw-file-preview half of
   `TableInspector.jl` (`TablePreview` / `inspect_table` — reading and detecting a delimited file).
   Deps: `DataBrowserAPI`, `DataBrowserProfiling`, plus `CSV`, `DataFrames`, `BetterFileWatching`. The
   table *viewer* (`InspectorTable` / `merge_item_tables`) is not a source concern — it goes to Plots
   at step 7; the ImGui panel `src/Gui/TableInspector.jl` goes to GUI at step 8.
5. **`DataBrowserCache`** — move `Cache.jl` + `Cache/`. Deps: `DataBrowserAPI`,
   `DataBrowserProfiling`, plus `DuckDB`, `DBInterface`, `DataFrames`, `SHA`. *Boundary to enforce:*
   cache identity is the **project** fingerprint (environment + data identity + parameters), never the
   workspace.
6. **`DataBrowserCore`** — move `ItemIndex.jl` + `ItemIndex/`, `Workspace.jl` + `Workspace/`, and the
   engine-interface methods left behind from `Project.jl`, and lift the `WorkGraph` scheduler out of
   `Workspace/Operations.jl` + `Processing.jl` into its own module. Deps: API, Sources, Cache,
   Annotations, Profiling — **not** GLMakie. *Boundary to enforce:* callbacks receive the workspace as
   run-time context, not as an owning instance.
7. **`DataBrowserPlots`** — move the GLMakie rendering half of `Visualization.jl`, plus the table
   viewer's data model (`InspectorTable` / `merge_item_tables`). Deps: `DataBrowserAPI`, `GLMakie`.
   The viewer currently hardwires `DataFrame` and GLMakie; generalizing it to arbitrary payloads is
   Stage 2 generic-plotter work, not this mechanical move.
8. **`DataBrowserGUI`** — move `Browser.jl` + `Browser/` + `Gui/`. Deps: API, Core, Plots,
   Annotations, plus `CImGui`, `GLFW`, `ModernGL`, `NativeFileDialog`, `Observables`.
9. **Umbrella** — what remains of `MeasurementBrowser.jl` is the include-and-export shell. Rename the
   package and module to `DataBrowser`, depend on `DataBrowserGUI`, re-export the public API, wire
   defaults, delete `Precompile.jl`'s monolith assumptions and the `MB_BENCH_ENGINE_ONLY` branch.

## Keeping it green between steps

Until a package is extracted, its code stays where it is and the main package still includes it, so
there is never a half-moved module. The only churn per step is: the extracted files, their new
`Project.toml`, and the `import ..Foo` → `using DataBrowserFoo` swaps at the call sites. Tests run
against the same fixtures throughout; a step that can't stay green means its boundary is wrong and gets
revised before moving on.

## Open questions to settle during extraction

- **Exact API/Core line for the item model.** `AbstractDataItem` is clearly API; `DataItem` /
  `ItemRecord` are Core. Confirm nothing in a shared signature forces an engine type up into API.
- **`WorkGraph` module shape.** Whether it comes out as an internal Core module or its own `lib/`
  package depends on whether anything but Core needs it; default to a Core-internal module until a
  second consumer appears.
- **Annotations identity type.** Which API type annotations key on (item id vs. collection path)
  falls out of the data-model work in [data-model-generalization.md](data-model-generalization.md).
