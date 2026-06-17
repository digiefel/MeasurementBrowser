# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

For the architectural model, read [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md) first — it is the
current-state source of truth. This file is the operational quick-reference.

## Commands

```bash
# Install dependencies
julia --project -e 'using Pkg; Pkg.instantiate()'

# Run tests
julia --project -e 'using Pkg; Pkg.test()'

# Precompile without launching UI
julia --project -e 'using Pkg; Pkg.precompile()'
```

`Annotations` (`src/Annotations`) is the only local subpackage; instantiate it too when its deps
change: `julia --project=src/Annotations -e 'using Pkg; Pkg.instantiate()'`.

There is no generic launcher script. The app is started from **project code**: define a project,
open a workspace on a data root, and open the browser. See "Launching" below.

## Architecture

MeasurementBrowser is a **project-agnostic** Julia GUI engine for browsing and analyzing
measurement data. It combines GLMakie plotting with a CImGui UI. The package ships **no
domain knowledge** — no measurement kinds, readers, analyzers, or plotters. All of that lives in
*project code* (external scripts) registered through the project API. RuO2 ferroelectric/semiconductor
work and the TASE project are projects built on top, not part of the package.

### Modules (root `src/MeasurementBrowser.jl`, in include order)

| Module | File | Owns |
|---|---|---|
| `Projects` | `Projects.jl` | The `Project` struct + recipe types + the `AbstractDataItem` contract + interface stubs. Defined early so other modules can name the concrete type; its methods live in `Project.jl`. |
| `ItemIndex` | `ItemIndex.jl` (+ `SourceFiles.jl`, `Scanning.jl`) | `SourceFile`, the internal `ItemRecord`, the concrete `DataItem`, the hierarchy, and filesystem scanning. |
| `Cache` | `Cache.jl`, `Cache/ProjectCache.jl` | Fingerprint-keyed HDF5 cache of the scan + item data. |
| `Workspace` | `Workspace.jl` (+ `DataAccess.jl`, `Operations.jl`) | One open data root: index, selection, loaded/cached data access, background jobs. |
| `Visualization` | `Visualization.jl` | Plot engine; `RegistryPlot{Kind}` adapter bridging registered plots to dispatch. |
| `TableInspector` | `TableInspector.jl` | Tabular preview of a delimited file. |
| (project methods) | `Project.jl` | `define_project`, the `register_*` functions, the engine-interface impls, and `Project` serialization. |
| `Browser` | `Browser.jl` (+ `Browser/*`, `Gui/*`) | CImGui frontend: tree/plot/table panels, state, persistence, tags, Makie embedding. |

### Project API

Project code builds a `Project` by mutation, pointing each registration at plain callbacks:

```julia
project = define_project("TASE"; description="…")

register_item!(project, :iv;
    detect  = file -> occursin("iv", lowercase(file.filename)),  # Bool; first match wins
    read    = file -> DataFrame(...),                            # whole file, parsed once
    entries = (file, data) -> [DataItem(kind=:iv, collection=[...], parameters=...)],  # one per item
    process = (item, data) -> data,                              # optional; default passthrough
    stats   = (item, processed) -> Dict{Symbol,Any}(...),        # optional
    label   = item -> "…")                                       # optional

register_collection_stat!(project; kinds=[:iv],                  # cross-item fold over a collection
    compute_stats = group -> …)                                  # fills each record's `stats` in place

register_plot!(project, :iv; label="I–V", setup=…, draw=…)       # one plot recipe per kind
```

The `entries` callback enumerates the items in a file and is where a single file expands into multiple
entries (e.g. one item per fatigue cycle, distinguished by `parameters`); the engine fills each item's
`filepath`/`timestamp` from the `SourceFile` and mints `unique_id` from filepath + kind + params.
`entries` returns the package's `DataItem` (the **recipe API**) — or, to go beyond it, a project's own
`AbstractDataItem` subtype carrying typed fields and its data (the **type API**); the engine derives
the internal record from either via the contract. Re-calling a `register_*` with the same key replaces
the recipe in place, so editing and re-running a line updates a live project.

### Data Flow

1. **Scan:** `scan_source` walks the data root (every visible file except package sidecars is a
   candidate), and for each file the project's first matching `detect` selects a recipe; files no
   recipe claims are skipped. `interpret_file` runs `read → entries`, derives each item's internal
   `ItemRecord` via the contract, and folds `process → stats` into the **same pass** (so the parse is
   freed per file and scan memory stays bounded). Merges source-root `device_info.txt` metadata into
   each record's `collection_metadata`.
2. **Hierarchy:** records organize into a `Hierarchy` tree by their `collection::Vector{String}`.
   Cross-item `register_collection_stat!` folds run after the full scan.
3. **Cache:** the whole scan (including the `Project`'s recipes) is serialized into a fingerprint-keyed
   HDF5 cache; a fresh cache restores the hierarchy quickly while scanning continues.
4. **View:** the engine materializes data-bearing items for the selection (`materialize_items`): the
   recipe API resolves data through `read_item_data`/`process_item_data` (memory → cache → source) and
   wraps it in a `DataItem`; the type API re-runs `read`/`entries` and returns the project's own items.
   Either way each item carries `item.data`. Registered plots render through the `RegistryPlot{Kind}`
   adapter, with `setup(ws, items)` / `draw(ws, items, fig)` reading `item.data`.

### Key Types (`ItemIndex.jl`, `Projects.jl`)

- `AbstractDataItem` — the contract every item answers: `item_id`, `item_label`, `kind`, `collection`,
  `parameters`, `stats`, `item_data`, `read_data`, plus optional `process`/`cacheable`. The engine is
  written against this, never a concrete item type.
- `DataItem <: AbstractDataItem` — the normal item the package ships and `register_item!`'s `entries`
  produces: `kind`, `collection`, `parameters`/`stats` dicts, and `data` (as `item.data`).
- `ItemRecord` — **internal**, data-less metadata record the hierarchy/scan/cache store: ids,
  `clean_title`, `kind`, `timestamp`, `collection`, `collection_metadata`, `parameters`, `stats`. Never
  named by project code, never a field of any item; the engine converts item ↔ record via the contract.
- `SourceFile` — one physical file: path, timestamp, fingerprint, and the records interpreted from it.
- `Hierarchy` / `HierarchyNode` — the browsable collection tree.
- `Project` / `ItemRecipe` / `CollectionStatRecipe` / `PlotRecipe` — the registry model.

## Launching

```julia
using MeasurementBrowser
# … define_project + register_* …
ws = open_workspace(project, "/path/to/data")
open_browser(ws)
```

Project scripts live outside the repo (e.g. the TASE_SNS and v2 RuO2 analysis projects). They add
`MeasurementBrowser` as a dependency, register their item kinds, and launch the browser.

## Coding Conventions

- Julia 1.12; 4-space indentation; ~100-char lines
- `snake_case` functions/variables, `UpperCamelCase` types, `ALL_CAPS` constants
- Don't catch errors that should be fixed; let failures surface
- Docstrings for public APIs; brief comments only where non-obvious
- Short imperative commit messages matching existing history style
- Pre-pre-alpha: prefer clean replacement over compatibility shims (see ARCHITECTURE.md engineering rules)

## Testing

Tests in `test/` with fixture CSVs under `test/fixtures/`. Projects under test are built inline with
`define_project` + `register_*` (see `test/test_project.jl`, `test/test_scan_profile.jl`). For
plot/GUI changes, assert metadata/labels and that figure creation doesn't error rather than
pixel-perfect checks.

## Notes

- UI requires an OpenGL-capable environment (GLMakie + CImGui).
- Scanning large measurement trees is slow; point at a small subfolder or `test/` fixtures during
  development.
- **Detection order matters:** recipes are checked in registration order and the first `detect`
  returning `true` wins — register more specific patterns (e.g. `pund_fatigue`) before general ones
  (`pund`).
- Docs split: `docs/*.md` are current-state reference; `docs/plans/*.md` are forward-looking design.
