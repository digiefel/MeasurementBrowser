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
item data. It combines GLMakie plotting with a CImGui UI. The package ships **no
domain knowledge** — no item kinds, readers, analyzers, or plotters. All of that lives in
*project code* (external scripts) registered through the project API. RuO2 ferroelectric/semiconductor
work and the TASE project are projects built on top, not part of the package.

### Modules (root `src/MeasurementBrowser.jl`, in include order)

| Module | File | Owns |
|---|---|---|
| `Projects` | `Projects.jl` | The low-level source contract (`AbstractDataSource`/`AbstractDataSourceItem`/`AbstractDataItem` + interface stubs), the `Project` recipe struct, and recipe types. Methods live in `Project.jl`. |
| `ItemIndex` | `ItemIndex.jl` (+ `Scanning.jl`) | The internal `ItemRecord`, the concrete `DataItem`, the hierarchy, and source-neutral scanning. |
| `DirectorySource` | `DataSources/DirectorySource.jl` | Directory traversal, `SourceFile`, file fingerprints, sidecar exclusion, timestamps, `metadata.txt`, and `open_workspace(project, root_path)`. |
| `Cache` | `Cache.jl`, `Cache/ProjectCache.jl` | Source-identity-keyed DuckDB cache of the scan + item data. |
| `Workspace` | `Workspace.jl` (+ `DataAccess.jl`, `Operations.jl`) | One open source: index, selection, loaded/cached data access, background jobs. |
| `Visualization` | `Visualization.jl` | Plot engine; `RegistryPlot{Kind}` adapter bridging registered plots to source plot dispatch. |
| `TableInspector` | `TableInspector.jl` | Tabular preview of a delimited file. |
| (recipes) | `Project.jl` | `define_project`, the `register_*` functions, project-side source interpretation/loading methods, and serialization. |
| `Browser` | `Browser.jl` (+ `Browser/*`, `Gui/*`) | CImGui frontend: tree/plot/table panels, state, persistence, tags, Makie embedding. |

### Two-layer API

The engine is written against a **low-level source contract** (`AbstractDataSource` →
`AbstractDataSourceItem` → `AbstractDataItem`; see [docs/api.md](docs/api.md)). The **high-level
callback API** below (`define_project` + `register_*`) is a convenience on top of that contract:
`DirectorySource` discovers files, and project methods interpret/load them. The callback API is the
exported surface; the low-level types are internal for now (`MeasurementBrowser.name`), staged for a
future submodule.

### Project API (high-level callback convenience)

Project code builds a `Project` by mutation, pointing each registration at plain callbacks:

```julia
project = define_project("TASE"; description="…")

register_item!(project, :iv;
    detect  = file -> occursin("iv", lowercase(file.filename)),  # Bool; first match wins
    read    = file -> DataFrame(...),                            # whole file, parsed once
    entries = (file, data) -> [DataItem(kind=:iv, collection=[...], parameters=..., data=data)],
    process = item -> DataItem(item, clean(item.data)),          # optional; default passthrough
    stats   = item -> Dict{Symbol,Any}(...),                     # optional; background analysis
    label   = item -> "…")                                       # optional

register_collection_stat!(project; kinds=[:iv],                  # cross-item fold over a collection
    compute_stats = group -> …)                                  # returns a Dict stored on the collection node

register_plot!(project, :iv; label="I–V", setup=…, draw=…)       # one plot recipe per kind
```

The `entries` callback enumerates the items in a file and is where a single file expands into multiple
entries (e.g. one item per fatigue cycle, distinguished by `parameters`). When an entry omits `id`,
the directory adapter mints one from the source-item path, kind, and params. Entries attach the raw
per-item payload as `item.data`; optional `process(item)` returns the item that stats and views
receive, and optional `stats(item)` runs after indexing in workspace background analysis.
`entries` returns the package's `DataItem` (the **recipe API**) — or, to go beyond it, a project's own
`AbstractDataItem` subtype carrying typed fields and its data (the **type API**); the engine derives
the internal record from either via the contract. Re-calling a `register_*` with the same key replaces
the recipe in place, so editing and re-running a line updates a live project.

### Data Flow

1. **Scan:** `scan_source!` calls `source_items(source)` for the discovered units. For each unit,
   `data_items(project, source, source_item)` interprets it into data items (the adapter applies the first
   matching `detect`'s `read → entries`). The engine derives each item's internal `ItemRecord` via the
   contract; per-source-item failures are recorded and scanning continues. The file-backed adapter
   provides collection parameters from `metadata.txt`; the index inherits them through the hierarchy
   and merges them into each record's effective `parameters`.
2. **Analysis:** after the tree exists, a workspace background job materializes items with
   `read → entries → process`, computes per-item `stats(item)`, then runs
   `collection_stats(project, source, collection, items)` (the callback form is
   `register_collection_stat!`) for collection nodes.
3. **Cache:** the `SourceScan` is stored in a source-id-keyed DuckDB cache; a fresh cache
   restores the hierarchy quickly while scanning continues.
4. **View:** for the selection, the engine reloads items via
   `load_data_item(project, source, source_item_id, id)` (memory → cache → origin); each carries `item.data`.
   Registered plots render through the internal plot dispatch, with
   `setup(ws, items)` / `draw(ws, items, fig)` reading `item.data`.

### Key Types (`ItemIndex.jl`, `Projects.jl`)

- `AbstractDataSource` / `AbstractDataSourceItem` / `AbstractDataItem` — the low-level contract: a
  source owns lifecycle + discovery (`source_items`), a source item is one discovered unit
  (`data_items`, `load_data_item`), a data item is one logical browser object (`id`, `item_label`,
  `kind`, `collection`, `parameters`, `stats`, `item_data`, + optional `process`/`cacheable`/
  `fingerprint`). The engine is written against these, never a concrete type.
- `DataItem <: AbstractDataItem` — the normal item the recipe API's `entries` produces: `kind`,
  `collection`, `parameters`/`stats` dicts, and `data` (as `item.data`).
- `SourceFile <: AbstractDataSourceItem` — the built-in file-backed source item: path, filename,
  timestamp, fingerprint. A bare `SourceFile` has no universal `data_items`.
- `ItemRecord` — **internal**, data-less record the hierarchy/scan/cache store: source-item identity
  (`source_item_id`/label/fingerprint/path/timestamp), item identity (`id`, `item_label`,
  `kind`), `collection`, `parameters`, `stats`, `item_fingerprint`. Never named
  by source/project code, never a field of any item; the engine converts item ↔ record via the contract.
- `SourceScan` — internal scan result (source-neutral; source-level identity lives once on the scan).
  `Hierarchy` / `HierarchyNode` — the browsable collection tree.
- `Project` / `ItemRecipe` / `CollectionStatRecipe` / `PlotRecipe` — the callback registry model;
  `DirectorySource` supplies the built-in directory-backed source.

## Launching

```julia
using MeasurementBrowser
# … define_project + register_* …
ws = open_workspace("/path/to/data"; project)   # or open_workspace(mysource) for a low-level source
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
- Scanning large source trees is slow; point at a small subfolder or `test/` fixtures during
  development.
- **Detection order matters:** recipes are checked in registration order and the first `detect`
  returning `true` wins — register more specific patterns (e.g. `pund_fatigue`) before general ones
  (`pund`).
- Docs split: `docs/*.md` are current-state reference; `docs/plans/*.md` are forward-looking design.
