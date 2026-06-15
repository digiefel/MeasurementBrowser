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
| `Projects` | `Projects.jl` | The `Project` struct + recipe types + interface stubs. Defined early so other modules can name the concrete type; its methods live in `Project.jl`. |
| `MeasurementIndex` | `MeasurementIndex.jl` (+ `SourceFiles.jl`, `Scanning.jl`) | `SourceFile`, `MeasurementInfo`, `DeviceInfo`, the hierarchy, and filesystem scanning. |
| `Cache` | `Cache.jl`, `Cache/ProjectCache.jl` | Fingerprint-keyed HDF5 cache of the scan + measurement data. |
| `Workspace` | `Workspace.jl` (+ `DataAccess.jl`, `Operations.jl`) | One open data root: index, selection, loaded/cached data access, background jobs. |
| `Visualization` | `Visualization.jl` | Plot engine; `RegistryPlot{Kind}` adapter bridging registered plots to dispatch. |
| `TableInspector` | `TableInspector.jl` | Tabular preview of measurement data. |
| (project methods) | `Project.jl` | `define_project`, the `register_*` functions, the engine-interface impls, and `Project` serialization. |
| `Browser` | `Browser.jl` (+ `Browser/*`, `Gui/*`) | CImGui frontend: tree/plot/table panels, state, persistence, tags, Makie embedding. |

### Project API

Project code builds a `Project` by mutation, pointing each registration at plain callbacks:

```julia
project = define_project("TASE"; description="…")

register_measurement!(project, :iv;
    detect       = file -> occursin("iv", lowercase(file.filename)),  # Bool; first match wins
    read         = file -> DataFrame(...),                            # whole file, parsed once
    measurements = (file, data) -> [MeasurementInfo(...)],            # one entry per measurement
    process      = (mi, data) -> data,                               # optional; default passthrough
    stats        = (mi, processed) -> Dict{Symbol,Any}(...),         # optional
    label        = mi -> "…")                                        # optional

register_device_stat!(project; measurement_kinds=[:iv],              # cross-measurement fold
    compute_stats = group -> …)                                      # fills each mi.stats in place

register_plot!(project, :iv; label="I–V", setup=…, draw=…)           # one plot recipe per kind
```

The `measurements` callback sets each measurement's device identity (`DeviceInfo`) and is where a
single file expands into multiple entries (e.g. one `MeasurementInfo` per fatigue cycle, distinguished
by `parameters`); the engine mints each `unique_id` from filepath + kind + params. Re-calling a
`register_*` with the same key replaces the recipe in place, so editing and re-running a line updates
a live project.

### Data Flow

1. **Scan:** `scan_source` walks the data root, and for each file the project's first matching
   `detect` selects a recipe. `interpret_file` then runs `read → measurements`, and folds each
   measurement's `process → stats` into the **same pass** (so the parsed table is freed per file and
   scan memory stays bounded). Merges source-root `device_info.txt` metadata.
2. **Hierarchy:** measurements organize into a `MeasurementHierarchy` tree by `DeviceInfo.location`.
   Cross-measurement `register_device_stat!` folds run after the full scan.
3. **Cache:** the whole scan (including the `Project`'s recipes) is serialized into a fingerprint-keyed
   HDF5 cache; a fresh cache restores the hierarchy quickly while scanning continues.
4. **View:** views request data via `read_measurement_data`/`process_measurement_data` (memory →
   cache → source). Registered plots render through the `RegistryPlot{Kind}` adapter.

### Key Types (`MeasurementIndex.jl`, `Projects.jl`)

- `MeasurementInfo` — one logical measurement: ids, `clean_title`, `measurement_kind::Symbol`,
  `timestamp`, `device_info`, `parameters` (known at interpret time), `stats` (computed later).
- `DeviceInfo` — `location::Vector{String}` (the tree path) + `parameters` metadata dict.
- `SourceFile` — one physical file: path, timestamp, header summary, fingerprint, measurements.
- `MeasurementHierarchy` / `HierarchyNode` — the browsable device tree.
- `Project` / `MeasurementRecipe` / `DeviceStatRecipe` / `PlotRecipe` — the registry model.

## Launching

```julia
using MeasurementBrowser
# … define_project + register_* …
ws = open_workspace(project, "/path/to/data")
open_browser(ws)
```

Project scripts live outside the repo (e.g. the TASE_SNS and v2 RuO2 analysis projects). They add
`MeasurementBrowser` as a dependency, register their measurement kinds, and launch the browser.

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
