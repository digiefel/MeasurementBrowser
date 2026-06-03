# MeasurementBrowser Architecture

> **Onboarding doc.** Start here before touching code. Skim → click into the focused docs below for the area you're working on. Keep this file in sync when structure changes.

## What this is

A Julia GUI app for browsing and analyzing ferroelectric and semiconductor device measurements (RuO2 test chips). CImGui drives the UI; GLMakie renders plots embedded in the UI via a custom integration layer. CSV files on disk are the data; source-root text files are metadata; the HDF5 cache is generated acceleration for project loading and plotting.

## Three-stage pipeline

```
                ┌─────────────────────────────────────────────────────┐
disk (CSVs +    │  scan      →   hierarchy   →   plot                 │
metadata)       │  --------     -----------     ------                │
                │  walk fs      build tree      figure_for_file(s)    │
                │  parse        per-leaf        Makie figure embedded │
                │  filenames    measurements    in CImGui via         │
                │  merge meta                   MakieIntegration      │
                └─────────────────────────────────────────────────────┘
```

1. **Scan** ([scanning](scanning.md)): `scan_source` walks a folder of CSVs, parses each filename to a `MeasurementInfo`, merges in `device_info.txt` metadata, and inserts results into a tree.
2. **Hierarchy** ([data-model](data-model.md)): a variable-depth `MeasurementHierarchy` whose leaves hold measurements. Stable identity is the slash-joined `device_path_key`.
3. **Cache** ([cache](cache.md)): generated HDF5 data can repopulate the browser tree and serve normal plot jobs without reparsing every CSV.
4. **Plot** ([gui](gui.md), [plotting](plotting.md)): the GUI's tree panel drives a selection vector; selection feeds cached or source-loaded plot data into Makie figures. Figures render asynchronously via a `PlotJob` queue.

## Module map

```
MeasurementBrowser.jl (root)
  ├── DeviceParser.jl       — filename parsing, hierarchy types, device_info.txt
  ├── ProjectCache.jl       — HDF5 cache identity, load, status, build, update
  ├── PlotJobs.jl           — async figure rendering job queue
  ├── Gui.jl                — slim shell; includes files in src/Gui/        [see gui.md]
  │     └── Gui/
  │         ├── State.jl         — ui_state init + scan/cache/figure-script state machines
  │         ├── BadAndStyling.jl — bad-registry helpers, visibility, styling
  │         ├── PlotPanel.jl     — plot job queue, plot/combined plot windows
  │         ├── TreePanel.jl     — hierarchy tree + measurements panel + selection
  │         ├── InfoModal.jl     — info window, figure-script window, device-info modal
  │         └── Layout.jl        — frame loop, docking, menu bar, perf/project windows
  ├── MakieIntegration.jl   — GLMakie ↔ CImGui screen-per-title bridge
  └── projects/RuO2/        — project-specific filename interpreters & cache payload schema

Subpackages (path-deps; each has its own Project.toml):
  src/DataLoader/    — CSV readers (read_fe_pund, read_iv_sweep, …)
  src/DataAnalysis/  — analyzers (analyze_pund, analyze_tlm_combined, …)
  src/DataPlotter/   — figure_for_file, figure_for_files; depends on Loader+Analysis
  src/Annotations/   — on-disk annotation formats (Coords, Layout, Tags, Notes)
```

When you change a subpackage, run `Pkg.instantiate()` in its directory too — see [CLAUDE.md](../CLAUDE.md).

## On-disk metadata

All metadata lives at the **source root** alongside the CSVs (not inside the repo). Lean, hand-editable text. See [storage](storage.md) for formats.

- `device_info.txt` — per-device parameters with path-fragment matching (area, thickness, etc.).
- `tags.txt` — tag catalog and per-path assignments.

The HDF5 cache is generated and lives outside the source root. See [cache](cache.md).

## Conventions cheatsheet

- Julia 1.12; 4-space indentation; ~100-char lines.
- `snake_case` functions/variables, `UpperCamelCase` types, `ALL_CAPS` constants.
- Don't catch errors; let failures surface.
- Update `detect_measurement_kind` ordering carefully — specific patterns before general ones.
- Tests in `test/` with fixture CSVs under `test/fixtures/`.
- Figure-creation tests assert metadata/labels, not pixels.

## Deep-dive docs

| Topic | When to read |
|---|---|
| [data-model.md](data-model.md) | You're touching `DeviceParser`, hierarchy traversal, or paths/IDs. |
| [gui.md](gui.md) | You're adding/modifying a panel, window, or interaction. |
| [storage.md](storage.md) | You're adding a new metadata file or changing a format. |
| [cache.md](cache.md) | You're touching HDF5 cache identity, loading, writing, status, or plot payloads. |
| [annotations.md](annotations.md) | You're touching `src/Annotations/` — coords, layout, tags, notes. |
| [figure_scripts.md](figure_scripts.md) | You're working on the figure-script export feature. |

In-flight design and roadmaps live under [plans/](plans/), separate from the current-state docs above.

## Doc maintenance rule

**Any change that affects structure described in these docs must update them in the same commit.** That includes:

- Adding/removing/renaming files in `src/`.
- Adding new on-disk metadata formats.
- Changing public APIs of subpackages.
- Restructuring `Gui.jl` panels or `ui_state` keys.

If you're not sure whether a change is "structural enough" to document, it probably is. Err on the side of updating.

These docs describe **current state only**. No history, no "used to be", no references to in-flight work or open plans. If a change supersedes an old structure, replace the description rather than annotating the transition. Past states belong in `git log`; in-flight feature plans belong in plan files.
