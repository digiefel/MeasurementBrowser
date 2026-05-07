# MeasurementBrowser Architecture

> **Onboarding doc.** Start here before touching code. Skim → click into the focused docs below for the area you're working on. Keep this file in sync when structure changes.

## What this is

A Julia GUI app for browsing and analyzing ferroelectric and semiconductor device measurements (RuO2 test chips). CImGui drives the UI; GLMakie renders plots embedded in the UI via a custom integration layer. CSV files on disk are the data; everything else is metadata or computed.

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

1. **Scan** ([scanning](scanning.md)): `scan_directory` walks a folder of CSVs, parses each filename to a `MeasurementInfo`, merges in `devices_info.txt` metadata (path-prefix matching), and inserts results into a tree.
2. **Hierarchy** ([data-model](data-model.md)): a variable-depth `MeasurementHierarchy` whose leaves hold measurements. Stable identity is the slash-joined `device_path_key`.
3. **Plot** ([gui](gui.md), [plotting](plotting.md)): the GUI's tree panel drives a selection vector; selection feeds `figure_for_file` (single) or `figure_for_files` (combined). Figures render asynchronously via a `PlotJob` queue.

## Module map

```
MeasurementBrowser.jl (root)
  ├── DeviceParser.jl       — filename parsing, hierarchy types, devices_info.txt
  ├── BadRegistry.jl        — load/save bad_measurements
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
  └── projects/RuO2/        — project-specific filename interpreters & cache schema

Subpackages (path-deps; each has its own Project.toml):
  src/DataLoader/    — CSV readers (read_fe_pund, read_iv_sweep, …)
  src/DataAnalysis/  — analyzers (analyze_pund, analyze_tlm_combined, …)
  src/DataPlotter/   — figure_for_file, figure_for_files; depends on Loader+Analysis
  src/Annotations/   — on-disk annotation formats (Coords, Layout, Tags, Notes)
```

When you change a subpackage, run `Pkg.instantiate()` in its directory too — see [CLAUDE.md](../CLAUDE.md).

## On-disk metadata

All metadata lives at the **source root** alongside the CSVs (not inside the repo). Lean, hand-editable text. See [storage](storage.md) for formats.

- `devices_info.txt` — per-device parameters with hierarchical path matching (area, thickness, etc.).
- `bad_measurements` — flat list of paths/IDs flagged bad.

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
| [annotations.md](annotations.md) | You're touching `src/Annotations/` — coords, layout, tags, notes. |
| [figure_scripts.md](figure_scripts.md) | You're working on the figure-script export feature. |

## Doc maintenance rule

**Any change that affects structure described in these docs must update them in the same commit.** That includes:

- Adding/removing/renaming files in `src/`.
- Adding new on-disk metadata formats.
- Changing public APIs of subpackages.
- Restructuring `Gui.jl` panels or `ui_state` keys.

If you're not sure whether a change is "structural enough" to document, it probably is. Err on the side of updating.

These docs describe **current state only**. No history, no "used to be", no references to in-flight work or open plans. If a change supersedes an old structure, replace the description rather than annotating the transition. Past states belong in `git log`; in-flight feature plans belong in plan files.
