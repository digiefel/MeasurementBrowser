# MeasurementBrowser Architecture

> **Onboarding doc.** Start here before touching code. Skim → click into the focused docs below for the area you're working on. Keep this file in sync when structure changes.

## What this is

MeasurementBrowser is a Julia working environment for measurement projects. It should be strictly
better than opening a REPL, finding files, loading them, extracting useful tables, computing values,
and writing figure code by hand. The app makes that workflow interactive: open a project, browse the
measurement structure, select measurements or devices, inspect data-derived values, and switch plots
or views without repeating the same parsing work.

CImGui drives the UI; GLMakie renders plots embedded in the UI via a custom integration layer. CSV
files on disk are the source data; source-root text files are metadata; the generated cache lets the
package reuse project structure and measurement data without making cache storage part of project
code.

Project code should describe the project in the same terms a script would use: which measurements a
source file contains, how to load the data for those measurements, and how to present that data. The
package owns scanning, cache storage, background jobs, and UI state.

The intended experience is live and composable. A project can stay open while source files are added
or changed, and the browser should update without forcing the user back through startup or manual
reload steps. Views should be able to follow selections, devices, or matching rules, so a plot or
inspection tool can continue to show the relevant data as the source tree changes. Built-in
visualizers should handle common inspection tasks, while project code adds only the interpretation
and presentation details that are specific to the experiment.

Keep this boundary small. It should be possible to evolve project code toward hot reloading,
interpreted execution, or an out-of-process implementation without changing what a project author is
trying to express: source files become measurements, measurements provide reusable data, and reusable
data feeds inspection and presentation tools.

## Core Flow

```
                ┌─────────────────────────────────────────────────────┐
disk (CSVs +    │  scan      →   hierarchy   →   data      →   view   │
metadata)       │  --------     -----------     ----          ----     │
                │  walk fs      build tree      cache or      Makie  │
                │  interpret    per-leaf        source        figure │
                │  source       measurements    load          table  │
                │  files        + stats                      summary │
                └─────────────────────────────────────────────────────┘
```

`scan_source` walks source CSVs, asks the active project which measurements each file contains, merges
source-root metadata, computes project stats, and builds the hierarchy described in
[data-model.md](data-model.md). The cache described in [cache.md](cache.md) can restore that
hierarchy quickly while source scanning continues in the background.

When a view needs measurement data, it goes through `read_measurement_data(project, measurements)`.
That package-owned path returns cached dataframe data when it is valid and falls back to the
project's `load_source_data` when the cache is missing or stale. GUI selection, background jobs, and
Makie embedding are described in [gui.md](gui.md).

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
  └── projects/RuO2/        — project-specific source interpretation, data loading, and plotting

Subpackages (path-deps; each has its own Project.toml):
  src/DataAnalysis/  — analyzers (analyze_pund, analyze_tlm_combined, …)
  src/Annotations/   — on-disk annotation formats (Coords, Layout, Tags, Notes)
```

`src/DataLoader/` still exists for old generic CSV readers, but it is deprecated. New source readers
belong in project code.

When you change an active subpackage, run `Pkg.instantiate()` in its directory too.

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
| [cache.md](cache.md) | You're touching HDF5 cache identity, loading, writing, status, or measurement data. |
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
