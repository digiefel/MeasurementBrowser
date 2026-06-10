# MeasurementBrowser Architecture

> **Onboarding doc.** Start here before touching code. This file describes the architectural model
> that is not obvious from the filesystem. It should not be a file inventory.

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

## Ownership Boundaries

The most important boundary is between project meaning and package machinery.

Project code owns:

- mapping source files to logical measurements
- loading source data into reusable tables
- computing project-specific stats or processed data
- defining project-specific visualizers when generic ones are not enough

Package code owns:

- opened-project state
- cache identity, freshness, storage, and repair
- source scanning and background jobs
- hierarchy, selection, and GUI state
- generic visualizers, plot windows, workflow persistence, and figure composition

Project code should not know whether data came from memory, cache, or source files. Package code
should not know the experimental meaning of a project file beyond the project-facing functions it
calls.

The package expresses that boundary through focused internal modules. `Project` defines the methods
implemented by a measurement project. `MeasurementIndex` owns source-file records, logical
measurements, hierarchy construction, and filesystem scanning. `Cache` owns generated HDF5 state.
`Workspace` owns cache-aware measurement data access. `Visualization` defines the shared plotting
operations used by projects and the browser. `MeasurementBrowser` exports the small project-facing
API while keeping these modules internal.

## On-disk metadata

All metadata lives at the **source root** alongside the CSVs (not inside the repo). Lean, hand-editable text. See [storage](storage.md) for formats.

- `device_info.txt` — per-device parameters with path-fragment matching (area, thickness, etc.).
- `tags.txt` — tag catalog and per-path assignments.

The HDF5 cache is generated and lives outside the source root. See [cache](cache.md).

## Engineering Rules

This project is pre-pre-alpha. Prioritize simple, principled code and development velocity over
backwards compatibility. Do not keep compatibility paths, fallback behavior, or legacy symbols
unless there is a clear current need. When an old shape is wrong, replace it cleanly and update the
docs and tests to match.

- Julia 1.12; 4-space indentation; ~100-char lines.
- `snake_case` functions/variables, `UpperCamelCase` types, `ALL_CAPS` constants.
- Don't catch errors that should be fixed; let failures surface clearly.
- Update `detect_measurement_kind` ordering carefully — specific patterns before general ones.
- Tests live in `test/` with fixture CSVs under `test/fixtures/`.
- Figure-creation tests assert metadata and figure construction behavior, not pixels or labels.

## Deep-dive docs

Current-state reference docs:

| Topic | When to read |
|---|---|
| [data-model.md](data-model.md) | You're touching source-file records, measurements, hierarchy traversal, or paths/IDs. |
| [gui.md](gui.md) | You're adding/modifying a panel, window, or interaction. |
| [storage.md](storage.md) | You're adding a new metadata file or changing a format. |
| [cache.md](cache.md) | You're touching HDF5 cache identity, loading, writing, status, or measurement data. |
| [annotations.md](annotations.md) | You're touching `src/Annotations/` — coords, layout, tags, notes. |
| [figure_scripts.md](figure_scripts.md) | You're working on the figure-script export feature. |

Future-facing docs and roadmaps live under [plans/](plans/), separate from the current-state docs
above. Start with [plans/README.md](plans/README.md) for the plan map and
[plans/workspace-vision.md](plans/workspace-vision.md) for the long-term GUI/API/workflow vision.

## Doc maintenance rule

Update current-state docs when a change affects the model a reader cannot infer by grepping files:

- data identity or hierarchy semantics
- cache/source correctness rules
- source-root metadata formats
- public package or subpackage API contracts
- GUI state ownership or cross-panel interaction rules

Do not update architecture docs just because files moved, were renamed, or were split. The
filesystem already records that. If a doc becomes a stale migration note or repeats another doc,
delete it or fold the useful idea into the owning doc.
