# MeasurementBrowser Architecture

> **Onboarding doc.** Start here before touching code. This file describes the architectural model
> that is not obvious from the filesystem. It should not be a file inventory.

## What this is

MeasurementBrowser is a Julia working environment for data projects. It should be strictly
better than opening a REPL, finding files, loading them, extracting useful tables, computing values,
and writing figure code by hand. The app makes that workflow interactive: open a project, browse the
source structure, select collections or items, inspect data-derived values, and switch plots or views
without repeating the same parsing work.

CImGui drives the UI; GLMakie renders plots embedded in the UI via a custom integration layer. CSV
files on disk are one possible source; source-root text files can provide metadata; the generated
cache lets the package reuse project structure and item data without making cache storage part of
project code.

Project code should describe the project in the same terms a script would use: which logical items a
source item contains, how to load data for those items, and how to present that data. The
package owns scanning, cache storage, background jobs, and UI state.

The intended experience is live and composable. A project can stay open while source files are added
or changed, and the browser should update without forcing the user back through startup or manual
reload steps. Views should be able to follow selections, collections, or matching rules, so a plot or
inspection tool can continue to show the relevant data as the source tree changes. Built-in
visualizers should handle common inspection tasks, while project code adds only the interpretation
and presentation details that are specific to the experiment.

Keep this boundary small. It should be possible to evolve project code toward hot reloading,
interpreted execution, or an out-of-process implementation without changing what a project author is
trying to express: source items become logical items, loaded item data feeds inspection and
presentation tools.

## Core Flow

```
                 ┌──────────────────────────────────────────────────────────┐
AbstractDataSource│ scan          →  hierarchy  →  data       →  view        │
(root, DB, stream)│ ----------       ----------    ----           ----        │
                 │ source_items     build tree    load_data_item  Makie      │
                 │ + data_items     per-leaf       (memory →       figure /   │
                 │ per source item  items+stats    cache → source) table      │
                 └──────────────────────────────────────────────────────────┘
```

The engine is written against the **low-level source contract** ([api.md](api.md)): an
`AbstractDataSource` owns lifecycle and discovery, an `AbstractDataSourceItem` is one discovered unit,
and an `AbstractDataItem` is one logical browsable item. `scan_source!` calls `source_items(source)`,
interprets each unit with `data_items(project, source, source_item)`, records per-source-item failures,
computes collection-node stats via `collection_stats(project, source, collection, items)`, and builds the hierarchy described in
[data-model.md](data-model.md). The cache ([cache.md](cache.md)) restores that hierarchy quickly while
scanning continues in the background.

When a view needs item data, the engine reloads the selected items via
`load_data_item(project, source, source_item_id, item_id)` (memory → cache → source). Each item carries
`item.data`. GUI selection, background work, and Makie embedding are described in [gui.md](gui.md).

The **high-level callback API** (`define_project` + `register_*`, the exported convenience surface)
uses the built-in `DirectorySource <: AbstractDataSource`: the source walks a data root into
`SourceFile`s, while project methods apply the recipes' `detect`/`read`/`entries`/`stats` and
`read`/`process` loading. Nothing downstream of the contract can tell a callback project from a
hand-written project/source pair.

## Ownership Boundaries

The most important boundary is between source meaning and package machinery.

A source owns:

- discovering source items
- source-specific metadata, such as `DirectorySource` loading `metadata.txt`

A project/source implementation owns:

- interpreting each source item into logical data items
- loading one item's data on demand (`load_data_item`)
- computing per-item and per-collection stats
- defining project-specific visualizers when generic ones are not enough

The workspace owns:

- the open source and its identity
- the progressively populated item index
- selection identities
- cache identity, freshness, storage, and repair
- scanning, cache work, progress, errors, and cancellation
- loaded items already materialized in memory

The browser owns windows, controls, filters, and temporary rendering state. Annotations store
user-authored tags, notes, coordinates, and spatial positions (currently next to the cache, keyed by
`source_id`; see [cache.md](cache.md)). Other package modules own generic visualizers, workflow
persistence, and figure composition.

Source code should not know whether data came from memory, cache, or the origin. Package code should
not know the meaning of a source item beyond the contract methods it calls.

The package expresses that boundary through focused internal modules. `Projects` defines the
`AbstractDataItem` contract and the `Project` recipe type. `DataSources/DirectorySource.jl` owns the
built-in directory source, `SourceFile`, file fingerprints, directory traversal, sidecar exclusion,
and `metadata.txt` collection metadata. `ItemIndex` owns the internal `ItemRecord`, the concrete
`DataItem`, hierarchy construction, and scanning. `Cache` owns generated HDF5 state. `Workspace` owns
one project/source pair, its index, selection, cache, loaded data, and background work.
`Visualization` defines the shared plotting operations. `Browser` owns typed frontend state and
CImGui rendering. `MeasurementBrowser` exports the small high-level API while keeping most low-level
source contract names internal.

`open_workspace(project, source)` or `open_workspace(project, root_path)` creates that owner and
immediately starts cache loading and scanning. Source code does not manage its cache, jobs, or browser
state.

## On-disk metadata

Directory-backed collection metadata and annotations are lean text files (see [storage](storage.md)
for formats):

- `metadata.txt` — `DirectorySource` source-root parameters with path-fragment matching.
- `tags.txt` — tag catalog and per-path assignments.

`metadata.txt` belongs to `DirectorySource`; sources without that file simply have no collection
metadata. Annotation files currently live next to the cache, keyed by `source_id`, so a non-filesystem
source still has somewhere to persist user-authored notes/tags/layout. The HDF5 cache itself is
generated and lives outside the source.

## Engineering Rules

This project is pre-pre-alpha. Prioritize simple, principled code and development velocity over
backwards compatibility. Do not keep compatibility paths, fallback behavior, or legacy symbols
unless there is a clear current need. When an old shape is wrong, replace it cleanly and update the
docs and tests to match.

- Julia 1.12; 4-space indentation; ~100-char lines.
- `snake_case` functions/variables, `UpperCamelCase` types, `ALL_CAPS` constants.
- Don't catch errors that should be fixed; let failures surface clearly.
- Detection runs the project's recipes in registration order (first matching `detect` wins) — register specific patterns before general ones.
- Tests live in `test/` with fixture CSVs under `test/fixtures/`.
- Figure-creation tests assert metadata and figure construction behavior, not pixels or labels.

## Deep-dive docs

Current-state reference docs:

| Topic | When to read |
|---|---|
| [api.md](api.md) | You're writing project code: registering item kinds, stats, or plots, or a custom `AbstractDataItem`. |
| [data-model.md](data-model.md) | You're touching source-file records, items, hierarchy traversal, or paths/IDs. |
| [gui.md](gui.md) | You're adding/modifying a panel, window, or interaction. |
| [storage.md](storage.md) | You're adding a new metadata file or changing a format. |
| [cache.md](cache.md) | You're touching HDF5 cache identity, loading, writing, status, or item data. |
| [annotations.md](annotations.md) | You're touching annotations: coordinates, layout, tags, or notes. |

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
