# DataBrowser Architecture

> **Onboarding doc.** Start here before touching code. This file describes the architectural model
> that is not obvious from the filesystem. It should not be a file inventory.

## What this is

DataBrowser is a Julia working environment for data projects. It should be strictly
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
source item → interpret → logical data → process → analyze → collection process/analyze → views
                 │             │             │
                 └─ index      └─ DuckDB     └─ DuckDB + item metadata
```

The engine is written against the **low-level source contract** ([api.md](api.md)): an
`AbstractDataSource` owns lifecycle and discovery, an `AbstractDataSourceItem` is one discovered unit,
and an `AbstractDataItem` is one logical browsable item. `scan_source!` discovers source items
(streaming via `on_item` when the source supports it), compares each fingerprint with the cache
fingerprint table, and queues new or stale items during the walk without publishing yet — so
directory discovery stays fast. After the walk, a publisher task submits `SourceChanges` batches so
`SOURCE_INTERPRET` starts immediately and overlaps the remaining enqueue without holding the walk
on `publish_lock` / ingest. Removals and source-metadata reconcile happen in a final batch when
that publisher drains. Source watchers submit the same batch type;
`DirectorySource` watches the directory tree recursively and reports file creates, edits, renames,
deletes, and metadata changes through that path. Interpretation workers call
`data_items(project, source, source_item)`, put interpreted data into the memory cache, and publish
their own completions: each finishing worker takes the workspace publish lock, rejects stale
revisions, atomically publishes replacement records into `WorkspaceIndex`, sends semantic record
writes/deletes to `ProjectCache`, queues follow-up work, and notifies the idle condition. Nothing
polls: the engine progresses on completion events alone, blocking callers
(`request_processed_items`, `wait_workspace_idle!`) sleep until a publication wakes them, and the
GUI frame only calls `refresh_status!` to rebuild the display snapshot. The graph holds only live
jobs (`:waiting`, `:queued`, `:running`): each node tracks who waits on it (`dependents`) and how
many dependencies remain (`pending`). Finished work lives in the cache `result_states` ledger (or
interpret/source tables for `SOURCE_INTERPRET`). `WorkspaceIndex` owns completed records,
hierarchy, metadata, failures, and selections. The five work kinds — `SOURCE_INTERPRET`,
`ITEM_PROCESS`, `ITEM_ANALYZE`, `COLLECTION_PROCESS`, `COLLECTION_ANALYZE` — run on one fixed
long-lived pool. Processing is selection-driven by default: a live scan interprets every changed
source item (so the tree fills) but only interprets — `ITEM_PROCESS`, `ITEM_ANALYZE`, and collection
work are scheduled on demand when items are requested (`request_processed_items`). `background_processing`
opts into eagerly processing and analyzing every item as it is interpreted; background and selected
work then share the same work node, and selection only raises priority and joins the existing
revision. Either way, interpret-time invalidation always clears stale `result_states` rows so the
next read recomputes. Memory-only cache sessions never schedule unselected background processing on
reopen gap-fill (they have no persisted index to fill from).
Collection-node work runs afterward from published item-analysis results.
The cache ([cache.md](cache.md)) restores the previous hierarchy quickly while scanning continues.

When a view needs item data, delivery consults the live graph then the cache: a live node blocks
until it publishes; a `RESULT_READY` row loads the payload; absence enqueues work. Otherwise the
processing job loads interpreted data from the in-memory interpreted store or reuses the normal
`data_items` path as source fallback, then runs `process`.
Registration views receive processed data with effective metadata. Type-based views receive the
project's concrete `AbstractDataItem` values.
GUI selection, background work, and Makie embedding are described in [gui.md](gui.md).

The **high-level callback API** (`define_project` + `register_*`, the exported convenience surface)
uses the built-in `DirectorySource <: AbstractDataSource`: the source walks a data root into
`SourceFile`s, while the registration adapter applies `detect`/`read`, resolves the returned
`DataItem` descriptions, then applies `process`/`analyze` through the same engine. Registered callbacks operate
on ordinary data and metadata. A hand-written project/source pair returns concrete
`AbstractDataItem` values and keeps those types through processing and visualization.

## Ownership Boundaries

The most important boundary is between source meaning and package machinery.

A source owns:

- discovering source items
- its own lifecycle and change watching
- source-specific vocabulary, such as the noun used for discovered source items
- source-specific parameter input, such as `DirectorySource` loading `metadata.txt`

A project/source implementation owns:

- interpreting each source item into logical data items
- processing one interpreted item
- computing per-item and per-collection metadata (item/collection `analyze`, collection `process`)
- defining project-specific visualizers when generic ones are not enough

The workspace owns:

- the open source and its identity
- the progressively populated item index
- selection identities
- cache identity, freshness, storage, and repair
- scanning, cache work, progress, errors, and cancellation
- work dependency graph state and source fallback

The browser owns windows, controls, filters, and temporary rendering state. Annotations store
user-authored tags, notes, coordinates, and spatial positions (currently next to the cache, keyed by
`source_id`; see [cache.md](cache.md)). Other package modules own generic visualizers, workflow
persistence, and figure composition.

The Performance window always exposes bounded project callback timings, plot timings, failures,
memory, and rendering diagnostics. A workspace opened with internal profiling enabled also owns a
manual structured trace session for engine spans, process counters, and optional CPU samples. Starting
that session during a build first cancels and drains the current work, then records one clean rebuild.
See [profiling.md](profiling.md). No diagnostic is persisted in the generated cache.

Source code should not know whether data came from memory, cache, or the origin. Package code should
not know the meaning of a source item beyond the contract methods it calls.

The package expresses that boundary through focused packages. `DataBrowserAPI` defines the source
and item contracts plus the registration-based `Project`. `DataBrowserSources` owns the built-in
directory source, `SourceFile`, file fingerprints, directory traversal, sidecar exclusion, and
`metadata.txt` collection input. Its registration adapter is one high-level API over the same source
contract used by typed projects. `DataBrowserCore` owns interpretation, hierarchy construction,
workspace state, selection, and background work. `DataBrowserCache` owns generated DuckDB state.
`DataBrowserPlots` owns plotting operations and Makie integration. `DataBrowserGUI` owns typed
frontend state and CImGui rendering. `DataBrowser` exports the public API and wires the default
packages together.

`open_workspace(project, source)` or `open_workspace(project, root_path)` creates that owner and
opens the source, starts cache loading and scanning, and attaches its watcher. Source code does not
manage its cache, jobs, or browser state.

## On-disk metadata

Directory-backed collection metadata and annotations are lean text files (see [storage](storage.md)
for formats):

- `metadata.txt` — `DirectorySource` source-root collection metadata with path-fragment matching,
  projected through the package-owned collection values created by the registration adapter.
- `tags.txt` — tag catalog and per-path assignments.

`metadata.txt` belongs to `DirectorySource`; sources without that file simply create registered
collections with empty metadata. Annotation files currently live next to the cache, keyed by `source_id`, so a non-filesystem
source still has somewhere to persist user-authored notes/tags/layout. The DuckDB cache itself is
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
| [api.md](api.md) | You're writing project code with registration callbacks or typed sources/items. |
| [data-model.md](data-model.md) | You're working with identity, metadata, collections, or pipeline semantics. |
| [gui.md](gui.md) | You're adding/modifying a panel, window, or interaction. |
| [storage.md](storage.md) | You're adding a new metadata file or changing a format. |
| [cache.md](cache.md) | You're touching DuckDB cache identity, loading, writing, status, or item data. |
| [profiling.md](profiling.md) | You're measuring project callbacks, engine spans, build counters, CPU hotspots, or crash traces. |
| [annotations.md](annotations.md) | You're touching annotations: coordinates, layout, tags, or notes. |

Future-facing docs and roadmaps live under [plans/](plans/), separate from the current-state docs
above. Start with [vision.md](vision.md) for the product and architecture north star,
[plans/roadmap.md](plans/roadmap.md) for the staged plan, and [plans/README.md](plans/README.md) for
the plan map.

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
