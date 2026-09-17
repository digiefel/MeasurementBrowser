# DataBrowser Roadmap

This roadmap / changelog follows releases of the main `DataBrowser` package.
Checked items already exist on `main`; completed work links to the pull request or commit that
introduced it.

## 0.1.0 — Public API foundation

Establish the package and data-engine foundation on which the interactive application can grow. The
focus is a headless workspace engine, clear package boundaries, public project APIs for ordinary
callbacks and concrete item types, and examples that exercise those APIs without reaching into
workspace internals.

- [x] Replace the old measurement model with source-owned items, workspaces, and hierarchy.
  ([#2](https://github.com/digiefel/MeasurementBrowser/pull/2))
- [x] Consolidate parameters and computed statistics into one metadata pipeline.
  ([#3](https://github.com/digiefel/MeasurementBrowser/pull/3))
- [x] Rebuild the workspace around event-driven publication, a live work graph, and buffered cache
  writes. ([#4](https://github.com/digiefel/MeasurementBrowser/pull/4))
- [x] Watch directory sources continuously and expose cheap workspace/cache status to the GUI.
  ([#5](https://github.com/digiefel/MeasurementBrowser/pull/5))
- [x] Split the monolith into the `DataBrowser` package family.
  ([#6](https://github.com/digiefel/MeasurementBrowser/pull/6))
- [x] Make `DataBrowserGUI` the lightweight GUI host and `DataBrowserPlots` its default Makie
  extension. ([#7](https://github.com/digiefel/MeasurementBrowser/pull/7))
- [x] Move the cache below Core and make tabular cache storage use the Tables.jl interface.
  ([#8](https://github.com/digiefel/MeasurementBrowser/pull/8))
- [x] Add the public documentation book and runnable examples.
  ([`912a3c9`](https://github.com/digiefel/MeasurementBrowser/commit/912a3c90b42f057e7c0096b00d7aa290768d858a))
- [x] Make collections first-class hierarchical values with stable occurrence IDs, compact internal
  keys, cache reopen, collection metadata, and one record/index model shared by the engine and GUI.
- [x] Route registered and typed items through one post-interpretation identity/default
  normalization path. `_registered_item` should adapt callback output into `RegisteredDataItem`, not
  mint final item ids before the common path.
- [x] Remove the unused data-item fingerprint contract and storage. Delete
  `fingerprint(::AbstractDataItem)`, `ItemRecord.item_fingerprint`, and the persisted item-fingerprint
  column; retain source-item fingerprints as the change tokens compared across workspace openings.
- [x] Give source items, data items, and collections the same `id(value)` and `label(value)` public
  interface. Replace `source_item_id`, `source_item_label`, and `item_label` cleanly, without
  compatibility aliases; keep `source_id` and `source_label` for the workspace's data source.
- [x] Add a package-owned `source_item_key::Int64`, persist the source-item id-to-key mapping once,
  and use the integer key for internal and database references just as `item_key` and
  `collection_key` are used. Keep source-item ids as the stable public identity and measure the
  resulting database and loaded-index memory change.
- [x] Resolve `label(value)` once from each live source item, data item, and collection during
  interpretation, persist that display value on the corresponding record, and make
  `label(record)` return it without materializing a payload or running project code in the UI.
- [x] Make typed materialization explicit: reopening restores records, not arbitrary user-defined
  instances; a valid cached processed payload is delivered without rerunning `process`, and an
  `reconstruct(::Type{T}, id, data, metadata)` method rebuilds a concrete item from stored values.
  The default returns data already of type `T`; otherwise it declines reconstruction and the
  engine reruns the required source stages. Collections reconstruct from their stored type, id,
  and metadata. The internal `attach_record` hook remains separate debt tracked below.
- [x] Remove the registration-only `item isa RegisteredDataItem` payload-cache gate, and the
  item-level `cacheable` predicate with it. Persistence becomes the payload's supported shape
  alone; rehydration into a user type becomes `reconstruct` dispatch, keeping cached payload
  delivery distinct from rebuilding a typed item for multiple dispatch. (Moved up from the engine consolidation release: the
  `DataBrowserRecipes` extraction below cannot leave Cache registration-free while this gate
  stands.)
- [x] Replace `HierarchyNode` and the hierarchy-owned object graph with a package-owned
  `CollectionRecord` plus collection-parent, child, and membership indexes owned by
  `WorkspaceIndex`. Collection records expose a durable occurrence ID, compact internal key, and
  resolved label. Treat the browser tree as one projection of that index, not as the data model.
- [x] Replace the monolithic typed `data_items` hook with the staged typed contract of
  [typed-pipeline.md](typed-pipeline.md): explicit `read` and `entries` stages, item and collection
  `process`/`analyze`, project-aware stage forms with context-free defaults, and the source as an
  argument to `read` only. Give the stages distinct work/cache boundaries without a compatibility
  shim.
- [x] Extract the registration dialect into `DataBrowserRecipes`, a package built purely on the
  type API and re-exported by the `DataBrowser` umbrella. Move the recipe-holding `Project` there,
  declare `AbstractProject` in `DataBrowserAPI`, make `register_item!` and
  `register_collection_analysis!` adapter methods over the shared stage contract, and remove every
  registration reference from Core, Cache, Sources, and the GUI packages. Ship `register_csv!` as
  the first premade recipe.
- [x] Clean-up and rename/file organization pass of DataBrowserAPI. `interface.jl` is gone: its
  declarations were the project contract and now sit in `project_contract.jl` beside
  `AbstractProject`, with the stage contract in `stage_contract.jl`.
- [x] Remove the custom internal tracing system; keep workspace diagnostics and use Julia's standard tools for scoped profiling. ([#11])

## 0.2.0 — Engine consolidation

Pause feature expansion for a bounded technical-debt pass. The September 2026 audit found the debt
concentrated in three places: Core reading Cache internals, a duplicated memory cache backend, and
GUI state that mirrors workspace state. Phases are ordered by dependency. Each phase ends with a
green test suite; engine phases also end with a benchmark run compared against the committed
baseline.

Decisions taken during the audit:

- The memory-only cache backend (`MemoryCacheDB`, `AbstractCacheDB`) is removed. `cache=false`
  opens an in-memory DuckDB through the one remaining code path. Granular caching (per stage, per
  kind), non-replayable streaming sources, and disk-budget control are payload *policies* consulted
  by that one backend; they are designed when pipeline inspection (0.14.0) needs them, not before.
- `query_items(workspace, sql)` becomes internal until the query model of 0.5.0 replaces it.
- Two timing systems stay, on purpose: `Browser.@timed` is always-on GUI frame diagnostics;
  `DataBrowserProfiling` is opt-in pipeline profiling.
- Files a user may want to inspect with other tools — tags, notes, layout, later the project file —
  live in the project root. Only the cache lives in the depot (and that may change too).
- The public API has two tiers, both public and both documented: **project authors** use
  `define_project`, `register_*!`, and the workspace operations; **extension authors** use the stage
  contract, the source/item/collection contracts, `reconstruct`, GUI extension hooks, and plot kinds.
  Every exported name belongs to exactly one tier.
- Recipes consumes the same public type API as other extensions. Its convenience must not depend
  on privileged engine access. Keep the type API small and usable on its own; consolidate duplicate
  representations and execution paths rather than adding abstractions to conceal their differences.

### 0.2.1 Repairs

- [x] Run the validated performance benchmark after unit tests and write `bench/status.txt`.
- [ ] Delete dead code: `Workspace.jl` imports of the nonexistent `resolve_type`/`type_name` (two
  precompile warnings); the no-op `reconcile_source_metadata_cache!(…; collections=…)` call in
  `publish_work_success!` and its unused `collections`/`refresh_hierarchy` keywords;
  `rebuild_workspace_hierarchy!`, `cancel_analysis!`, `cancel_cache!`, `reset_work_graph!`,
  `start_cache!`, `stop_cache!`, `_flush_operation`, `_flush_rows`, `_callback_name`,
  `plot_kind_symbol`, `BrowserState.project_locked`, `BrowserState.project_preference`, Core's
  unused `using DataBrowserAnnotations`. Drop unused declared deps (`Statistics`, `Tables` in GUI;
  `DataFrames` in Plots) after verifying.
- [ ] Fix docstrings that describe behavior the code does not have: `set_cache_memory_limit!`
  (workspace) is not live; `load_cache_index` overlays uncommitted buffers; `cache_stage_summary`
  counts queued, not persisted; `query_items` exposes item columns too; `close_cache_db!` ordering;
  `finish_debug_timings!` still records in-flight sections. Remove or create the missing
  `docs/cache.md` and `docs/profiling.md` targets.

### 0.2.2 Cache boundary and one backend

- [x] Unify the cache payload and item reconstruction contract. Core extracts item data; Cache
  stores that payload with the same meaning in memory and on disk. Distinguish a cache miss from
  a stored `nothing`, and check payload availability without loading large tables. Disk writes go
  directly to the payload store, without eligibility predicates or a memory fallback. Always call
  the public reconstruction method with current item metadata, including metadata produced by
  item process and analyze. Preserve the explicit source fallback without running process twice.
- [x] Separate `SOURCE_READ` and `SOURCE_INTERPRET` execution and cache results. Retain read
  results in a bounded memory cache, reuse them during interpretation replay, and distinguish source
  failures, invalidation, and progress by stage. Persist source fingerprints for failed and empty
  interpretations as well as successful ones.
- [x] Finish removing Core's reads of Cache fields. Stage, fingerprint, metadata, and payload
  access uses Cache's exported interface. Core owns its workspace idle deadline; Cache keeps its
  flush-deadline helper private. Remove unused imports, including the concrete cache backend.
- [ ] Let Cache own lifecycle setup: remove no-op start/stop calls and internalize identity/header
  writes. Keep open, close, explicit flush and storage diagnostics meaningful.
- [ ] Delete `MemoryCacheDB` and `AbstractCacheDB`; `cache=false` opens DuckDB `":memory:"`. Removes
  roughly ten duplicated method families whose semantics had already diverged (memory recorded
  result failures in `failures`, disk did not).
- [ ] Give stage outputs exact publication references and preserve metadata by producing stage.
  Separate recorded outcome from payload availability. Reuse valid downstream results after
  upstream eviction; reconstruct resident and restored values with the same metadata.
- [ ] Publish payloads, metadata and completion as one accepted update across all stages. Coordinate
  buffered overlays and transactional flush, check revisions at acceptance, and validate output
  publications on restore. Include metadata-only collection-process output; remove the split
  stage-specific storage/completion calls.
- [ ] Replace one source's output and membership changes through one Cache operation. Core supplies
  replacements and invalidations; Cache retires storage. Cover fewer/empty outputs, moved members,
  shared collections and selection preservation for surviving public IDs.
- [ ] Fix the concurrency issues this exposes: `query_view_signature` is recreated without a lock;
  `cache_pending_counts` locks ten stores separately and returns an incoherent snapshot;
  `open_cache_db` does not stop already-started stores and flush tasks when a later constructor
  fails; `close_cache_db!` keeps only the first exception.

Target interface for the remaining cache work: lifecycle/flush and storage diagnostics; bulk restore;
exact-reference payload reads/checks; complete result writes; source identity, removal and output
replacement; and explicit invalidation supplied by Core. Core normalizes records and extracts
payloads. Cache receives storage data, never a live workspace or user source object. Cache owns
encoding, buffering and persistence, not dependency discovery or recomputation policy. Prepare
expensive encoding outside the workspace lock, then check the existing job revision and accept the
prepared update during publication. Discovery, successful empty interpretation and deletion remain
separate operations. Outcome validity and payload availability remain separate facts.

The order is lifecycle cleanup, one backend, output references, coordinated result publication,
then source replacement. Retain the existing index and scheduler interfaces. Source metadata belongs
to Sources; code fingerprinting belongs to cache identity. Validate storage changes with fixed
workloads, completed disk writes, cache reads and process memory measurements.

### 0.2.3 One index, one status

- [ ] Evaluate the work graph against cancellation, invalidation, priority, streaming, and
  collection edge cases; finish with a bounded tuning pass or an explicit redesign.
- [ ] Remove `WorkspaceIndex.source`. It is a full `SourceScan` snapshot rebuilt by
  `refresh_workspace_source!` (copies the whole collection index and sorts every item per batch),
  and its only two readers ask `isa SourceScan`. Derive that boolean from the scan state.
  `SourceScan` remains the Cache load result only.
- [ ] Add `WorkspaceIndex(source_id)`; delete the two hand-built seven-argument constructions.
- [ ] Name item metadata layers after the pipeline stages on `ItemRecord` and remove the parallel
  `WorkspaceIndex.item_metadata` dict, so each item's interpret and analyze layers live in one place.
- [ ] One status model. `scan.state` (7 symbols), `cache_state` (8), `cache.operation` (4), and the
  `status.label in ("Fresh", "Loaded", "Errors")` check in `workspace_status` collapse into one
  enum and one `busy` derivation; `source_scan_running`, `cache_work_running`,
  `engine_work_running`, `workspace_busy` become one function. `analysis_errors` becomes a typed
  collection, not `Dict{Union{String,Int64},String}`.
- [ ] Small Core cleanups: `indexed_collection_path` and `collection_value` drop the `Workspace`
  argument they ignore; the `ItemRecord` copy constructor stops `deepcopy`ing metadata by default;
  the scan's `current::Dict{String,Any}` fingerprint map gets a type; `InspectorTable` and
  `merge_item_tables` move out of Core (only GUI and Plots use them).
- [ ] Use compact integer item keys in SQL tables and other measured hot paths while retaining
  stable logical item identities at the project boundary.

### 0.2.4 Public API tiers

- [ ] Extend the common reconstruction contract to collection analysis metadata. Collection
  reconstruction currently receives only `own_metadata`, without analysis results. User types
  decide which results to retain and how to represent them. Remove the registration adapter's
  reliance on the internal `attach_record`/`ItemRecord` path for labels and collection paths;
  it must consume the same public contract as any other extension. Item payload and metadata
  delivery belongs to 0.2.2.
- [ ] Preserve member metadata changes returned by collection `process`, including changes with
  unchanged payloads. `run_collection_process` currently selects outputs only when `item_data`
  changes by identity, and collection completion publishes no member metadata. Cover cumulative
  device history: successive measurements need different wakeup/fatigue counts even when their
  waveform payloads are unchanged.
- [ ] Sort every exported name into the project-author or extension-author tier and stop exporting
  the rest: `items_for_file`, `SourceFile`, the concrete `Project`, `gui_timings`,
  `reset_timings!`, `plot_kinds`, `CollectionRecipe`, `ItemRecipe`, `NoMatch`,
  `RegisteredReadResult`. Export what extension authors need and cannot reach today:
  `GuiExtension`, `register_gui_extension!`, `NamedCollection`.
- [ ] Generate the umbrella's `using`/`export` lists from one table instead of two hand-kept blocks.
- [ ] One selector concept — `ItemRecord`, id, or current selection — implemented once and shared by
  `select_items!`, `materialize_items`, `read_item_data`; delete the runtime-typed
  `select_items!(::AbstractVector)` fallback. The Cache-side `read_payload` name is handled by
  the payload contract work in 0.2.2.
  `query_items(ws)` becomes `item_ids(ws)`; `query_items(ws, sql)` becomes internal.
- [ ] Add the workspace accessors the GUI currently reaches for by field: `items(ws)`,
  `item(ws, id)`, `collections(ws)`, `selection(ws)`, `errors(ws)`, and an `on_change(ws, f)`
  subscription replacing the `status_dirty` poll.
- [ ] Recipes: rename `register_collection_analysis!` to `register_collection!` (it registers
  `process` too); look the recipe up once per stage chain instead of three linear scans.
- [ ] Strengthen internal module boundaries and import hygiene across the package family.
- [ ] Run every example entirely through the documented public APIs and remove any remaining public
  callback dependency on cache, index, scheduler, or browser values.

### 0.2.5 GUI as one caller of the API

- [ ] Move the 58 `workspace.<field>` reads across ten GUI and Plots files onto the Phase 3
  functions. This is the first concrete step toward the shared command layer of 0.6.0.
- [ ] One shared helper for selected-item table materialization; `TableInspector.jl:19-106` and
  `TablePlotPanel.jl:67-113` are the same code.
- [ ] One hierarchy projection per frame. `_render_hierarchy_tree_panel` is 347 lines and walks the
  collection tree five times per frame; cache the prepared projection between invalidations and
  rebuild item-panel rows only when selection, visibility, tags, or item state change.
- [ ] Remove duplicated GUI state: filters held both as Julia strings and ImGui pointers and synced
  every frame; the main plot window special-cased instead of being a `PlotViewState`;
  `PlotState.kind_by_item` beside `PlotViewState.plot_kind`; `PROJECT_PLOT_RECIPES` with no removal
  path; `PlotsExtension.reset!` not clearing `MAKIE_CONTEXT`.
- [ ] Tighten the extension boundary: `PlotsExtension` calls `Browser._project_visible_selection`
  and `Browser._items_for_ids`; make those tier-B functions. Reduce the twelve-hook protocol to
  what the one extension uses; type `save_view`/`load_view!` instead of `Dict{String,Any}`.
- [ ] Split `BrowserState` (40 fields) into selection/filter, windows, diagnostics, persistence,
  and lifecycle state. Last, because the items above shrink it first.

### 0.2.6 Cache internals

- [ ] One `MetaVType` registry replacing the six parallel type maps (`_meta_vtype`,
  `_vtype_julia`, `_meta_vtype_sql`, `_encode_wide_value`, `_decode_wide_value`,
  `_duckdb_sql_type_maybe`).
- [ ] A `with_connection(db) do … end` helper replacing thirteen hand-written connect/try/finally
  sites.
- [ ] Group `CacheDB` (25 fields) into index stores, metadata stores, payload store, and runtime
  state. Inline `load_cache_index_body` and `_report_loaded_cache_index`.
- [x] Use one `PipelineStage` for scheduling and cache results, ordered from source read through
  collection analyze. Source read and interpretation have independent outcomes.
- [ ] Settle on one of `cache`/`cachedb`/`cache_db` during the interface changes.
- [ ] Cache identity gains a project fingerprint. Vision §12 requires project definition + data +
  parameters; today it is project *name* + source id, so editing project code reuses stale results.
  Changed code must invalidate success and failure automatically. Investigate Revise integration
  and a stable persisted code/environment fingerprint; begin conservatively at project scope.
  Mismatch follows eager/requested recomputation policy, not a warning-only path.
- [ ] Define and test persistence at every expensive pipeline boundary (discovery, read, entries,
  item and collection process/analyze). A valid persisted stage satisfies downstream work without
  rerunning earlier user code.
- [ ] Support expensive or non-repeatable sources (simulations, compressed inputs, streams) through
  source-owned durable handles or persisted interpreted outputs; add `replayable(::AbstractDataSource)`
  when the first such source arrives.
- [ ] Profile the DuckDB flush path at millions of rows and remove the dominant avoidable cost.
- [ ] Document cache contracts in source docstrings and the generated public documentation.

### 0.2.7 Sources, Annotations, Recipes

- [ ] Audit source fingerprinting and document exactly what each change token invalidates across
  live updates and workspace reopen.
- [ ] Keep `default_collection_path(source, source_item)` source-owned. Replace
  `annotate_collection_path` with `source_collection_metadata(source, paths)`: a batch of public
  collection-identity paths in, owned metadata dictionaries per path level out. DirectorySource
  owns sidecar parsing and matching; Core owns composition and invalidation. Store declared and
  source metadata separately and replace the source contribution on interpretation, refresh and
  reopen, so removed keys reveal the declared values beneath them. Snapshot source metadata under
  its lock without calling user reconstruction there. Test both DirectorySource and another source.
- [ ] Split `DirectorySource` into scan, `metadata.txt`, and watcher files. The watcher rescans
  the whole tree and fingerprints every file three times on each event; `.git` is skipped by the
  watcher but traversed by the scan; `readdir` order is not sorted; `metadata_lock` is held while
  calling user `reconstruct`; two concurrent `watch_source` calls race and a crashed watcher task is
  reported as "already watching".
- [ ] `metadata.txt`: decide whether values may be quoted; validate metadata value types at parse
  time instead of failing downstream.
- [ ] Annotations live in the project root. The GUI currently derives the annotation root from the
  cache path under `DEPOT_PATH` (`Browser/Operations.jl:23-26`), against the vision. The GUI uses
  `item_annotation_key`/`ancestor_annotation_keys` instead of rebuilding keys; enforce or drop the
  "item ids and collection paths never overlap" claim; one parse-error type for the three stores.
- [ ] Profiling: delete or test the sampling profiler API (`start_sampling!`, `stop_sampling!`,
  `cancel_sampling!`; no callers).

### 0.2.8 Tests

- [x] Package-owned contract suites using only each package's dependencies; individually selectable
  files, with unchanged successful dependency closures reused.
- [x] Full verification includes all subpackages and the umbrella benchmark as an application smoke
  test. Commands live in `test/README.md`; measurements are defined in the benchmark code.

### 0.2.9 Benchmarks

- [x] Measure clean package precompilation when source dependencies change; reuse its compiled
  output for tests and application workloads.
- [x] Measure process-to-browser startup and saved-cache reopen with a real GUI and default plots.
- [x] Measure engine indexing overhead, cache throughput, materialization latency and process memory
  using prepared type-API inputs.

### 0.2.10 Documentation

- [ ] `ARCHITECTURE.md` promises two diagrams per package and has them for Sources only. Add Core
  and Cache after Phases 1–2, GUI after Phase 4.
- [ ] `docs/api.md`: the two tiers from Phase 3, one line per name.
- [ ] `vision.md` links to `cache.md`, `plans/plotting-api-design.md`, `plans/spatial-browser.md`,
  `plans/project-persistence.md`; none exist. Fix or remove.
- [ ] `AGENTS.md`: require a bench smoke before committing engine changes.

## 0.3.0 — Tags

Make tags a dependable, machine-interpretable way to classify items and collections. Start by
checking the current `tags.txt` loading path and the state of existing files, then complete the API
and GUI around arbitrary tags. The existing `bad` behavior becomes one ordinary use of the same tag
system.
Every change keeps `tags.txt` human-readable and recoverable independently of the application.

- [x] Implement plain-text tag definitions, item/collection assignments, inherited lookup, and the
  current `bad` bridge. ([#2](https://github.com/digiefel/MeasurementBrowser/pull/2))
- [ ] Provide API operations to create, edit, remove, assign, and unassign tags on items and
  collections.
- [ ] Provide GUI controls for the same operations, including multi-selection.
- [ ] Apply tags consistently to colours, visibility, and the selected item set.

## 0.4.0 — Notes

Make notes a complete human-facing memory and context feature, separate from tags and plotting.
Notes belong to items and collections, remain readable as plain text, and are edited primarily
through the API and GUI.

- [x] Implement the current plain-text note sections and their read/write API.
  ([#2](https://github.com/digiefel/MeasurementBrowser/pull/2))
- [ ] Settle how notes attach to current items and collections.
- [ ] Add API operations for reading, creating, editing, and removing notes.
- [ ] Add GUI views for reading and editing notes in the context of the current item or collection.
- [ ] Preserve notes across workspace reopen and source refresh.

## 0.5.0 — Find, filter, and view items

Turn the current hierarchy-only browser into several coordinated views over the same item set. Julia
code gets concrete database queries; GUI filtering produces live selections without exposing query
language. Clicking, Ctrl-clicking, filtering, grouping, and switching views all operate on the same
selection model.

- [x] Query committed effective metadata through DuckDB.
  ([#4](https://github.com/digiefel/MeasurementBrowser/pull/4))
- [ ] Define one query model over metadata and tags, with concrete results for Julia callers and live
  results for GUI views.
- [ ] Add a visual filter builder to the GUI.
- [ ] Add a flat item table with selectable metadata columns alongside the hierarchy view.
- [ ] Support flattening, grouping, sorting, and filtering without copying item state into the GUI.
- [ ] Show when computed statistics used by a filter are still being populated.
- [ ] Persist useful item-view and filter state with the project.

## 0.6.0 — Shared application API

Give Julia code and the GUI the same application capabilities through ordinary public functions.
The API covers operations on projects, workspaces, selections, queries, views, and background work;
the GUI becomes one caller of those operations rather than a second implementation.

- [ ] Define the public operations needed to drive a workspace from Julia without browser state.
- [ ] Move GUI actions onto those operations and keep GUI-only state limited to rendering and local
  interaction.
- [ ] Make repeated project setup, workspace operations, and view changes idempotent where users
  naturally rerun code during development.
- [ ] Expose progress, failures, cancellation, and results in forms usable by both interactive and
  programmatic callers.
- [ ] Keep the REPL usable while workspaces and GUI windows remain open.

## 0.7.0 — Figure composer

Build the DBPlots figure composer before the interactive plot builder. Its first inputs are axes and
figures created by project code. Users can arrange those existing components into larger figures and
edit the composition through both Julia and the GUI.

- [ ] Define composable DBPlots objects for axes, panels, and complete figures.
- [ ] Accept project-created Makie axes and figures without taking ownership away from project code.
- [ ] Add GUI operations for creating layouts, inserting components, moving them, resizing them, and
  removing them.
- [ ] Support linked axes and shared presentation settings across composed panels.
- [ ] Keep composed figures live against their item selections where requested.

## 0.8.0 — Plot builder

Add the fast path from selected data to an individual axis. The plot builder is a DBPlots visualizer
with a focused toolbar: choose data, map columns or dimensions, select a plot form, and adjust its
presentation. Every axis it creates can be inserted into the figure composer.

- [ ] Build axes from the current item selection and from filtered item views.
- [ ] Provide data mapping controls for common tabular X/Y plots.
- [ ] Provide line and scatter plots with editable series, axes, labels, scales, and styling.
- [ ] Represent builder state through the shared application API rather than private widget state.
- [ ] Insert a built axis into an existing composed figure without recreating it by hand.

## 0.9.0 — Python API

Make DataBrowser usable as a Python package over JuliaCall. Python is a full application interface:
it manages or attaches to the Julia environment, calls the shared application API, and presents data
through natural Python containers.

- [ ] Add an installable Python package in this repository with Julia and package bootstrapping.
- [ ] Wrap projects, workspaces, items, selections, queries, tags, notes, and background operations
  with Python-facing types and errors.
- [ ] Convert tabular and array data naturally to pandas and NumPy without unnecessary copies where
  the runtimes permit it.
- [ ] Support Python readers, processors, and analyzers through explicit callback bridges.
- [ ] Drive figure composition and the plot builder from Python.
- [ ] Prototype the Python-native visualization path and settle the roles of Matplotlib, lightweight
  ImPlot views, and Julia-owned DBPlots visualizers.

## 0.10.0 — Figure annotations

Add figure annotations as DBPlots-owned, editable objects inside composed figures. They describe a
figure rather than an item: arrows, text, regions, fit labels, and other axis-relative additions.

- [ ] Define figure-annotation objects and their coordinate systems.
- [ ] Add, select, edit, move, style, and remove figure annotations through Julia and the GUI.
- [ ] Attach annotations to axes, data coordinates, or layout coordinates as appropriate.
- [ ] Integrate fit results and region selections without coupling them to item tags or notes.
- [ ] Preserve figure annotations with composed-figure state.

## 0.11.0 — Core visualizers

Make common inspection fast without requiring Makie or project-specific plot code. Expand the base
GUI with lightweight visualizers, using ImPlot where it provides the right interaction and keeping
the table inspector as a first-class view.

- [ ] Add a lightweight line/scatter visualizer for simple numeric tables and vectors.
- [ ] Add concise scalar, metadata, and nested-value inspectors.
- [ ] Let users open several independent visualizer windows over different live selections.
- [ ] Register every built-in visualizer through the same GUI extension and window surfaces used by
  other first-party features.

## 0.12.0 — Array and image visualizers

Add dedicated exploration for multidimensional data. These visualizers use the plot builder and
figure composer where useful while retaining controls for dimensions, slicing, colour, and image
presentation.

- [ ] Add two-dimensional array and image visualizers.
- [ ] Persist supported dense N-dimensional array payloads and verify warm delivery without rereading
  or reinterpreting their source items.
- [ ] Add heatmaps for gridded arrays and suitable tables.
- [ ] Add dimension and slice controls for higher-dimensional arrays.
- [ ] Add editable colour limits, scales, colormaps, and image display settings.
- [ ] Offer every applicable visualizer for a value rather than assigning each value one exclusive
  shape.

## 0.13.0 — Comparison and summary visualizers

Add visualizers for comparing many items and groups. This release builds on filtered item views and
the plot builder so overlays and summaries remain attached to the selections that created them.

- [ ] Add enhanced table plots and multi-trace overlays.
- [ ] Group and style traces by metadata, tags, and collection.
- [ ] Add parameter-sweep and pivoted-table views.
- [ ] Add summaries and histograms across selections and groups.
- [ ] Add fit views for common models, beginning with linear fits.

## 0.14.0 — Pipeline inspection

Let users inspect the values produced throughout a project pipeline. The workspace supplies the
chosen value through its normal data-access path; inspectors and visualizers continue to consume the
payload they are given rather than learning special pipeline logic.

- [ ] Expose the results of `read`, `entries`, `process`, `analyze`, collection `process`, and
  collection `analyze` through the shared application API.
- [ ] Use the public API vocabulary for the GUI control and skip operations absent from a project.
- [ ] Coalesce sibling items when moving to a source-level result and group items when moving to a
  collection-level result.
- [ ] Show metadata added by analysis operations beside the corresponding data.
- [ ] Reuse cached results and retain uncached intermediate values within explicit memory bounds.
- [ ] Add the global GUI control and make the current inspectors and visualizers respond to the
  selected pipeline result without per-visualizer pipeline code.

## 0.15.0 — Project configuration

Consolidate the existing project-local state into `dbproject.toml`. Entry code remains responsible
for starting the project: Julia or Python code loads or creates the config, defines the project, and
opens the workspace. Reopening means rerunning the same entry code.

- [x] Persist tree, item, filter, and plot-view state in the current `databrowser.toml`.
  ([#2](https://github.com/digiefel/MeasurementBrowser/pull/2),
  [#7](https://github.com/digiefel/MeasurementBrowser/pull/7))
- [ ] Define the TOML schema for sources, cache settings, GUI views, and extension state.
- [ ] Load or create the config from Julia and Python entry code.
- [ ] Edit project settings through the API and GUI and write them back predictably.
- [ ] Place `tags.txt` and the notes file beside `dbproject.toml`.
- [ ] Restore useful GUI and visualizer state when the same entry code reopens the project.
- [ ] Keep cache reuse and source changes correct across repeated openings.

## 0.16.0 — Command-line interface

Add a real shell interface over the shared application API. It runs project entry code, obtains the
workspace that code creates, and performs the same data and project operations available from Julia,
Python, and the GUI. Startup strategy is chosen from measured Julia startup and reuse behavior.

- [ ] Define how a CLI invocation identifies and runs project entry code.
- [ ] Add commands for inspecting projects and workspaces, querying and selecting items, running
  processing, and exporting data or figures.
- [ ] Provide structured output suitable for shell pipelines alongside readable interactive output.
- [ ] Measure startup cost and implement the simplest process-reuse or precompilation strategy that
  makes repeated commands practical.
- [ ] Keep long-running workspace work observable and cancellable from the shell.

## 0.17.0 — Integration and release preparation

Exercise the complete application as one system, keep performance visible, and make installation and
first use straightforward. This version is for repairing the seams found when the Julia API, GUI,
Python API, CLI, project configuration, and visual tools are used together on real projects.

- [ ] Install DataBrowser as an ordinary Julia package without hand-assembling its component
  environment.
- [ ] Run the documented examples and at least one substantial real project through the supported
  interfaces.
- [ ] Validate first-party add-on packages can provide types, processing, visualizers, windows, and
  project scaffolding while project environments pin them normally.
- [ ] Keep browse-while-building, warm reopen, item views, filtering, plotting, and figure editing
  responsive at realistic scales.
- [ ] Complete user documentation for installation, project setup, GUI use, Julia, Python, and the
  CLI.
- [ ] Resolve the cross-interface inconsistencies and reliability failures found by that use.

## 1.0.0 — Complete DataBrowser application

DataBrowser 1.0 is a cohesive scientific data application built around live Julia projects. A user
can install the package, define a project in ordinary code, browse and organize its items, inspect
the pipeline, query and filter metadata, process data, build and compose figures, and return to the
same project state. The same application capabilities are available through Julia, the GUI, Python,
and the command line.

- [ ] Tags and notes are dependable project features with complete API and GUI use.
- [ ] Queries and visual item views support concrete programmatic results and live GUI selections.
- [ ] The workspace engine remains responsive and observable while sources and project code change.
- [ ] The figure composer, plot builder, figure annotations, and generic visualizers cover common
  interactive analysis without project-specific GUI code.
- [ ] Pipeline results are inspectable through the API and GUI using the project's own operation
  vocabulary.
- [ ] Julia, GUI, Python, and CLI callers share the same application capabilities.
- [ ] Entry code can create or reopen `dbproject.toml`, its project state, tags, notes, cache, and
  saved views.
- [ ] Installation, documentation, examples, diagnostics, and realistic performance checks support
  normal use of the application.

## Longer-term goals

### Packaged application and shareable project bundles

Provide a packaged application that can open a reproducible project bundle containing its Julia
environment, entry code, configuration, data, tags, notes, and optional cache. The bundle becomes the
shareable unit, gives the packaged application one clear entry point, and may instantiate its Julia
environment through Pkg when opened.

### GUI-authored workflows

Represent a sequence of GUI actions as an editable, replayable equivalent of a script: open a
project, find data, process it, create views, edit figures, and export results.

### First-party analysis add-ons

Build substantial XPS, ellipsometry, semiconductor, ferroelectric, and other analysis add-ons from
the same extension surfaces. An add-on supplies reusable code and may provide project scaffolding
while each project keeps its own configuration and data.
