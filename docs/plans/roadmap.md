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
- [ ] Route registered and typed items through one post-interpretation identity/default
  normalization path. `_registered_item` should adapt callback output into `RegisteredDataItem`, not
  mint final item ids before the common path.
- [ ] Remove the unused data-item fingerprint contract and storage. Delete
  `fingerprint(::AbstractDataItem)`, `ItemRecord.item_fingerprint`, and the persisted item-fingerprint
  column; retain source-item fingerprints as the change tokens compared across workspace openings.
- [ ] Give source items, data items, and collections the same `id(value)` and `label(value)` public
  interface. Replace `source_item_id`, `source_item_label`, and `item_label` cleanly, without
  compatibility aliases; keep `source_id` and `source_label` for the workspace's data source.
- [ ] Add a package-owned `source_item_key::Int64`, persist the source-item id-to-key mapping once,
  and use the integer key for internal and database references just as `item_key` and
  `collection_key` are used. Keep source-item ids as the stable public identity and measure the
  resulting database and loaded-index memory change.
- [ ] Resolve `label(value)` once from each live source item, data item, and collection during
  interpretation, persist that display value on the corresponding record, and make
  `label(record)` return it without materializing a payload or running project code in the UI.
- [ ] Make typed materialization explicit and tested: reopening restores records, not arbitrary
  user-defined instances; a valid cached processed payload is delivered without rerunning
  `process`, an opt-in `construct(::Type{T}, data, metadata)` method rebuilds a concrete item from
  record and payload with its rederived identity validated, and a missing typed item without one
  is recreated through `source_items` → `read` → `entries`. Obtain a typed collection value from
  `collection(item)` only after materializing its owning item.
- [x] Replace `HierarchyNode` and the hierarchy-owned object graph with a package-owned
  `CollectionRecord` plus collection-parent, child, and membership indexes owned by
  `WorkspaceIndex`. Collection records expose a durable occurrence ID, compact internal key, and
  resolved label. Treat the browser tree as one projection of that index, not as the data model.
- [ ] Replace the monolithic typed `data_items` hook with the staged typed contract of
  [typed-pipeline.md](typed-pipeline.md): explicit `read` and `entries` stages, item and collection
  `process`/`analyze`, project-aware stage forms with context-free defaults, and the source as an
  argument to `read` only. Give the stages distinct work/cache boundaries without a compatibility
  shim.
- [ ] Extract the registration dialect into `DataBrowserRecipes`, a package built purely on the
  type API and re-exported by the `DataBrowser` umbrella. Move the recipe-holding `Project` there,
  declare `AbstractProject` in `DataBrowserAPI`, make `register_item!` and
  `register_collection_analysis!` adapter methods over the shared stage contract, and remove every
  registration reference from Core, Cache, Sources, and the GUI packages. Ship `register_csv!` as
  the first premade recipe.
- [ ] Clean-up and rename/file organization pass of DataBrowserAPI. // e.g. what's interface.jl?? 
- [ ] Remove built-in profiling and consolidate/document proper debugging and profiling.
- [ ] Run every example entirely through the documented public APIs and remove any remaining public
  callback dependency on cache, index, scheduler, or browser values.

## 0.2.0 — Tags

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

## 0.3.0 — Notes

Make notes a complete human-facing memory and context feature, separate from tags and plotting.
Notes belong to items and collections, remain readable as plain text, and are edited primarily
through the API and GUI.

- [x] Implement the current plain-text note sections and their read/write API.
  ([#2](https://github.com/digiefel/MeasurementBrowser/pull/2))
- [ ] Settle how notes attach to current items and collections.
- [ ] Add API operations for reading, creating, editing, and removing notes.
- [ ] Add GUI views for reading and editing notes in the context of the current item or collection.
- [ ] Preserve notes across workspace reopen and source refresh.

## 0.4.0 — Find, filter, and view items

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

## 0.5.0 — Engine consolidation and scale

Pause feature expansion for a bounded architecture and performance pass. Measure current behavior,
make explicit decisions about the work graph and layer boundaries, and remove recurring per-frame or
per-row costs before the application API and plotting surface grow substantially.

- [ ] Evaluate the current work graph against cancellation, invalidation, priority, streaming, and
  collection edge cases; finish with either a bounded tuning pass or an explicit redesign.
- [ ] Clarify ownership between the workspace, index, project cache, database, and write buffers;
  remove duplicated state and layer-skipping call paths.
- [ ] Use compact integer item keys in SQL tables and other measured hot paths while retaining stable
  logical item identities at the project boundary.
- [ ] Audit source fingerprinting and document exactly what each source-provided change token
  invalidates across live updates and workspace reopen.
- [ ] Remove the registration-only `item isa RegisteredDataItem` payload-cache gate. Choose disk
  caching through `cacheable(stage, value)` dispatch on the stage's output value and the payload's
  supported shape, while keeping cached payload delivery distinct from reconstructing a
  user-defined typed item through `construct` for multiple dispatch.
- [ ] Define and test persistence at every expensive pipeline boundary: source discovery, read,
  entries, item processing and analysis, and collection processing and analysis. A valid persisted
  stage must satisfy downstream work without rerunning earlier user code; supported core payload
  shapes cache automatically, while custom typed values may opt into rehydration with a
  `construct` method.
- [ ] Provide the workspace-keyed cached constructor entry point
  `(::Type{T})(ws, key...) where {T<:AbstractDataItem}` so concrete user items rehydrate from the
  cache on demand — threaded, query- and selection-driven — and settle its key grammar.
- [ ] Support expensive or non-repeatable sources such as simulations, compressed inputs, and
  streams through source-owned durable handles/snapshots or persisted interpreted outputs. Memory
  eviction must not rerun a simulation or consume a stream again when a valid durable result exists.
- [ ] Collapse workspace busy and idle decisions into one engine model and one watcher snapshot.
- [ ] Rebuild item-panel rows only when selection, visibility, tags, or relevant item state changes.
- [ ] Cache hierarchy preparation and visible-collection results between invalidations.
- [ ] Profile the DuckDB flush path at millions of rows and remove the dominant avoidable cost.
- [ ] Strengthen internal Julia module boundaries and explicit import hygiene inside the package
  family.

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
