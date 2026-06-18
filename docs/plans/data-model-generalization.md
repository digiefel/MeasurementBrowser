# DataBrowser: General Data Model

## Goal

Grow MeasurementBrowser into **DataBrowser**: an open, scriptable, data-focused IDE rather than a
measurement-only viewer. The same app should browse measurements, images, fitted models, and
arbitrary datasets; organize them by meaning rather than by folder; and present each through the
visualizers and figures in [workspace-vision.md](workspace-vision.md).

The data model we ship today is the right shape with the wrong vocabulary and two hard couplings.
Almost every field of `MeasurementInfo` is already domain-neutral. The parts that block other kinds
of data are exactly two:

1. the data is pinned to `DataFrame`, and
2. an item's place in the tree is a stored, scan-time field (`DeviceInfo.location`) rather than a
   derived view.

This doc defines how the core types generalize while keeping the concreteness that makes the current
model legible.

---

## Principles

These three decide every ambiguous call below.

1. **Generality at the base, clarity at the leaves.** The abstract umbrella type `AbstractDataItem`
   (and `Collection`) is allowed to be bland — *no one reads it in normal use*. Its meaning is a small
   documented interface, not the noun. The clarity lives in the concrete types **projects** define
   for their own domain (an RuO2 measurement, a device) — or in well-chosen `kind`s and `parameters`
   on the package's normal `DataItem`. The package ships `AbstractDataItem`, the interface, and the
   concrete `DataItem`; it ships **no** domain types.
2. **Two first-party APIs, neither mandatory.** Both register a kind through `register_item!`; they
   differ only in what the `entries` callback returns — your own `AbstractDataItem` subtype (type API)
   or the package's `DataItem` (recipe API). The recipe API is a *factory* for a `DataItem`, not a
   competitor. The public contract is the `AbstractDataItem` interface (methods).
3. **The tree is a view, not an identity.** Where an item sits in the hierarchy is *computed by a
   `collect` function*, not baked into the item. The filesystem tree is just the default collection.

---

## What stays

This is a vocabulary + extensibility change, not a rewrite. Unchanged:

- the **project / package split** — project code interprets and presents items; package code owns
  cache, jobs, selection, views, persistence ([workspace-vision.md](workspace-vision.md));
- the **scan pipeline** (`collect → interpret → index`) and per-kind scan profiling;
- the **hierarchy as a fundamental browsing structure** — now realized as *the default view*;
- **stable identity, `parameters`, and `stats`** as the per-item metadata buckets
  ([measurement-parameters-and-stats.md](measurement-parameters-and-stats.md)).

---

## Core types

```julia
abstract type AbstractDataItem end   # the contract every item satisfies
abstract type Collection end         # a semantic, metadata-bearing node in the tree (a device, a dataset)

# The interface IS the contract — a subtype indexes if it answers these. A type implements them
# however it likes: the package's DataItem reads dicts; a project subtype reads its own typed fields.
id(::AbstractDataItem)::String              # stable identity
item_label(::AbstractDataItem)::String           # display title
kind(::AbstractDataItem)::Symbol                 # coarse tag for icons / UI bucketing / register_plot! key
collection(::AbstractDataItem)::Vector{String}   # canonical placement (see "The tree is a view")
parameters(::AbstractDataItem)::Dict{Symbol,Any} # metadata known at interpret time
stats(::AbstractDataItem)::Dict{Symbol,Any}      # values computed later
item_data(::AbstractDataItem)                    # -> Any: the loaded data (DataFrame, image, …)
read_data(::AbstractDataItem)                    # engine-internal loader
# Drawing is NOT part of this interface. Items are drawn by plots registered through
# `register_plot!` (see plotting-api-design.md), whose `draw` reads item.data.
# Optional, with sensible defaults:
process(::AbstractDataItem, data)                # default: passthrough
cacheable(::AbstractDataItem)::Bool              # default: false (cache deferred; DataFrame path unchanged)
```

The package ships `AbstractDataItem`, the interface, the concrete `DataItem`, and the internal
`ItemRecord`. Domain-named concrete types are **project** code, not shipped:

```julia
# package-provided:
struct DataItem <: AbstractDataItem … end       # the normal item: parameters/stats dicts + data
# (+ abstract AbstractDataItem / Collection, the interface, and the internal ItemRecord — see below)

# project-provided via the type API — illustrative, NOT shipped by the package:
struct PundMeasurement <: AbstractDataItem … end   # metadata as typed fields it owns
struct Device          <: Collection        … end  # area_um2, t_HZO_nm, …
```

This generalizes today's package types: metadata-only `MeasurementInfo` becomes the internal
`ItemRecord`; `DeviceInfo`'s path folds into `collection` and its node metadata into the `Collection`
abstraction. Neither a measurement nor a device *type* is package-provided anymore.

`kind` is a coarse tag (icon, UI bucket, plot-registry key), **not** the dispatch key it is today —
the type carries the real meaning. Recipe-API projects still get a `Symbol`-keyed experience through
the package's `DataItem`.

Item **identity is file + kind + params**, never its tree position. `Collection` is the typed,
metadata-bearing entity attached to a *meaningful* node — a device with area/thickness, a dataset
with provenance. Plain intermediate path segments are just strings; you only reach for a `Collection`
when a node carries semantics or metadata, rather than forcing every segment to be a heavy object.

## Two representations: internal record vs. the items you see

The bulk index holds *every* item; storing each one's data there is the unbounded-memory regression
we already fixed. So metadata and data live in two forms, bridged by the engine:

- **`ItemRecord` — completely internal, never exported, never named by project code.** The data-less
  metadata form the hierarchy, scan, and cache store: id, label, `kind`, `collection`,
  `collection_metadata`, `parameters`, `stats`. Today's metadata-only `MeasurementInfo` is an
  `ItemRecord` in all but name. `register_*` callbacks never see it, and it is **not** a field of any
  item — it is a separate, parallel type the engine converts to and from.
- **`AbstractDataItem` instances — what callbacks see.** The loaded, data-bearing values: the
  package's `DataItem`, or a project's own subtype. They carry `item.data`, materialized **only for
  the viewed selection**.

The bridge is engine-owned and uses **the contract, never a shared field**:

```
index : ItemRecord built from any item by calling the contract
        → ItemRecord(id(x), item_label(x), kind(x), collection(x), parameters(x), stats(x))
view  : a DataItem built from a stored record + freshly loaded data
        → DataItem(record, data); item.data is a plain field on the result
```

Consequences:

- **The contract is the interface, never a field.** A subtype indexes if it implements the methods —
  there is nothing like `record::ItemRecord` to declare or forget. A missing method is a clear "you
  didn't implement `kind`," not a struct-shape mismatch.
- **`ItemRecord` is invisible to extenders.** Only code modifying the engine meets it; writing a new
  item kind never does, and no item embeds it.
- **Memory is bounded by construction.** `ItemRecord` has no data slot; only the selected handful
  ever become data-bearing items.
- **No divergence for the recipe path.** A `DataItem` built from a record at view time shares the
  record's `parameters`/`stats` dicts, so engine-computed stats (from `register_collection_stat!`)
  are already present. Project subtypes that go fully typed reproduce their fields at view time via
  the project's loader and own their own metadata.

---

## Two first-party APIs

Both are first-party and supported, and they share **one entry point** — `register_item!(project,
:kind; detect, read, entries, …)`. The difference is purely **what `entries(file, data)` returns**:
the package's `DataItem` (recipe API) or your own `AbstractDataItem` subtype (type API). Both satisfy
the same interface, so nothing downstream can tell them apart. (A subtype cannot be scanned by
existence alone — going from files to items needs `detect`/`read`/`entries` — so the type API still
registers a kind; it does not rely on reflection over `subtypes(AbstractDataItem)`.)

| | Type API | Recipe API |
|---|---|---|
| You write | `struct CVImage <: AbstractDataItem` + interface methods, and `entries` returns `CVImage`s | `entries` returns the package's `DataItem`s |
| Meaning via | the subtype's interface methods | the generic `parameters` / `stats` dicts + optional `process`/`stats`/`label` callbacks |
| Fields | typed (`width`, `height`, …) | the generic `parameters` / `stats` dicts (no declared schema) |
| Used by | the package's own model, power users, plugins | quick projects, one-offs, scripts |
| Produces | the subtype value directly (carries its own data) | a `DataItem` whose data the engine resolves via `read`/`process` |

The package defines the model itself through the **type API** — the interface and the concrete
`DataItem` — and projects or plugins reach for the same type API when they want typed kinds. This is
the intentional version of what `AbstractProject` gestured at. Subtyping is **never mandatory**: the
recipe API is the lightweight path most projects use, and recipe authors put their data in
`parameters`/`stats` dicts (no schema — if you want typed named fields, that is precisely what the
type API is for).

---

## The tree is a view, produced by `collect`

An item's place in the hierarchy is produced by a **`collect` function** `item -> Vector{String}`,
not stored as identity. This single mechanism serves three needs that otherwise look like separate
features:

- **Default collection = the filesystem path** (relative to root). Zero project code: point at a
  folder, get the folder tree.
- **Project override = a simple function.** RuO2 supplies one that parses the filename into
  `chip/layer/device/geometry` — the same parsing it does today, relocated from "set `location`
  inside `measurements`" to "the project's `collect` function."
- **Live UI re-view = swap among built-in `collect_by` axes.** By folder, by `kind`, by date, by any
  `parameter` — instant, because it only re-collects in-memory items. No rescan, no scripting, no tags.

### Canonical vs. view collection

Two roles, kept distinct:

- **Canonical collection** — the project's default (filesystem, or RuO2's device path). Its output is
  stored as each item's `collection` at scan time. This is the **identity / metadata anchor**:
  `metadata.txt` path matching, tag assignments, and annotations all key off `collection`.
- **View collections** — any built-in or project-provided `collect_by` axis the UI applies for
  *presentation*. Computed on demand from item fields; they never touch the stored canonical
  `collection`.

So the hierarchy is "the default view" concretely: the canonical collection is the stable skeleton;
the UI can re-project items along other axes freely without disturbing identity or metadata.

### Tags (later)

Tags are a **cross-cutting axis** — flat, many-to-many, extrinsic ("wakeup-study", "redo"). They are
**not a second tree system**; they feed the *same* view mechanism. The failure mode to avoid is
bespoke tag-tree code duplicating the hierarchy. First use of tags is a **filter** over the canonical
tree (cheapest, most useful); tag-as-`collect_by`-axis comes after. Tags stay **distinct from
`collection`** in the data model — a structural path and a flat label set behave differently, and
merging them into one undifferentiated list is exactly the generality-induced vagueness Principle 1
forbids.

---

## Data and cache

The data type goes from `DataFrame` to `Any`. The load↔draw contract is owned locally by each item
kind. GLMakie renders images natively (`image!`), so an `:image` kind needs only a `read_data`; it is
drawn by the plot registered through `register_plot!` (see
[plotting-api-design.md](plotting-api-design.md)), whose `draw` reads `item.data` — the plot panel is
unchanged.

A view (table inspector, plot, anything) just tries to consume `item.data` and lets dispatch decide —
a 2-D array is table-viewable as much as it is image-viewable. A view never pre-blocks an item by
`kind` or payload type.

Because a loaded item carries its own `data`, plot/inspect callbacks take **items**, not a parallel
data array — `draw(workspace, items, figure)` reading `item.data`, instead of today's
`draw(workspace, measurements, processed, figure)` zipped by index. (`payload` as a name is gone; the
public surface is `item.data`, and the engine-internal loader is `read_data`.)

### Caching stays a trait — deferred, and the transition must stay safe

Caching is a **trait**, not a wrapper type (a `Cacheable{T}` wrapper would force wrap/unwrap at every
boundary):

```julia
cacheable(::Type) = false                   # arbitrary data: re-read from disk on demand
cacheable(::Type{DataFrame}) = true         # tabular data: HDF5-cached as today
store_data(io, x)                           # implemented for cacheable data types
load_data(io, ::Type{T})
```

The hard requirement is that **the cache stays usable across the transition**: the existing DataFrame
payload cache keeps working exactly as today (it is the `cacheable(DataFrame) = true` case), and any
new, non-tabular data type simply falls back to a disk re-read until it opts in. Image files are
random-access and cheap to re-read, so "not cached yet" is never a correctness problem, only a
performance one. An explicit **opt-out** disables caching regardless of type:
`register_item!(…; cache=false)` or `cacheable(::MyItem) = false`. Wiring `store_data`/`load_data`
for new data types can land whenever it is convenient — it does not block the item/data split.

### Cache invalidation under Revise — the minimal mechanism

The only invalidation requirement: **edit an analysis function, and the cache recomputes.** Achieve
it by adding the **project source file's fingerprint** to the cache key, alongside the data-file
fingerprints already used. Editing and saving a `process`/`stats` function changes that file's hash,
so the next `open_workspace`/rescan misses and recomputes. That is the entire mechanism — no
content-hashing of closures, no per-item version keys, no compatibility handling. (A project
defined directly in the REPL has no file to fingerprint; rescan manually — it's the throwaway case.)

---

## Source layer: `AbstractDataSource` (the low-level foundation)

The item/record split above is necessary but not sufficient: it says what a browsable item *is*, but
not where items come from or how their origin is opened, scanned, invalidated, and closed. The source
layer makes that explicit and promotes it to the first-party low-level API the rest of the engine — and
the callback API — is built on. Three abstract types:

```text
AbstractDataSource      → lifecycle + discovery   (dataset root, DB query, instrument, stream, service)
AbstractDataSourceItem  → one discovered unit      (file, row, run id, channel, image-stack member)
AbstractDataItem        → one logical browser item  (unchanged contract above)
```

The middle type earns its keep: scanning a source and interpreting one discovered unit are separate
concerns, and the unit is the natural grain of scan progress, failure reporting, and invalidation —
which the current code already treats `SourceFile` as, implicitly.
`data_items(project, source, source_item)` maps one unit to zero/one/many data items;
`load_data_item(project, source, source_item_id, id)` reloads one with its payload.

**No item registration in the low level.** A workspace still starts from a user-created project plus a
configured source value (`open_workspace(project, mysource)`), never by walking `subtypes`. The
exported callback API (`define_project` + `register_*`) uses `DirectorySource` for directory
discovery and project methods for recipe interpretation. `open_workspace` always takes the project
first; the second argument is either a data root or an explicit `AbstractDataSource`.

### Resolved decisions

1. **Directory metadata belongs to `DirectorySource`; annotations remain a cache-adjacent stopgap.**
   `metadata.txt` is loaded by the directory source. Tags/notes/layout still need a generalized
   source-owned storage capability later.
2. **The phase-3b "type API via `register_item!`" is superseded** by the source-value path — it is the
   cleaner realization of "the type API does not call `register_item!`". Rework, not extend.
3. **`collection_stats(project, source, collection, items)` is a low-level hook available everywhere**, stored
   on the `HierarchyNode` (node-level `stats`), not folded into member records. `register_collection_stat!`
   is its callback form. Bangless: it is a getter returning a `Dict`; the engine routine that writes the
   node is the mutator (`add_collection_stats!`).
4. **`load_data_item` takes id strings and re-looks-up / re-reads.** Accepted for now (FIXME-marked):
   the alternative of threading a source-item handle through the index is deferred.
5. **`ItemRecord` carries only source-*item* identity; source-*level* identity lives once on
   `SourceScan`** — no per-record duplication of `source_id`/`source_label`.
6. **The low-level types are not exported yet** — reachable as `MeasurementBrowser.name`, staged for a
   dedicated submodule. The exported surface stays the conservative high-level set; `PlotKind` stays
   internal too.
7. **`!` follows Julia's argument-mutation convention**, not "has side effects": `close_source!`
   mutates the source → bang; `open_source`/`source_items`/`data_items`/`load_data_item`/`watch_source`
   (default `nothing`) return values or observe → no bang.

## Anti-scope (what we are deliberately *not* doing)

- No general "view engine" up front — built-in collections + a filter, not a query language.
- No multi-parent *storage*. One canonical `collection` per item; multi-membership is a view/tag
  concern layered on top later.
- No data type hierarchy — data is `Any`; the contract is local to each kind.
- No required fields on item types — the contract is the `AbstractDataItem` interface, not a struct
  shape.
- No mandatory subtyping, no declared schema for the generic `DataItem`.
- No per-item cache version keys or content hashing.

---

## Relationship to existing plans

This doc owns the **core data-model generalization**: the item/collection type system, the two-API
contract, `collection` + the `collect` model, the data/cache trait, and cache-on-Revise. It does
**not** own:

- visualizers, figures, workflows, annotations, GUI/API parity — [workspace-vision.md](workspace-vision.md);
- spatial navigation and annotation storage — [spatial-browser.md](spatial-browser.md);
- the public meaning of `parameters` vs `stats` —
  [measurement-parameters-and-stats.md](measurement-parameters-and-stats.md).

When this lands, the current-state [../data-model.md](../data-model.md) and
[../ARCHITECTURE.md](../ARCHITECTURE.md) become the reference for the new types.

---

## Sequencing

Each step is independently testable and ends at a clean, restart-and-run state.

1. **Rename pass.** `measurement_kind → kind`; scope the word "recipe" to the callback path. Mechanical.
   *(The package rename `MeasurementBrowser → DataBrowser` is the most invasive move — it breaks
   `using MeasurementBrowser` for the external TASE/v2-RuO2 projects — so it is its own deliberate
   step, likely paired with the next release rather than folded in here.)*
2. **Interface + record/item split + `item.data`.** Introduce `AbstractDataItem`/`Collection` and the
   interface as the contract; rename today's metadata-only `MeasurementInfo` to the internal
   `ItemRecord` the hierarchy stores; make the recipe path build the package's `DataItem`; materialize
   data-bearing items (with `item.data`) for the viewed selection via the engine bridge. Switch
   plot/inspect callbacks to `(workspace, items, figure)` reading `item.data` — the parallel data
   array goes away. Fold `DeviceInfo`'s path into `collection` and its metadata into the `Collection`
   representation. No package-provided domain types — projects define their own subtypes or use the
   generic `DataItem`.
3. **Collection model.** Replace stored `DeviceInfo.location` with a derived `collection` produced by
   a `collect` function: default = filesystem, project-overridable. Keep canonical vs. view
   distinction; re-key `metadata.txt`/annotations off `collection`.
4. **Data `Any` + cache trait.** Loosen the data type from `DataFrame` to `Any`; add `cacheable`/
   `store_data`/`load_data` and the `cache=` opt-out; add the project-source fingerprint to the cache
   key. Land whenever convenient — the DataFrame cache keeps working untouched until then.
5. **First non-measurement kind.** Add an `Image` kind end-to-end *as an example project / test
   fixture* (not shipped in the package) — the proof the model holds for non-tabular data.
6. **Live re-collecting + tags.** Built-in `collect_by` axes selectable in the UI; tag filtering. Its
   own design pass, coordinated with the annotation plans.

---

## Deferred (decide when the step arrives)

- Arbitrary custom view axes beyond the built-ins — owner is likely "the project declares them,"
  decided at step 6.
- Whether `Collection` metadata (device area/thickness) keeps the current lookup-time path-prefix
  inheritance or becomes explicit `Collection` objects — settle during step 3.
</content>
