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

1. **Generality at the base, clarity at the leaves.** Abstract umbrella types (`DataItem`,
   `DataGroup`) are allowed to be bland — *no one reads them in normal use*. Their meaning is a small
   documented interface, not the noun. The clarity lives in the concrete types **projects** define
   for their own domain (an RuO2 measurement, a device) — or in well-chosen `kind`s and `parameters`
   on a `GenericItem`. The package ships the abstract base, the interface, and `GenericItem`; it
   ships **no** domain types.
2. **Two first-party APIs, neither mandatory.** A type API (subtype + dispatch) and a recipe API
   (register callbacks). The recipe API is a *factory* for the type API, not a competitor. The public
   contract is the `DataItem` interface (methods), never a required field.
3. **The tree is a view, not an identity.** Where an item sits in the hierarchy is *computed by a
   grouping function*, not baked into the item. The filesystem tree is just the default grouping.

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
abstract type DataItem end     # one browsable, viewable thing backed by a source
abstract type DataGroup end    # a semantic, metadata-bearing node in the tree (a device, a dataset)

# The interface IS the contract — a DataItem subtype indexes if it answers these, with no required
# fields. `data` is the one stored field the *loaded* item adds; see "Two representations".
item_id(::DataItem)::String                 # stable identity
item_label(::DataItem)::String              # display title
kind(::DataItem)::Symbol                    # coarse tag for icons / UI bucketing
group_path(::DataItem)::Vector{String}      # canonical placement (see "The tree is a grouping")
parameters(::DataItem)::Dict{Symbol,Any}    # metadata known at interpret time
stats(::DataItem)::Dict{Symbol,Any}         # values computed later
read_data(::DataItem)                       # -> Any: load this item's data (DataFrame, image, …)
render!(ax, ::DataItem)                     # draw it — reads item.data
# Optional, with sensible defaults:
process(::DataItem, data)                   # default: passthrough
cacheable(::DataItem)::Bool                 # default: derived from the data-type trait
```

The package ships the abstract base plus **one** concrete item — `GenericItem`, the recipe-API
backing — and the internal `ItemRecord`. Domain-named concrete types are **project** code, not
shipped:

```julia
# package-provided:
struct GenericItem <: DataItem … end          # recipe-API backing; holds an ItemRecord + data
# (+ abstract DataItem / DataGroup, the interface, and the internal ItemRecord — see below)

# project-provided via the type API — illustrative, NOT shipped by the package:
struct PundMeasurement <: DataItem  … end
struct Device          <: DataGroup … end       # area_um2, t_HZO_nm, …
```

This generalizes today's package types: metadata-only `MeasurementInfo` becomes the internal
`ItemRecord`; `DeviceInfo`'s path folds into `group_path` and its node metadata into the `DataGroup`
abstraction. Neither a measurement nor a device *type* is package-provided anymore.

`kind` is a coarse tag (icon, UI bucket), **not** the dispatch key it is today — the type carries
the real meaning. Recipe-API projects still get a `Symbol`-keyed experience through `GenericItem`.

Item **identity is file + kind + params**, never its tree position. `DataGroup` is the typed,
metadata-bearing entity attached to a *meaningful* node — a device with area/thickness, a dataset
with provenance. Plain intermediate path segments are just strings; you only reach for a `DataGroup`
when a node carries semantics or metadata, rather than forcing every segment to be a heavy object.

## Two representations: index record vs. loaded item

Data must never sit in the bulk index — `MeasurementHierarchy` holds *every* item, and storing each
one's data there is the unbounded-memory regression we already fixed. So there are two
representations, bridged by the engine:

- **`ItemRecord` — internal, never exported.** What the hierarchy, scan, and cache store: generic
  metadata only (id, label, `kind`, timestamp, `group_path`, `parameters`, `stats`, tags) — **no
  data** — plus engine-only bookkeeping. Today's metadata-only `MeasurementInfo` is already an
  `ItemRecord` in all but name.
- **`DataItem` — public, with a real `item.data`.** The loaded, renderable form, materialized **only
  for the viewed selection**.

The bridge is engine-owned and **uses the interface, not any field**:

```
index : ItemRecord built from a DataItem by calling its interface methods
        → ItemRecord(item_id(x), kind(x), group_path(x), parameters(x), stats(x), …)
view  : DataItem materialized from a record + its loaded data, via the project's loader
        → item.data is a plain field on the result
```

Consequences:

- **The contract is the interface, never a field.** A `DataItem` subtype indexes if it implements the
  methods — there is nothing like `record::ItemRecord` to declare or forget. A missing method is a
  clear "you didn't implement `kind`," not a struct-shape mismatch.
- **`ItemRecord` is invisible to extenders.** Only code modifying the engine meets it; writing a new
  item kind never does.
- **Memory is bounded by construction.** The index type has no data slot; only the selected handful
  ever become data-bearing items.
- **`GenericItem` avoids duplication privately.** The one concrete item the package owns holds an
  `ItemRecord` internally (reusing the metadata field list) plus `data` — an implementation choice,
  not a public requirement. Project-defined item kinds stay free-form.

---

## Two first-party APIs

Both are first-party and supported. They are not competing models — the recipe API produces a
`GenericItem` that satisfies the same interface the type API implements directly, so nothing
downstream can tell them apart.

| | Type API | Recipe API |
|---|---|---|
| You write | `struct CVImage <: DataItem` + interface methods | `register_item!(project, :image; read=…, render=…)` |
| Meaning via | dispatch — `render!(ax, ::CVImage, …)` | stored callbacks |
| Fields | typed (`width`, `height`, …) | the generic `parameters` / `stats` dicts (no declared schema) |
| Used by | the package's own model, power users, plugins | quick projects, one-offs, scripts |
| Produces | the subtype value directly | a `GenericItem <: DataItem` forwarding to callbacks |

The package defines the model itself through the **type API** — the interface and `GenericItem` — and
projects or plugins reach for the same type API when they want typed, dispatch-driven kinds. This is
the intentional version of what `AbstractProject` gestured at. Subtyping is **never mandatory**: the
recipe API is the lightweight path most projects use, and recipe authors put their data in
`parameters`/`stats` dicts (no schema — if you want typed named fields, that is precisely what the
type API is for).

---

## The tree is a grouping

An item's place in the hierarchy is produced by a **grouping function** `item -> Vector{String}`,
not stored as identity. This single mechanism serves three needs that otherwise look like separate
features:

- **Default grouping = the filesystem path** (relative to root). Zero project code: point at a
  folder, get the folder tree.
- **Project override = a simple function.** RuO2 supplies one that parses the filename into
  `chip/layer/device/geometry` — the same parsing it does today, relocated from "set `location`
  inside `measurements`" to "the project's grouping function."
- **Live UI re-view = swap among built-in groupings.** By folder, by `kind`, by date, by any
  `parameter` — instant, because it only regroups in-memory items. No rescan, no scripting, no tags.

### Canonical vs. view grouping

Two roles, kept distinct:

- **Canonical grouping** — the project's default (filesystem, or RuO2's device path). Its output is
  stored as each item's `group_path` at scan time. This is the **identity / metadata anchor**:
  `device_info.txt` path matching, tag assignments, and annotations all key off `group_path`.
- **View groupings** — any built-in or project-provided grouping the UI applies for *presentation*.
  Computed on demand from item fields; they never touch the stored canonical `group_path`.

So the hierarchy is "the default view" concretely: the canonical grouping is the stable skeleton;
the UI can re-project items along other axes freely without disturbing identity or metadata.

### Tags (later)

Tags are a **cross-cutting axis** — flat, many-to-many, extrinsic ("wakeup-study", "redo"). They are
**not a second tree system**; they feed the *same* view mechanism. The failure mode to avoid is
bespoke tag-tree code duplicating the hierarchy. First use of tags is a **filter** over the canonical
tree (cheapest, most useful); tag-as-grouping-axis comes after. Tags stay **distinct from
`group_path`** in the data model — a structural path and a flat label set behave differently, and
merging them into one undifferentiated list is exactly the generality-induced vagueness Principle 1
forbids.

---

## Data and cache

The data type goes from `DataFrame` to `Any`. The load↔render contract is owned locally by each item
kind. GLMakie renders images natively (`image!`), so an `:image` kind needs only a `read_data` and a
`render!` that reads `item.data` — the plot panel is unchanged. `TableInspector` guards on
`item.data isa DataFrame` and does not open for non-tabular items.

Because a loaded `DataItem` carries its own `data`, plot/inspect callbacks take **items**, not a
parallel data array — `draw(workspace, items, figure)` reading `item.data`, instead of today's
`draw(workspace, measurements, processed, figure)` zipped by index. (`payload` as a name is gone; the
public surface is `item.data`, and the engine-internal loader is `read_data`.)

Caching is a **trait**, not a wrapper type (a `Cacheable{T}` wrapper would force wrap/unwrap at every
boundary):

```julia
cacheable(::Type) = false                   # arbitrary data: re-read from disk on demand
cacheable(::Type{DataFrame}) = true         # tabular data: HDF5-cached as today
store_data(io, x)                           # implemented for cacheable data types
load_data(io, ::Type{T})
```

- Data whose type implements the trait is cached; anything else falls back to disk read gracefully
  (image files are random-access and cheap to re-read).
- An explicit **opt-out** disables caching regardless of type: `register_item!(…; cache=false)` or
  `cacheable(::MyItem) = false`.

### Cache invalidation under Revise — the minimal mechanism

The only invalidation requirement: **edit an analysis function, and the cache recomputes.** Achieve
it by adding the **project source file's fingerprint** to the cache key, alongside the data-file
fingerprints already used. Editing and saving a `process`/`stats` function changes that file's hash,
so the next `open_workspace`/rescan misses and recomputes. That is the entire mechanism — no
content-hashing of closures, no per-item version keys, no compatibility handling. (A project
defined directly in the REPL has no file to fingerprint; rescan manually — it's the throwaway case.)

---

## Anti-scope (what we are deliberately *not* doing)

- No general "view engine" up front — built-in groupings + a filter, not a query language.
- No multi-parent *storage*. One canonical `group_path` per item; multi-membership is a view/tag
  concern layered on top later.
- No data type hierarchy — data is `Any`; the contract is local to each kind.
- No required fields on item types — the contract is the `DataItem` interface, not a struct shape.
- No mandatory subtyping, no declared schema for `GenericItem`.
- No per-item cache version keys or content hashing.

---

## Relationship to existing plans

This doc owns the **core data-model generalization**: the item/group type system, the two-API
contract, `group_path` + the grouping model, the data/cache trait, and cache-on-Revise. It does
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
2. **Interface + record/item split + `item.data`.** Introduce `DataItem`/`DataGroup` and the
   interface as the contract; rename today's metadata-only `MeasurementInfo` to the internal
   `ItemRecord` the hierarchy stores; make the recipe path build a `GenericItem`; materialize
   data-bearing `DataItem`s (with `item.data`) for the viewed selection via the engine bridge. Switch
   plot/inspect callbacks to `(workspace, items, figure)` reading `item.data` — the parallel data
   array goes away. Fold `DeviceInfo`'s path into `group_path` and its metadata into the `DataGroup`
   representation. No package-provided domain types — projects define their own subtypes or use
   `GenericItem`.
3. **Grouping model.** Replace stored `DeviceInfo.location` with a derived `group_path` produced by a
   grouping function: default = filesystem, project-overridable. Keep canonical vs. view distinction;
   re-key `device_info.txt`/annotations off `group_path`.
4. **Data `Any` + cache trait.** Loosen the data type from `DataFrame` to `Any`; add `cacheable`/
   `store_data`/`load_data` and the `cache=` opt-out; add the project-source fingerprint to the cache
   key.
5. **First non-measurement kind.** Add an `Image` kind end-to-end *as an example project / test
   fixture* (not shipped in the package) — the proof the model holds for non-tabular data.
6. **Live re-grouping + tags.** Built-in view groupings selectable in the UI; tag filtering. Its own
   design pass, coordinated with the annotation plans.

---

## Deferred (decide when the step arrives)

- Arbitrary custom view axes beyond the built-ins — owner is likely "the project declares them,"
  decided at step 6.
- Whether `DataGroup` metadata (device area/thickness) keeps the current lookup-time path-prefix
  inheritance or becomes explicit group objects — settle during step 3.
