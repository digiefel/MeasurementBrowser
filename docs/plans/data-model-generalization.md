# DataBrowser: General Data Model

## Goal

Grow MeasurementBrowser into **DataBrowser**: an open, scriptable, data-focused IDE rather than a
measurement-only viewer. The same app should browse measurements, images, fitted models, and
arbitrary datasets; group them by meaning (a device, a wafer, an ML dataset) rather than by folder;
and present each through the visualizers and figures described in
[workspace-vision.md](workspace-vision.md).

The data model we ship today (`MeasurementInfo`, `DeviceInfo`, the device hierarchy) is the right
shape but the wrong vocabulary: almost every field is already domain-neutral, and the parts that are
not — the `DataFrame` payload and the directory-shaped tree — are the two things that block other
kinds of data. This doc defines how the core types generalize without losing the concreteness that
makes the current model legible.

## What stays

This is a vocabulary and extensibility change, not a rewrite. These do not change:

- the **project / package split**: project code identifies and interprets items; package code owns
  cache, jobs, selection, views, and persistence (see [workspace-vision.md](workspace-vision.md));
- the **scan pipeline** (`collect → interpret → index`) and per-kind scan profiling;
- the **hierarchy** as a fundamental browsing structure — items live in a navigable tree;
- **stable identity, parameters, and stats** as the per-item metadata buckets
  ([measurement-parameters-and-stats.md](measurement-parameters-and-stats.md)).

## Principle: generality at the base, clarity at the leaves

The risk in this work is trading a concrete, communicative `DeviceInfo` for a vague `DataGroup` and
calling it progress. We avoid that by being explicit about where generality is allowed to live:

- **Abstract types are general by necessity.** An umbrella over "a measurement and an image and an
  ML dataset" cannot have a domain-specific name. Its *meaning* comes from a small, documented
  interface — not from the noun.
- **Concrete first-party types carry the domain weight.** `DeviceInfo`'s concept does not disappear;
  it becomes a concrete `Device <: DataGroup` with the same area/thickness fields it has today.
  `MeasurementInfo` becomes a concrete `Measurement <: DataItem`. Everyday code, first-party code,
  and the RuO2 project keep naming `Measurement` and `Device` — as concrete as before, minus the
  `Info` suffix noise.
- **The abstract umbrella is only seen when writing generic tools.** A visualizer that draws "any
  tabular item" dispatches on the interface; a project that just wants to browse PUND data never
  names `DataItem` at all.

So the abstract names (`DataItem`, `DataGroup`) are deliberately bland, and that is fine, because
**no one reads them in normal use**. If a less generic abstract name reads better we can revisit, but
the clarity budget is spent on the concrete subtypes, not the base.

## Two first-party APIs

Both extension styles are first-party and supported. They are not competing models — the recipe API
is a *factory* for the type API.

| | Type API | Recipe API |
|---|---|---|
| You write | `struct CVImage <: DataItem` + interface methods | `register_item!(project, :image; read=…, render=…)` |
| Dispatch carries meaning | yes — `render!(ax, ::CVImage, …)` | no — behavior is in stored callbacks |
| Typed fields | yes (`width`, `height`, …) | no — generic `parameters`/`stats` buckets |
| Used by | the package itself, power users, plugins | quick projects, one-offs, scripts |
| Produces | the subtype value directly | a `GenericItem <: DataItem` forwarding to callbacks |

The package's own first-party item kinds (`Measurement`, `Image`, …) are defined through the **type
API** — this is what `AbstractProject` gestured at, made intentional. Subtyping is **not mandatory**:
the recipe API stays the lightweight path, and a `GenericItem` satisfies the same interface so the
rest of the app cannot tell the difference.

## Core types (sketch)

```julia
abstract type DataItem end     # one browsable, viewable thing backed by a source
abstract type DataGroup end    # a semantic group an item belongs to

# The interface = the meaning. Required of every DataItem subtype:
item_id(::DataItem)::String                 # stable identity
item_label(::DataItem)::String              # display title
item_kind(::DataItem)::Symbol               # coarse tag for icons / UI bucketing
read_payload(::DataItem)                    # -> Any: DataFrame, image, model, …
render!(ax, ::DataItem, payload)            # how it draws into a figure/panel
groups(::DataItem)::Vector{DataGroup}       # semantic membership (see below)
# optional, with sensible defaults:
process(::DataItem, payload)                # default: passthrough
stats(::DataItem, payload)::Dict{Symbol,Any}
cacheable(::DataItem)::Bool                 # default: from payload type trait

# First-party concrete types keep the old clarity:
struct Measurement <: DataItem … end        # was MeasurementInfo
struct Device      <: DataGroup … end        # was DeviceInfo (area_um2, t_HZO_nm, …)
struct GenericItem <: DataItem … end         # recipe-API backing type
```

`item_kind` is a coarse tag (icon, UI grouping), no longer the dispatch key it is today — the type
carries the real meaning. Projects that use the recipe API still get a `Symbol`-keyed experience.

## Payload and cache

The payload type goes from `DataFrame` to `Any`. The read↔render contract is owned locally by each
item kind, exactly as recipes already do for tables. GLMakie renders images natively (`image!`),
so an `:image` kind needs only a reader and a `render!` — the plot panel is unchanged.

Caching becomes a **trait**, not a wrapper type (a `Cacheable{T}` wrapper would force
wrap/unwrap at every boundary):

```julia
cacheable(::Type) = false                   # arbitrary payloads: re-read from disk on demand
cacheable(::Type{DataFrame}) = true         # tabular payloads: HDF5-cached as today
store_payload(io, x)                        # implemented for cacheable payload types
load_payload(io, ::Type{T})
```

- A payload whose type implements the cache interface is cached; anything else **falls back to disk
  read gracefully** (image files are already random-access and cheap to re-read).
- An explicit **opt-out keyword** lets a user disable caching regardless of type:
  `register_item!(…; cache=false)` / `cacheable(::MyItem) = false`. Useful for volatile sources or
  when re-read is cheaper than cache churn.

The `TableInspector` stays DataFrame-specific: it guards on `payload isa DataFrame` and simply does
not open for non-tabular items.

## Semantic groups, tags, and views

Groups are **not** directory paths. An image can belong to a `Device` exactly as a measurement does,
or to an ML `Dataset`, regardless of where its file sits. This is already half-true today — the
device location comes from filename parsing, not `dirname` — and we make it explicit: an item
*declares* its group membership via `groups(::DataItem)`.

This opens three related capabilities, in increasing scope:

1. **Semantic, multi-membership grouping (tags).** An item may belong to several groups at once
   (a device *and* a dataset *and* user-applied tags). Tags were already planned as annotations
   ([workspace-vision.md](workspace-vision.md), [spatial-browser.md](spatial-browser.md)); here they
   become a first-class grouping axis, not just labels.
2. **Multiple browse views.** A *view* projects items into a browsable tree along a chosen axis —
   by device, by date, by tag, by dataset — with optional filters. The same items, browsed many
   ways.
3. **Live and scripted views.** Views are definable both from the Julia API and built live in the
   UI, consistent with GUI/API parity.

To stay grounded, the migration keeps **one primary tree** (the current device hierarchy, now built
from a declared grouping axis instead of an assumed directory path) and adds tags and alternate
views incrementally. True multi-parent navigation is a view-layer feature layered on top, not a
change to how a single item is stored.

## Open design questions

- The abstract base names (`DataItem`/`DataGroup`) — keep deliberately bland, or find a base noun
  that reads better without overpromising?
- Does `groups(::DataItem)` return a list (multi-membership from day one) or a single primary group
  with tags layered separately? Leaning list, primary-tree-first in the UI.
- How does a `GenericItem` expose typed-ish fields the recipe author cares about — purely through
  `parameters`/`stats`, or a small declared schema?
- Cache identity for non-tabular payloads: fingerprint-by-file is enough for images, but processed
  payloads (fits, derived arrays) need a content/version key — shared with the existing project
  cache fingerprinting or separate?
- Where the grouping axis is declared: on the project, on the item kind, or on the view.

## Relationship to existing plans

This doc owns the **core data-model generalization**: the item/group type system, the two-API
contract, the payload/cache trait, and the move to semantic groups. It does **not** own:

- visualizers, figures, workflows, annotations, and GUI/API parity — owned by
  [workspace-vision.md](workspace-vision.md);
- spatial navigation and device/measurement annotation storage —
  [spatial-browser.md](spatial-browser.md);
- the public meaning of `parameters` vs `stats` —
  [measurement-parameters-and-stats.md](measurement-parameters-and-stats.md).

When this generalization lands, the current-state [../data-model.md](../data-model.md) and
[../ARCHITECTURE.md](../ARCHITECTURE.md) become the reference for the new types.

## Suggested sequencing

Each step is independently testable and ends at a clean, restart-and-run state:

1. **Rename pass.** `MeasurementBrowser → DataBrowser`, `measurement_kind → kind`, scope the word
   "recipe" to the callback path. Mechanical; no behavior change.
2. **Abstract interface + `GenericItem`.** Introduce `DataItem`/`DataGroup` and the interface; make
   the existing recipe path build a `GenericItem`; redefine `Measurement`/`Device` as concrete
   subtypes. No behavior change — the abstraction wraps what exists.
3. **Payload `Any` + cache trait.** Loosen the payload type, add `cacheable`/`store_payload`/
   `load_payload` and the `cache=` opt-out, fall back to disk read for non-cacheable payloads.
4. **First non-measurement kind.** Add an `Image` item kind end-to-end as the proof the model holds.
5. **Tags and views.** Multi-membership grouping and alternate browse views — its own design pass,
   coordinated with the annotation plans.
