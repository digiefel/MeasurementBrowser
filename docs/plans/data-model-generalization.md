# DataBrowser: General Data Model

## Goal

Grow MeasurementBrowser into **DataBrowser**: an open, scriptable, data-focused IDE rather than a
measurement-only viewer. The same app should browse measurements, images, fitted models, and
arbitrary datasets; organize them by meaning rather than by folder; and present each through the
visualizers and figures in [workspace-vision.md](workspace-vision.md).

The data model we ship today is the right shape with the wrong vocabulary and two hard couplings.
Almost every field of `MeasurementInfo` is already domain-neutral. The parts that block other kinds
of data are exactly two:

1. the payload is pinned to `DataFrame`, and
2. an item's place in the tree is a stored, scan-time field (`DeviceInfo.location`) rather than a
   derived view.

This doc defines how the core types generalize while keeping the concreteness that makes the current
model legible.

---

## Principles

These three decide every ambiguous call below.

1. **Generality at the base, clarity at the leaves.** Abstract umbrella types (`DataItem`,
   `DataGroup`) are allowed to be bland — *no one reads them in normal use*. Their meaning is a small
   documented interface, not the noun. The clarity budget is spent on concrete first-party types
   (`Measurement`, `Device`, `Image`), which stay as concrete as `MeasurementInfo`/`DeviceInfo` are
   today.
2. **Two first-party APIs, neither mandatory.** A type API (subtype + dispatch) and a recipe API
   (register callbacks). The recipe API is a *factory* for the type API, not a competitor.
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

# The interface IS the meaning. Required of every DataItem subtype:
item_id(::DataItem)::String                 # stable identity
item_label(::DataItem)::String              # display title
kind(::DataItem)::Symbol                    # coarse tag for icons / UI bucketing
read_payload(::DataItem)                    # -> Any: DataFrame, image, model, …
render!(ax, ::DataItem, payload)            # how it draws into a figure/panel
# Optional, with sensible defaults:
process(::DataItem, payload)                # default: passthrough
stats(::DataItem, payload)::Dict{Symbol,Any}
cacheable(::DataItem)::Bool                 # default: derived from the payload-type trait
```

First-party concrete types carry the old clarity:

```julia
struct Measurement <: DataItem  …  end       # was MeasurementInfo
struct Device      <: DataGroup …  end       # was DeviceInfo (area_um2, t_HZO_nm, …)
struct GenericItem <: DataItem  …  end        # the recipe-API backing type
# later: struct Image <: DataItem … end
```

`kind` is a coarse tag (icon, UI bucket), **not** the dispatch key it is today — the type carries
the real meaning. Recipe-API projects still get a `Symbol`-keyed experience through `GenericItem`.

Each `DataItem` carries the buckets it has today — `parameters` (known at interpret time), `stats`
(computed later) — plus a stored **`group_path::Vector{String}`** (its canonical placement; see
"The tree is a grouping") and a `tags` set (future; see "Tags"). Its identity is **file + kind +
params**, never its tree position.

`DataGroup` is the typed, metadata-bearing entity attached to a *meaningful* node — a device with
area/thickness, a dataset with provenance. Plain intermediate path segments are just strings; you
only reach for a `DataGroup` when a node carries semantics or metadata. (This keeps the path simple
and `DataGroup` meaningful, rather than forcing every segment to be a heavy object.)

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
| Used by | the package itself, power users, plugins | quick projects, one-offs, scripts |
| Produces | the subtype value directly | a `GenericItem <: DataItem` forwarding to callbacks |

The package's own first-party kinds use the **type API** — this is the intentional, well-designed
version of what `AbstractProject` gestured at. Subtyping is **never mandatory**: the recipe API stays
the lightweight path, and recipe authors put their data in `parameters`/`stats` dicts (no schema —
if you want typed named fields, that is precisely what the type API is for).

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

## Payload and cache

The payload type goes from `DataFrame` to `Any`. The read↔render contract is owned locally by each
item kind, exactly as recipes already do for tables. GLMakie renders images natively (`image!`), so
an `:image` kind needs only a reader and a `render!` — the plot panel is unchanged. `TableInspector`
guards on `payload isa DataFrame` and does not open for non-tabular items.

Caching is a **trait**, not a wrapper type (a `Cacheable{T}` wrapper would force wrap/unwrap at every
boundary):

```julia
cacheable(::Type) = false                   # arbitrary payloads: re-read from disk on demand
cacheable(::Type{DataFrame}) = true         # tabular payloads: HDF5-cached as today
store_payload(io, x)                        # implemented for cacheable payload types
load_payload(io, ::Type{T})
```

- A payload whose type implements the trait is cached; anything else falls back to disk read
  gracefully (image files are random-access and cheap to re-read).
- An explicit **opt-out** disables caching regardless of type: `register_item!(…; cache=false)` or
  `cacheable(::MyItem) = false`.

### Cache invalidation under Revise — the minimal mechanism

The only invalidation requirement: **edit an analysis function, and the cache recomputes.** Achieve
it by adding the **project source file's fingerprint** to the cache key, alongside the data-file
fingerprints already used. Editing and saving a `process`/`stats` function changes that file's hash,
so the next `open_workspace`/rescan misses and recomputes. That is the entire mechanism — no
content-hashing of closures, no per-payload version keys, no compatibility handling. (A project
defined directly in the REPL has no file to fingerprint; rescan manually — it's the throwaway case.)

---

## Anti-scope (what we are deliberately *not* doing)

- No general "view engine" up front — built-in groupings + a filter, not a query language.
- No multi-parent *storage*. One canonical `group_path` per item; multi-membership is a view/tag
  concern layered on top later.
- No payload type hierarchy — payload is `Any`; the contract is local to each kind.
- No mandatory subtyping, no declared schema for `GenericItem`.
- No per-payload cache version keys or content hashing.

---

## Relationship to existing plans

This doc owns the **core data-model generalization**: the item/group type system, the two-API
contract, `group_path` + the grouping model, the payload/cache trait, and cache-on-Revise. It does
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
2. **Abstract interface + `GenericItem`.** Introduce `DataItem`/`DataGroup` and the interface; make
   the recipe path build a `GenericItem`; redefine `Measurement`/`Device` as concrete subtypes. No
   behavior change — the abstraction wraps what exists.
3. **Grouping model.** Replace stored `DeviceInfo.location` with a derived `group_path` produced by a
   grouping function: default = filesystem, project-overridable. Keep canonical vs. view distinction;
   re-key `device_info.txt`/annotations off `group_path`.
4. **Payload `Any` + cache trait.** Loosen the payload type; add `cacheable`/`store_payload`/
   `load_payload` and the `cache=` opt-out; add the project-source fingerprint to the cache key.
5. **First non-measurement kind.** Add an `Image` kind end-to-end — the proof the model holds.
6. **Live re-grouping + tags.** Built-in view groupings selectable in the UI; tag filtering. Its own
   design pass, coordinated with the annotation plans.

---

## Deferred (decide when the step arrives)

- Arbitrary custom view axes beyond the built-ins — owner is likely "the project declares them,"
  decided at step 6.
- Whether `DataGroup` metadata (device area/thickness) keeps the current lookup-time path-prefix
  inheritance or becomes explicit group objects — settle during step 3.
