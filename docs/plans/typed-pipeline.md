# Typed Pipeline Foundation

## Purpose

The type API is the fundamental DataBrowser project API. The registration API is a convenience
dialect written over ordinary `(data, metadata)` values. It must not have more pipeline stages,
caching capabilities, or engine integration than typed projects.

The current implementation violates that layering. Typed projects provide one monolithic
`data_items` method; registration separately exposes `detect`, `read`, `entries`, and item and
collection `process`/`analyze`. The engine branches on `item isa RegisteredDataItem`, the cache
gates disk persistence on the same test, and `ItemRecipe` stores abstract `Function` fields. This
leaves typed projects second-class and hides the read/entries cache boundary.

Replace both paths with one staged, multiple-dispatch engine whose only vocabulary is the typed
stage contract, and move the registration dialect into its own package, **`DataBrowserRecipes`**,
implemented entirely as a client of that contract. Remove `data_items` cleanly; there are no
compatibility shims in the pre-alpha API.

## The stage contract

Source discovery remains owned by the source (`source_items(source)`). The pipeline stages are
generic functions declared in `DataBrowserAPI`. The engine calls only the project-aware forms;
their defaults forward to the context-free forms projects usually implement:

```julia
read(project, source, item::AbstractDataSourceItem)         # -> the project's loaded value
entries(project, item::AbstractDataSourceItem, loaded)      # -> Vector{<:AbstractDataItem}
process(project, item::AbstractDataItem)                    # -> AbstractDataItem
analyze(project, item::AbstractDataItem)                    # -> Dict
process(project, collection::AbstractCollection, items)     # -> Vector{<:AbstractDataItem}
analyze(project, collection::AbstractCollection, items)     # -> Dict

read(source, item)                  # context-free forms with identity/empty defaults;
entries(item, loaded)               # the default entries wraps one loaded value into one item
process(item);  analyze(item)
process(collection, items);  analyze(collection, items)
```

`read` performs the expensive source operation once. `entries` expands the result into zero, one,
or many concrete data items without rereading the source. A typed project that owns its source
dispatches the context-free forms on its own types; a typed project reusing a shared source
(`DirectorySource`) dispatches the project-aware forms on its project type.

**The source appears exactly once, at `read`.** Every later stage is a pure function of values:
the source item (an address — fingerprintable, recordable), the loaded payload, the items. That
rule is the cache design. A valid cached stage result satisfies every downstream consumer without
the origin existing, which is what makes rerun-entries-without-reread, warm reopen, eviction
recovery, and cross-machine rehydration sound by construction rather than by discipline. A `read`
that returns a live handle (an HDF5 group, a DB cursor) is legitimate but uncacheable; that
project's durable boundary sits at `process`. `entries` still receives the source item because
identity is not payload: the loaded value stays purely the expensive data (cacheable, without
duplicated identity), while default sibling ids, labels, and collection placement derive from the
source item — which is exactly what the engine can supply without touching the origin.

There is no routing stage. Registration `detect` is not a universal typed stage: a typed source
already controls which source-item types it discovers, and the recipes adapter runs `detect`
inside its own `read` method (read runs for every changed source item anyway, and detect costs
filename-predicate time). Cheap kind classification for status surfaces stays a description query
(`detect_kind`), not a pipeline stage.

`read`'s return value and `entries`' `loaded` argument deliberately have no abstract supertype:
the loaded value is a private handoff between two stages of the same project, and the type
discipline is enforced by the user's concrete signature on the receiving end.

## Project contract

`DataBrowserAPI` gains `AbstractProject` and keeps the project declarations (`project_name`,
`project_description`, and eventually `project_fingerprint` for cache identity). The concrete
recipe-holding `Project` struct is the registration dialect's implementation of that contract and
moves to `DataBrowserRecipes`; per-scan profiling state moves out of the user-facing value with
it. Everything in the vision that depends on projects — persistence, `dbproject.toml`, provenance
— attaches to the contract, not to the recipe bag.

## DataBrowserRecipes

The registration dialect becomes its own package: `define_project`, `register_item!`, and
`register_collection_analysis!`, plus premade recipes built on the same public surface
(`register_csv!` first). It depends on `DataBrowserAPI` and `DataBrowserSources` (its `detect`
callbacks and adapter methods dispatch on `SourceFile`); the `DataBrowser` umbrella re-exports it
so `using DataBrowser` keeps working unchanged. It is the first first-party package
built purely on the type API — a living conformance test of the extension surface, the same role
`DataBrowserPlots` plays for the GUI extension surface.

The registration identity lives in a value type parameter:

```julia
struct RegisteredReadResult{K,D}
    data::D
    metadata::MetadataDict
end

struct RegisteredDataItem{K,D} <: AbstractDataItem
    id::String
    label::String
    data::D
    metadata::MetadataDict
    collection::Vector{String}
end
```

For `register_item!(project, :pund; ...)`, `K` is `:pund`; `RegisteredDataItem{:pund,DataFrame}`
is a distinct concrete type, and the registration name is not an untyped runtime field. There is
no `RegisteredSourceItem`: detection happens inside the adapter's `read`, so the `{K}` tag first
appears on the read result:

```julia
function DataBrowser.read(p::RecipeProject, source, f::SourceFile)
    recipe = detect_recipe(p, f)      # first matching detect wins, registration order
    recipe === nothing && return NoMatch()
    invoke_read(recipe, f)            # function barrier -> RegisteredReadResult{K,D}
end

DataBrowser.entries(p::RecipeProject, f::SourceFile, loaded::RegisteredReadResult{K}) where K
DataBrowser.process(p::RecipeProject, item::RegisteredDataItem{K}) where K
DataBrowser.analyze(p::RecipeProject, item::RegisteredDataItem{K}) where K
```

Each adapter method obtains the current recipe for `K`, invokes the adjacent user callback over
ordinary data and metadata, and returns the next package-owned carrier. `entries` on `NoMatch`
yields zero items. Re-registering `:pund` replaces the recipe value in the project; it does not
redefine methods, so there are no world-age or method-accumulation problems. Registration
collection callbacks use the same approach over aligned data/metadata vectors, while the engine
owns the concrete `CollectionRecord` and work key.

## Type-stable callback storage

Replace mutable recipes containing `Function` fields with immutable parametric recipes whose field
types are the concrete callback types:

```julia
struct ItemRecipe{Detect,Read,Entries,Process,Analyze,Label,Collection,Id}
    detect::Detect
    read::Read
    entries::Entries
    process::Process
    analyze::Analyze
    label::Label
    collection::Collection
    id::Id
end
```

Project storage is necessarily heterogeneous because different registrations have different
callback types. Dynamic selection happens once when looking up the recipe for a source item or a
`K`; a function barrier then receives the concrete recipe and specializes the complete stage call.
Do not spread abstract `Function` dispatch through per-item loops.

## Caching and rehydration

The cache stores package-owned records plus payloads in supported shapes, never arbitrary user
structs. Three separate questions decide what happens to a stage result. Today they are tangled
into one boolean, `cacheable(::AbstractDataItem)`, which usefully answers none of them.

1. **Can storage hold this shape?** The engine's own question, answered without the project by
   `_storable_table`. A project whose payload is not tabular opts in by dispatching
   `cacheable_data` on the payload type. No item-level predicate is involved.
2. **Can this value become the user's type again?** The type's question, answered by whether it
   implements `reconstruct` (below). Its presence *is* the declaration; there is no boolean.
3. **Should we bother, given the cost?** Nobody asks this today. If a project ever needs to
   suppress persistence for a cheap-to-recompute payload, the dialect can expose it per
   registration. Do not build it speculatively.

`cacheable(::AbstractDataItem)` is therefore **removed**, not generalized into
`cacheable(stage, value)`. It carries no information now: for typed items the
`isa RegisteredDataItem` gate discards it before it is consulted, and for registered items it is
defined as `cacheable_data(item.data)` — that is, `Tables.istable` — which `_storable_table`
re-checks as its own first line. Per-stage granularity would be illusory besides: only item and
collection `process` produce payloads at all, while `analyze` produces metadata, which is always
stored.

A cached value comes back in two forms: **payload delivery** (views, tables, queries — needs no
user type) and **materialized delivery** (running further user dispatch — needs the real type).
Only the second requires the project to say anything:

```julia
reconstruct(::Type{T}, data, metadata::Dict)::T
```

Package-owned, opt-in, and pure, gated by `hasmethod` — a generic function rather than a
constructor, so opting in is unambiguous and collision-free, and the package never defines methods
on a name a project owns. The rebuilt item re-derives `id`, `collection`, and `metadata`; the
engine validates them against the record and fails loudly on mismatch. Without a `reconstruct`
method the engine falls back to rerunning upstream stages — always correct, just slower, and
payload delivery keeps working meanwhile. `RegisteredDataItem` implements `reconstruct`
internally, which dissolves the registration-only cache gate into dispatch any item type can
implement. Rehydration must stay a pure function of cached content; values a type needs to rebuild
itself belong in its metadata, never in live workspace state.

The engine always takes the cheapest valid path: cached downstream payload → `reconstruct` →
rerun upstream stages. Warm reopen therefore delivers concrete user types at deserialization
speed, threaded across the worker pool, rebuilding only the items a query or selection actually
touches.

## Work and cache stages

Replace the monolithic `SOURCE_INTERPRET` work with distinct per-stage results. The work graph
models the real pipeline:

```text
source discovery
    → read
    → entries
    → item process
    → item analyze
    → collection process
    → collection analyze
```

Each stage gets its own profiling, failure, invalidation, memory-retention, and persistence
boundary. A valid downstream cache result prevents rerunning only the stages it actually subsumes,
and per-stage results are the hook the pipeline-inspection control (roadmap 0.14.0) reads from.

## Macro syntax

A macro may later provide a more declarative layout, but it is not needed to create types at
runtime and must not own different semantics. An `@register_item` block could expand to the same
`register_item!` recipe construction while preserving the current callback form for ordinary code,
Revise, generated projects, and programmatic registration. Do not use `eval` to generate named
structs and global methods per registration; parametric carriers provide distinct concrete types
without world-age or redefinition problems.

## Open questions

- The recipes `entries` item descriptor (today's `DataItem` construction) needs its final name and
  shape.
- `project_fingerprint`: what enters the project-definition fingerprint and when it invalidates.
- `reconstruct` needs a concrete type, but `ItemRecord.kind` is a `Symbol` — for typed items,
  `Symbol(nameof(typeof(item)))`. Resolving `:PUNDLoop` back to `PUNDLoop` needs a module to look
  in. Candidates: persist a module-qualified name; have the project declare its item types; or
  accept rerun-from-source whenever the type cannot be resolved.

## Migration sequence

1. Declare `AbstractProject`, the project-aware stage forms with context-free defaults, typed
   collection `process`/`analyze`, and `reconstruct`.
2. Create `DataBrowserRecipes`: move `Project`, recipes, and `register_*` there; introduce the
   parametric carriers and immutable parametric recipes; implement the adapter methods; re-export
   the dialect from the `DataBrowser` umbrella.
3. Route the engine exclusively through the project-aware stage forms and one shared post-entries
   record normalization path.
4. Split the work graph, result states, profiling, and failure reporting into per-stage
   boundaries.
5. Remove the `isa RegisteredDataItem` cache gate and `cacheable(::AbstractDataItem)` with it;
   persistence becomes the storable-shape check, and rehydration becomes `reconstruct` dispatch.
6. Remove `data_items` and every registration branch outside `DataBrowserRecipes`; update public
   examples and documentation so an equivalent typed and registered project visibly traverse the
   same stages.
7. Add `register_csv!` as the first premade recipe.

## Acceptance

- The typed API exposes every engine pipeline stage available to registration.
- `DataBrowserCore`, `DataBrowserCache`, `DataBrowserSources`, and the GUI packages contain no
  reference to recipes, carriers, or `register_*`; only `DataBrowserRecipes` (and the umbrella
  re-export) knows the dialect exists.
- Registration callbacks execute only through adapter methods over the same stage runners and work
  keys as typed methods; no generic engine or cache path branches on `isa RegisteredDataItem`.
- The source is an argument to `read` only; every post-read stage is a pure function of records,
  source items, and values.
- After recipe lookup, stage runners infer the concrete recipe, callback, and carrier types;
  focused inference tests cover the function barriers.
- Re-registering a name replaces callbacks without accumulating methods or registrations.
- Typed collection process/analyze dispatch receives the current project-created collection value,
  while cache/index/GUI state continues to use `CollectionRecord`.
- A cached processed payload is delivered to views without running user code; `reconstruct` (where
  defined) rebuilds concrete items without rerunning `read` or `process`, and the rebuilt item's
  rederived identity is validated against the record.
- No item-level `cacheable` predicate remains; a typed item and a registered item with equivalent
  payloads persist identically.
- Full tests, docs, and both equivalent example styles pass without a compatibility layer for
  `data_items`.
