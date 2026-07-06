# Data Model

The model has three core source/data layers and one internal record, in increasing specificity:

```text
AbstractDataSource      → owns lifecycle + discovery               (a dataset root, DB query, stream)
AbstractDataSourceItem  → one discovered unit                      (a file, row, run id, channel)
AbstractDataItem        → one logical browsable item               (what you index/select/plot)
ItemRecord              → internal, data-less index/cache record   (never seen by source/project code)
```

`DataSources/DirectorySource.jl` defines the built-in file-backed source item (`SourceFile`).
`ItemIndex` defines the internal `ItemRecord`, the concrete `DataItem`, and the collection hierarchy.
The project-facing API is in [api.md](api.md); the design rationale is in
[plans/data-model-generalization.md](plans/data-model-generalization.md).

## Source item vs. data item

A **source item** is one addressable thing a source discovers — the unit of scan progress, failure
reporting, and invalidation. A **data item** is one logical browsable object. The mapping is
one-to-many: `data_items(project, source, source_item)` interprets a single unit into zero, one, or
many data items. A `SourceFile` is a source item; the items a project reads out of it are data items.

## Two representations: record vs. item

Metadata and data live in two forms, bridged by the engine via the `AbstractDataItem` contract:

- **`ItemRecord`** — the **internal**, data-less metadata record the hierarchy, scan, and cache store.
  Source and project code never construct or name it.
- **`AbstractDataItem`** instances — what callbacks and views see: the package's `DataItem` or a
  source's own subtype. They carry interpreted data into the work graph and processed data into
  views. Later selections materialize the processed stage from DuckDB or repeat the required
  interpreted/processing work through the normal source path.

```julia
abstract type AbstractDataItem end   # contract: id, item_label, kind, collection, metadata,
                                     # item_data (+ optional process, cacheable, item_fingerprint)

struct DataItem <: AbstractDataItem  # the normal item; the recipe API's `entries` produces it
    id::String
    label::String
    kind::Symbol
    collection::Vector{String}
    metadata::Dict{Symbol,Any}       # one dict: parsed parameters + computed values, provenance merged
    data::Any                        # loaded item data, reachable as `item.data`
end

struct ItemRecord                    # internal metadata record (never seen by source/project code)
    # source-item identity — which discovered unit produced this, and how to reload it
    source_item_id::String
    source_item_path::Union{String,Nothing}
    source_item_timestamp::Union{DateTime,Nothing}
    # logical item identity + metadata
    id::String                         # stable within its source
    item_label::String
    kind::Symbol
    collection::Vector{String}       # ["RuO2test", "A9", "VI", "D1"] — canonical tree placement
    metadata::Dict{Symbol,Any}       # the entries layer: metadata from interpreting this item
    item_fingerprint::Any            # nothing → item data not persistently cacheable
end
```

An item carries **one `metadata` dict**. Internally the engine tracks provenance in layers — the
entries layer (`ItemRecord.metadata`, written at interpretation and restored on reload), inherited
collection metadata, and the computed layers from item `analyze` and collection `process` — but
project code only ever sees the merged result. Delivery is one unconditional merge: inherited
collection metadata ⊕ entries ⊕ computed layers, deeper/newer wins.

Source-*level* identity (`source_id`, `source_label`) is **not** copied onto every record — it lives
once on the `SourceScan`. A record carries only the source-*item* identity it needs to be reloaded.

## Scan Result

The completed scan is source-neutral. It stores source-level identity once and the hierarchy of
`ItemRecord`s. No records are stored on the public `SourceFile`:

```julia
struct SourceScan                     # the result of one full scan, source-neutral
    source_id::String
    source_label::String
    hierarchy::Hierarchy
    analysis_failures::Vector{ItemFailure}
end
```

Source-item fingerprints (for reopen invalidation) live only in the cache `source_items` table.
`scan_source!` loads them directly at reopen; a memory-only cache holds none, so every discovered
item re-interprets. Collection-metadata freshness at reopen diffs cached `source_collection_metadata` against current
`collection_metadata(source, path)` for each hierarchy path.

## Hierarchy

Records organize into a tree by their `collection::Vector{String}`. Each node keeps two dicts: the
source-provided `metadata` (re-applied wholesale on source updates) and the `analysis` output of
collection `analyze` (filled by workspace background analysis). The GUI shows them merged:

```julia
struct HierarchyNode
    name::String                       # one path segment
    kind::Symbol                       # :root | :level (intermediate) | :leaf (last segment)
    metadata::Dict{Symbol,Any}         # effective source-provided collection metadata at this node
    children::Vector{HierarchyNode}
    items::Vector{ItemRecord}          # records, populated only on :leaf nodes
    analysis::Dict{Symbol,Any}         # collection `analyze` output for this node
end

struct Hierarchy
    root::HierarchyNode
    all_items::Vector{ItemRecord}
    source_id::String
    index::Dict{Tuple{Vararg{String}}, HierarchyNode}   # path tuple → node
    has_collection_metadata::Bool
    skipped_count::Int
end
```

## Variable hierarchy depth

Depth is **not** uniform. One branch may be `RuO2test/A9/VI/D1` (4 segments) while another is
`RuO2test_A10/VI/25um/D1` (also 4) and another could be 3 or 5. Don't assign semantic meaning to a
fixed depth.

- All items attach at `:leaf` nodes (the last segment).
- Intermediate nodes are `:level`.
- A leaf is **the last segment of its path**, not "depth N from root."

## Identity / path keys

An item's `collection::Vector{String}` is its canonical placement. Two equivalent representations:

- **Tuple** (used for `index` lookup): `("RuO2test", "A9", "VI", "D1")`.
- **String** (used in on-disk files): `"RuO2test/A9/VI/D1"` — built by `collection_path_key(collection)`.
- **Round-trip**: `collection_path_tuple("RuO2test/A9/VI/D1")` parses and validates the string form.

The slash-joined string is the canonical identifier in any stored metadata. New metadata systems
should reuse it.

## Collection metadata

Sources may provide collection metadata for a collection path through
`collection_metadata(source, collection_path)`. The hierarchy owns inheritance: parent collection
metadata flows to child nodes. `ItemRecord.metadata` stores only the item's entries layer. The
workspace materializes each item's effective metadata by merging the containing hierarchy node with
the item-local values; item values win on key conflicts.

For `DirectorySource`, `metadata.txt` rows are keyed by slash-joined collection path fragments. The
source owns loading, watching, parsing, and excluding that configured file from source discovery.
Metadata changes re-merge effective item inputs without reinterpreting unchanged source files. The
on-disk format is documented in [storage.md](storage.md).

## The metadata pipeline

Each item runs a fixed staged pipeline, each stage fed the previous output:

```text
entries → item process → item analyze → collection process → collection analyze
```

`entries` seeds the item's initial metadata and data. `process` may change per-item data or metadata.
`analyze` returns a dict merged over the item's computed layer (new keys win). The collection stages
are gated on every member of a registered kind finishing its item stages, and group by collection
path plus registered kind:

- **collection `process`** is the down-flow: it rewrites each member (one output per input, same ids),
  overwriting per-item data or metadata for the whole collection. Returning the same `data` object
  skips the payload rewrite, so metadata-only recipes cost nothing extra.
- **collection `analyze`** folds the post-process members into one dict attached to the collection
  node's `analysis` only. It does not flow down to items.

Re-running a stage replaces that stage's previous output; the internal layering keeps stale keys from
accumulating.

## Freshness

Freshness comes from source data only. In-session it is event-driven: a source upsert/removal or a
collection-metadata change propagates through the work graph, bumping the affected keys' revisions
(plain per-key counters) and re-marking their results stale; recompute overwrites the stored rows.
At reopen the same mechanism replays from a source-fingerprint diff — the scan compares per-source-item
fingerprints and the source's collection-metadata fingerprints against the stored `SourceScan`, and
marks the subtrees downstream of anything changed. A result's validity is derived from its position
downstream of unchanged sources, not from a stored per-result claim. Callback code is not tracked yet:
an edited callback leaves stale stored results until the next recompute or rebuild.

## Virtual item expansion

A single source item may yield several data items when `data_items` returns more than one. They share
a source-item identity but must have distinct `id` values within the source. The recipe API mints
missing ids from the source-item path, kind, and the `metadata` that distinguishes siblings.

The source item is read once for the expansion. Its children enter the workspace work graph and may
run on different scheduler threads. A child may carry a view into the parsed data instead of a copy;
the package keeps the interpreted parent alive until processing finishes. Because siblings may run
concurrently, overlapping views must be treated as read-only. DataFrames restored from the cache are
independent ordinary `DataFrame`s.

Expansion is purely a source/project concern. For example, a source may split a multi-device sweep
into one item per device, or expand a fatigue file into one item per cycle (storing the cycle number
in `metadata`).

Item kinds and discovery rules belong to source/project implementations. The package stores the
resulting `Symbol` (`kind`) without assigning it project-specific meaning.

## Direct I-V data

DC current-voltage measurements use `:i` for current in amperes and `:v` for voltage in volts.
Voltage-driven sweeps, four-terminal measurements, breakdown measurements, and their visualizers
share these names so the same data can be reused across compatible views.
</content>
