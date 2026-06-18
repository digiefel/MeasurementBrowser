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
  source's own subtype. They carry `item.data`, materialized only for the viewed selection.

```julia
abstract type AbstractDataItem end   # contract: id, item_label, kind, collection, parameters,
                                     # stats, item_data (+ optional process, cacheable, item_fingerprint)

struct DataItem <: AbstractDataItem  # the normal item; the recipe API's `entries` produces it
    id::String
    label::String
    kind::Symbol
    collection::Vector{String}
    parameters::Dict{Symbol,Any}
    stats::Dict{Symbol,Any}
    data::Any                        # the payload, reachable as `item.data`
end

struct ItemRecord                    # internal metadata record (never seen by source/project code)
    # source-item identity — which discovered unit produced this, and how to reload it
    source_item_id::String
    source_item_fingerprint::Any            # nothing → not persistently cacheable
    source_item_path::Union{String,Nothing}
    source_item_timestamp::Union{DateTime,Nothing}
    # logical item identity + metadata
    id::String                         # stable within its source
    item_label::String
    kind::Symbol
    collection::Vector{String}       # ["RuO2test", "A9", "VI", "D1"] — canonical tree placement
    collection_metadata::Dict{Symbol,Any}   # merged from metadata.txt (area_um2, t_HZO_nm, …)
    parameters::Dict{Symbol,Any}     # acquisition settings known while interpreting the unit
    stats::Dict{Symbol,Any}          # values computed after the required context is available
    item_fingerprint::Any            # nothing → payload not persistently cacheable
end
```

Source-*level* identity (`source_id`, `source_label`) is **not** copied onto every record — it lives
once on the `SourceScan`. A record carries only the source-*item* identity it needs to be reloaded.

## Scan Result

The completed scan is source-neutral. It stores source-level identity once, source-item fingerprints
for invalidation, and the hierarchy of `ItemRecord`s. No records are stored on the public
`SourceFile`:

```julia
struct SourceScan                     # the result of one full scan, source-neutral
    source_id::String
    source_label::String
    source_item_fingerprints::Dict{String,Any}
    hierarchy::Hierarchy
    analysis_failures::Vector{ItemFailure}
end
```

## Hierarchy

Records organize into a tree by their `collection::Vector{String}`. Nodes carry their own
collection-level stats (filled by the source's `collection_stats` hook / `register_collection_stat!`):

```julia
struct HierarchyNode
    name::String                       # one path segment
    kind::Symbol                       # :root | :level (intermediate) | :leaf (last segment)
    children::Vector{HierarchyNode}
    items::Vector{ItemRecord}          # records, populated only on :leaf nodes
    stats::Dict{Symbol,Any}            # collection-level stats for this node
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

## Collection metadata (path-prefix matching)

`metadata.txt` rows are keyed by path; metadata applies to **all descendants whose collection is
prefixed by that key**, with longer (more specific) keys overriding shorter ones, and is merged into
each record's `collection_metadata`. The on-disk format is documented in [storage.md](storage.md).
Collection metadata is a capability of the file-backed source; non-file sources need not provide it.

Implication: hierarchical metadata is stored once at the highest applicable level. Inheritance is
"lookup-time" — children don't physically carry parent metadata in their structs.

## Collection stats

Cross-item stats over one collection node are computed by the
`collection_stats(project, source, collection, items)` hook after the hierarchy is built, and stored
on the node's `stats` dict. The high-level `register_collection_stat!` is the callback form of the
same hook. This is distinct from per-item `stats`: item stats describe one item; node stats describe a
group.

## Virtual item expansion

A single source item may yield several data items when `data_items` returns more than one. They share
a source-item identity but must have distinct `id` values within the source. The recipe API mints
missing ids from the source-item path, kind, and the `parameters` that distinguish siblings.

Expansion is purely a source/project concern. For example, a source may split a multi-device sweep
into one item per device, or expand a fatigue file into one item per cycle (storing the cycle number
in `parameters`).

Item kinds and discovery rules belong to source/project implementations. The package stores the
resulting `Symbol` (`kind`) without assigning it project-specific meaning.

## Direct I-V data

DC current-voltage measurements use `:i` for current in amperes and `:v` for voltage in volts.
Voltage-driven sweeps, four-terminal measurements, breakdown measurements, and their visualizers
share these names so the same data can be reused across compatible views.
</content>
