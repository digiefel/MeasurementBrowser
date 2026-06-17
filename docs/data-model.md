# Data Model

`ItemIndex` defines physical source files, the internal item record, the concrete data item, and the
collection hierarchy. Projects interpret indexed source files; package code stores and presents the
resulting items. The project-facing API is in [api.md](api.md); the design rationale is in
[plans/data-model-generalization.md](plans/data-model-generalization.md).

## Two representations: record vs. item

Metadata and data live in two forms, bridged by the engine via the `AbstractDataItem` contract:

- **`ItemRecord`** — the **internal**, data-less metadata record the hierarchy, scan, and cache store.
  Project code never constructs or names it.
- **`AbstractDataItem`** instances — what callbacks see: the package's `DataItem` or a project's own
  subtype. They carry `item.data`, materialized only for the viewed selection.

```julia
abstract type AbstractDataItem end   # the contract: item_id, item_label, kind, collection,
                                     # parameters, stats, item_data, read_data (+ optional process, cacheable)

struct DataItem <: AbstractDataItem  # the normal item; register_item!'s `entries` produces it
    unique_id::String
    label::String
    kind::Symbol
    collection::Vector{String}
    parameters::Dict{Symbol,Any}
    stats::Dict{Symbol,Any}
    data::Any                        # the payload, reachable as `item.data`
end

struct ItemRecord                    # internal metadata record (never seen by project code)
    unique_id::String                # stable identity for one logical item
    filename::String
    filepath::String
    clean_title::String
    kind::Symbol
    timestamp::Union{DateTime,Nothing}
    collection::Vector{String}       # ["RuO2test", "A9", "VI", "D1"] — canonical tree placement
    collection_metadata::Dict{Symbol,Any}   # merged from device_info.txt (area_um2, t_HZO_nm, …)
    parameters::Dict{Symbol,Any}     # acquisition settings known while interpreting the file
    stats::Dict{Symbol,Any}          # values computed after the required context is available
end

struct SourceFile
    unique_id::String
    filepath::String
    filename::String
    timestamp::Union{DateTime,Nothing}
    fingerprint::FileFingerprint
    measurements::Vector{ItemRecord}   # the records interpreted from this file
end

struct HierarchyNode
    name::String                       # one path segment
    kind::Symbol                       # :root | :level (intermediate) | :leaf (last segment)
    children::Vector{HierarchyNode}
    measurements::Vector{ItemRecord}   # records, populated only on :leaf nodes
end

struct Hierarchy
    root::HierarchyNode
    all_measurements::Vector{ItemRecord}
    root_path::String
    index::Dict{Tuple{Vararg{String}}, HierarchyNode}   # path tuple → node
    has_collection_metadata::Bool
    project::Project
    skipped_count::Int
end
```

(The hierarchy fields are still named `measurements`/`all_measurements`; they hold `ItemRecord`s.)

## Variable hierarchy depth

Depth is **not** uniform. One branch may be `RuO2test/A9/VI/D1` (4 segments) while another is `RuO2test_A10/VI/25um/D1` (also 4) and another could be 3 or 5. Don't assume "level 3 = device".

- All items attach at `:leaf` nodes (the last segment).
- Intermediate nodes are `:level`.
- A leaf is **the last segment of its path**, not "depth N from root."

## Identity / path keys

An item's `collection::Vector{String}` is its canonical placement. Two equivalent representations:

- **Tuple** (used for `index` lookup): `("RuO2test", "A9", "VI", "D1")`.
- **String** (used in on-disk files): `"RuO2test/A9/VI/D1"` — built by `collection_path_key(collection)`.
- **Round-trip**: `collection_path_tuple("RuO2test/A9/VI/D1")` parses and validates the string form.

The slash-joined string is the canonical identifier in any on-disk metadata file. New metadata systems should reuse it.

## Path-prefix matching

`device_info.txt` rows are keyed by path; metadata applies to **all descendants whose collection is prefixed by that key**, with longer (more specific) keys overriding shorter ones, and is merged into each record's `collection_metadata`.

The on-disk format is documented in [storage.md](storage.md). The scan applies matching parameters
before records enter the hierarchy.

Implication: hierarchical metadata is stored once at the highest applicable level. Inheritance is "lookup-time" — children don't physically carry parent metadata in their structs.

## Virtual item expansion

A single source file may yield several items when the project's `entries` callback returns a vector
with more than one entry. They share a filepath but get distinct `unique_id` values, which the engine
mints from filepath + kind + the `parameters` that distinguish siblings.

Expansion is purely a project concern. For example, a project may split a multi-device sweep into one
entry per device, or expand a fatigue file into one entry per cycle (storing the cycle number in
`parameters`).

Item kinds and filename rules belong to project implementations. The package stores the resulting
`Symbol` (`kind`) without assigning project-specific meaning to it.

## Direct I-V data

DC current-voltage measurements use `:i` for current in amperes and `:v` for voltage in volts.
Voltage-driven sweeps, four-terminal measurements, breakdown measurements, and their visualizers
share these names so the same data can be reused across compatible views.
