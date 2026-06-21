# Source Cache

The cache makes a workspace useful before a current source scan finishes and prevents repeated
source reads after data has already been interpreted. It is generated package data stored outside
the source and is never exposed to project code.

## Startup And Updates

Opening a workspace loads the existing DuckDB index while discovering the current source items. A
valid cached index can populate the previous hierarchy, parameters, stats, fingerprints, and failures
immediately.

Discovery compares each current source-item fingerprint with the cached one. An unchanged
fingerprinted source item reuses its cached records without reading or analyzing the source again;
new, changed, and unfingerprinted source items are interpreted normally. A fully unchanged source
returns the cached scan directly. Sources with collection parameters currently re-interpret every
source item because those parameters are folded into stored item records.

Current source items are processed concurrently. One worker owns one source item until it has:

1. produced its logical data items,
2. run `process` and item stats (across scheduler threads for a large expansion),
3. returned the completed records and cacheable loaded items.

The scan collector publishes records immediately and stages up to one source item per worker at a
time. Each worker-sized batch is committed independently. This bounds both loaded Julia data and
DuckDB's transaction-owned write memory; a canceled or interrupted build can reuse the source items
already committed. The bounded result queue prevents loaded items from accumulating for the whole
source.
After all source items finish, collection stats run from data-less items containing the completed
parameters and item stats.

## Identity

Cache identity is source-based. The deterministic id is derived from `source_id(source)`:

```text
<source-label>-<source-id-digest>
```

Cache files live under:

```text
DEPOT_PATH[1]/measurementbrowser/cache/<source-label>/<cache-id>.duckdb
```

A cache is accepted only when its source identity and schema version match.

## Invalidation

| Fingerprint | Scope | Missing (`nothing`) means |
|---|---|---|
| `fingerprint(source_item)` | records and cached item data produced by that source item | re-read it every scan and do not persist its item data |
| `fingerprint(data_item)` | cached data for that logical item | use only source-item invalidation |

Changing one source-item fingerprint deletes only the records and cached item data derived from that
unit. Parameter-only changes update metadata without discarding still-valid item data.

## DuckDB Contents

| Table | Purpose |
|---|---|
| `meta` | schema and source identity |
| `source_items` | source-item fingerprints, locations, and errors |
| `items` | data-less logical item records |
| `metadata` | typed item/collection parameters and stats |
| `item_data` | item fingerprints, native storage ids, and row counts |
| `dataframe_schemas` | ordered user column names for each DataFrame shape |
| `dataframe_<storage-id>` | native columnar rows for all cached DataFrames with one shape |

Standard `DataItem` values carrying any `AbstractDataFrame`, including row views returned by an
expanding `entries` callback, are registered with DuckDB and copied directly into shared columnar
tables. Cached reads return ordinary independent `DataFrame`s. Cache writes stage at most one source
item per scan worker at a time, so an expanded source is written once without allowing loaded data
to accumulate for the whole scan. Julia serialization is not involved. The public model remains
`item_data(item)` plus the existing `cacheable(item)` and `fingerprint(item)` hooks.
Built-in `DataItem`s carrying other data types remain source-backed until that type has a native
cache implementation.

## Reads And Concurrency

Item materialization checks, in order:

1. the workspace's bounded in-memory item cache,
2. valid cached item data in DuckDB,
3. `load_data_item` from the source.

One persistent writer connection serializes mutations. A separate persistent reader serves
interactive item-data reads from the last committed snapshot. Item-data reads are requested in one
joined query; scans use bounded worker-sized write transactions. The cache database has a 512 MiB
buffer-memory limit, so DuckDB evicts committed table blocks instead of growing toward its default
system-wide allowance as the cache becomes large. Julia objects and DuckDB allocations outside its
buffer manager remain visible in the process RSS separately.

## Status

The current status model compares cached and current source items as fresh, stale, new, deleted, or
failed. Failures remain attached to their source/item and do not invalidate unrelated data.
