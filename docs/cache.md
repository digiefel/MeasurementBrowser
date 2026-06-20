# Source Cache

The cache makes a workspace useful before a current source scan finishes and prevents repeated
source reads after data has already been interpreted. It is generated package data stored outside
the source and is never exposed to project code.

## Startup And Updates

Opening a workspace loads the existing DuckDB index while discovering the current source items. A
valid cached index can populate the previous hierarchy, parameters, stats, fingerprints, and failures
immediately.

Current source items are processed concurrently. One worker owns one source item until it has:

1. produced its logical data items,
2. run `process` and item stats,
3. returned the completed records and cacheable loaded items.

The scan collector publishes records immediately and commits up to one source item per worker in one
DuckDB transaction. Its bounded result queue prevents completed loaded batches from accumulating for
the whole source.
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
| `fingerprint(source_item)` | records and cached item data produced by that source item | do not persist its item data |
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
| `payloads` | currently serialized cacheable loaded items, validated by fingerprints |

The `payloads` name is internal storage terminology. The public model is `item_data(item)` plus the
existing `cacheable(item)` and `fingerprint(item)` hooks.

## Reads And Concurrency

Item materialization checks, in order:

1. the workspace's bounded in-memory item cache,
2. valid cached item data in DuckDB,
3. `load_data_item` from the source.

One persistent writer connection serializes mutations. A separate persistent reader serves
interactive item-data reads from committed snapshots. Item-data reads are requested in one joined
query; batch writes use one transaction.

## Status

The current status model compares cached and current source items as fresh, stale, new, deleted, or
failed. Failures remain attached to their source/item and do not invalidate unrelated data.
