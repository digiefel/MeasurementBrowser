# Source Cache

The cache is generated package data. Project code supplies source interpretation, processing, and
statistics; the package decides where and when their reusable results are stored.

## Location And Identity

Each project has one predictable DuckDB file:

```text
DEPOT_PATH[1]/measurementbrowser/<project-name>/cache.duckdb
```

The project name must be one safe path component and is used unchanged. The cache records the project
name and `source_id(source)`. Opening the same project name for a different source fails rather than
creating another cache. Old hash-named caches are ignored.

## Stored State

DuckDB stores source-item fingerprints, data-less item records, typed parameters and statistics,
failures, and two independent item-data stages:

- `interpreted`: data returned by `data_items` before `process`;
- `processed`: data returned by `process` and consumed by views.

The `item_data` key is `(item_id, stage)`, so the two forms cannot replace each other. Cacheable
`AbstractDataFrame` values use shared native DuckDB tables grouped by compatible ordered columns and
types. Other data types remain source-backed.

`cacheable(item)` is evaluated independently for the interpreted and processed values. The built-in
`DataItem` DataFrame path is cacheable; the low-level default is false. A false result removes any old
entry for that stage while preserving the item record.

## Refresh

Opening a workspace loads the last committed index first, then discovers the current source items.
An unchanged non-null source fingerprint reuses its records and cached stages. New, changed, and
unfingerprinted source items run the one interpretation path:

```text
data_items → normalize records/data → optional interpreted write → processing queue
```

Small groups of up to four ready source results share a transaction. This keeps progress steps small
without paying for one DuckDB transaction per item. Each committed group is independently reusable
after cancellation or interruption.

When interpreted data is absent, selection and background work use the same source fallback: find the
current source item and call `data_items` again. Every logical item returned by that source pass is
normalized and made available to its waiting processing job. There is no separate reload callback.

## Processing And Reads

The workspace owns one bounded processing queue. Background work and selected items use the same job;
selection promotes an existing waiting job instead of starting another path. Before queuing warm
background work, one metadata query finds all valid processed entries, so cached items are not checked
one at a time.

Processing reads interpreted data from DuckDB or source fallback, calls `process`, evaluates
`cacheable` on the processed item, and writes up to four simultaneously ready results together. Item
statistics follow processing and are stored separately. Collection statistics run after processing
settles.

DuckDB is the only package-level shared cache. There is no Julia object LRU. An active plot or
inspector owns the processed items it currently displays; a later selection reads DuckDB or repeats
the required upstream work. A plot selection is materialized once and the same item objects are
passed to setup and drawing.

One persistent writer connection serializes mutations. A separate persistent reader serves committed
item data. DuckDB's buffer pool is limited to 1 GiB by default through
`CACHE_MEMORY_LIMIT_MIB`; `set_cache_memory_limit!` changes the default or an open workspace.

## Invalidation And Recovery

Changing a source-item fingerprint replaces only that source item's records, both data stages, item
statistics, and dependent collection results. A missing fingerprint is treated as changed on every
refresh. There is no automatic project-code identity yet; project interpretation, processing,
statistics, or cache-policy changes require **Rebuild Cache**.

Corrupt or schema-incompatible generated caches are rebuilt automatically because no useful user
choice exists. A project/source identity conflict is different: it fails clearly because silently
replacing another source's valid cache would lose useful state.

Failures remain attached to the smallest source item, logical item, or collection that failed. They
do not invalidate unrelated committed work.
