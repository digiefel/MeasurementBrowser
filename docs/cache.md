# Source Cache

Scanning a source is inherently slow: the engine must enumerate every source item, interpret each one
into data items, and finish analysis before the final tree is known. The cache stores the result of
that work so the previous tree can appear immediately while a new scan checks the source.

The source remains authoritative. The cache is generated package data, stored outside the source, and
never exposed to source or project code.

## Startup

Opening a workspace starts cache loading and source scanning together.

Cache loading reads one compact index containing the previous `SourceScan`. It can populate the
workspace index immediately with the previous hierarchy, parameters, stats, skipped-item count,
fingerprints, and analysis failures.

The scan independently streams current items into the same workspace index and eventually produces
the authoritative `SourceScan`. The engine compares that scan with the cached index and updates the
cache when source items were added, changed, deleted, or produced a different analysis result.

## Identity

Cache identity is **source-based**, not root/project-based. Each source has one deterministic cache
id derived from `source_id(source)`:

```text
<source-label>-<source-id-digest>
```

Cache files live under:

```text
DEPOT_PATH[1]/measurementbrowser/cache/<source-label>/<cache-id>.h5
```

Scan/index ownership is keyed by `source_id` + `source_fingerprint`. A cache is accepted only when the
identity fields and the schema version match. For the file-backed `RegisteredProjectSource`,
`source_id` is the normalized root path, so this reproduces the previous root-based behavior exactly.

## Fingerprints and what they invalidate

| Fingerprint | Scope | Missing (`=== nothing`) means |
|---|---|---|
| `source_fingerprint(source)` | whole scan + all payloads | skip persistent scan/payload caching for this source |
| `source_item_fingerprint(item)` | records + payloads from one source item | skip persistent payload caching for that item's data |
| `item_fingerprint(item)` | one item's payload | skip persistent payload caching for that item |

A changed `source_item_fingerprint` invalidates only the records and payloads derived from that source
item. Item payloads are keyed by the full tuple:

```text
source_id + source_fingerprint + source_item_id + source_item_fingerprint + item_id + item_fingerprint
```

This generalizes the previous file-fingerprint behavior: for a file-backed source the source-item
fingerprint is the file's `stat` fingerprint, so changing or deleting one file invalidates exactly its
payloads while a non-file or live source (no fingerprints) simply runs without a persistent cache.

## Contents

The HDF5 file contains:

| Entry | Purpose |
|---|---|
| `schema_version` attribute | Rejects old cache layouts before deserialization |
| `/index` | One uncompressed serialized `ProjectCacheIndex` for fast startup |
| `/items/<source-item>/<item>` | Lazily cached payloads returned by the source's data load |

`ProjectCacheIndex` contains the completed `SourceScan` and analysis errors grouped by source item.

Payloads are grouped by source item, so changing or deleting one unit invalidates its payload group
regardless of how many logical items it contains.

## Annotations and collection metadata

> **TODO (source-agnostic storage, deferred):** annotations (`tags.txt`, `notes.txt`, `layout.txt`)
> and collection metadata (`device_info.txt`) are currently stored **next to the cache**, keyed by
> `source_id`, rather than at a filesystem source root. This is the simplest thing that works for a
> generic `AbstractDataSource` (a DB query or stream has no root). The regression is that these files
> are no longer hand-editable at the data folder. The intended fix is to make annotation storage a
> *source capability* — a file-backed source exposes hand-editable source-root files; other sources
> persist next to the cache — so this is a stopgap, not the final design.

## Freshness

A cached source item matches the source only when its `source_item_fingerprint` matches (for a file,
normalized path + size + mtime). The final scan comparison reports:

| Count | Meaning |
|---|---|
| fresh | cached source item matches the current scan and has no analysis error |
| stale | fingerprint or analysis result changed |
| new | present only in the current scan |
| deleted | present only in the cache |
| errors | current source items whose analysis failed |

Analysis failure does not remove items from the tree. It is stored as a valid cache difference so one
failed source item cannot invalidate the rest.

## Updates

A normal update uses the completed `SourceScan` already owned by the workspace. It never starts
another scan. It replaces the compact index, preserves payloads for unchanged source items, and
deletes payload groups for changed, deleted, or newly failing ones. A rebuild replaces the whole HDF5
file. Index writing does not load item data; direct and processed payloads are cached only when a
caller requests them.

## Item data

Item materialization loads payloads through the source's `load_data_item`. The workspace keeps
loaded items in memory and the cache layer can persist payloads with sufficient source-item and item
fingerprints.

Before reading or writing a payload, the engine checks the current fingerprints against the in-memory
cache index; stale payloads are never returned. A source item or item with a `nothing` fingerprint is
loaded fresh every time and never persisted. The loaded cache index belongs to its workspace — there
is no process-wide active cache. All HDF5 access uses one process-wide lock, so background cache
updates cannot overlap plot reads or lazy payload writes.
</content>
