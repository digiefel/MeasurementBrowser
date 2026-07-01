# Make CacheBuffer a generic mechanism; ProjectCache owns storage

## Context

The cache must make project work effectively free after it has already been completed. Submitting a
cache write should normally change memory and return. Reopening a project should reconstruct the
completed index without loading payload bodies. Reading one payload should issue only the physical
data query once its small location is known.

The cache also has to remain useful as the data-access layer grows:

- point reads must be fast;
- native tabular data must remain queryable by DuckDB;
- a future query engine must be able to select items by parameters or statistics and then read or
  aggregate their compatible payload columns without materializing every item in Julia;
- payload storage must remain practical at approximately one billion physical rows;
- opening and scanning projects containing approximately 100 GB of source data must not require
  reading cached payload bodies;
- vectors, matrices, and user-defined cacheable representations will be added later without exposing
  SQL or database concepts to project authors.

The current physical DataFrame layout is intended to support both point reads and global queries.
`item_data` is a small directory from an item to its physical location. DataFrames with the same
ordered column names and element types share one native DuckDB table. A compact sequence identifies
one item inside that table, so long item IDs are not repeated on every physical row.

One processed payload therefore occupies a tightly coupled table family:

```text
item_data                         one pointer per cached item
dataframe_schemas                 one description per physical DataFrame shape
dataframe_<shape>                 native rows for every item with that shape
```

These tables must have one owner because inserting, replacing, or deleting a payload changes them
together. A single payload buffer owns the complete family and commits one payload mutation batch in
one transaction. This is the deliberate exception to the usual one-buffer/one-table mapping:

> Every cache table has exactly one owning buffer. A buffer may own several tables only when they are
> one inseparable result store.

Ordinary buffers still own exactly one table.

The existing `CacheBuffer` implementation is not the design to preserve. It contains callback-based
reads, copied snapshots, revision reconciliation, reset state, and overlapping flush helpers. The
rewrite should remove that machinery rather than rename or reproduce it.

## Ownership invariants

1. `Workspace` and the rest of the application use semantic Cache operations only. They never run
   SQL, open a DuckDB connection, name cache tables, or access an individual buffer.
2. `CacheDB` owns the DuckDB database handle and every buffer. It owns no general read, write, or
   maintenance connection.
3. Every disk-backed `CacheBuffer` owns one persistent read connection and one persistent write
   connection. Those connections are never handed to callers.
4. Every cache table has one buffer owner. No second buffer or maintenance path may access it.
   ProjectCache's row-type-specific `read` and `_flush_to_db!` methods operate only as that owning
   buffer's implementation and use only its connections.
5. `CacheBuffer.jl` names no concrete ProjectCache row, table, column, cache stage, or payload format.
6. ProjectCache owns row types, SQL storage mappings, payload layout, EAV encoding, cache-index
   assembly, and semantic cascades.
7. `meta` is the only unbuffered data exception because it is tiny. Schema creation and the closing
   checkpoint are explicit setup/teardown operations, not alternate data-access paths.

## Buffer model

`CacheBuffer{R}` is live when constructed. A disk-backed buffer opens its two connections and starts
its background flush task before the constructor returns. A memory-only buffer has neither
connection nor flush task. `close!` permanently closes either kind; a buffer is never started,
stopped, paused, restarted, or run.

The buffer owns:

- `queued`: the mutations currently being accepted;
- `writing`: the detached batch currently being written;
- append/edit/delete coalescing inside `queued`;
- optional row-count capacity and backpressure;
- its read and write connections;
- its background flush task;
- memory-first reads;
- permanent closure.

There is no copied snapshot, mutation revision, buffer revision, reset revision, snapshot eviction,
or flushing state.

### Atomic queue transfer

At the start of a flush, the background task holds the buffer condition lock only long enough to:

```text
writing = queued
writing_rows = queued_rows
queued = new empty dictionary
queued_rows = 0
```

It performs no iteration, encoding, SQL, or connection work while holding that lock. New writes can
immediately enter the new `queued` dictionary while DuckDB writes `writing`.

Coalescing examines only `queued`. `writing` is the previous batch and is never merged back into the
new batch:

- append `row2`, then edit `row2` before transfer: `queued` contains the latest value as one append;
- transfer starts, then edit `row2`: `writing` retains the append and the new `queued` contains an
  edit for the next transaction.

Reads resolve the newest logical value in this order:

```text
queued → writing → committed payload location when configured → DuckDB
```

For a complete ordinary-table read, DuckDB rows are overlaid first with `writing`, then with `queued`.
Deletes remove the corresponding value. The newer queue always wins.

### Capacity

Capacity is counted in retained payload rows, not estimated bytes.

- ProjectCache supplies `row_limit` at each buffer construction.
- Bookkeeping buffers use `row_limit = nothing`.
- Disk-backed payload buffers block only when the active `queued` batch cannot accept another value.
- The batch already in `writing` does not count against backpressure because it no longer blocks
  producers from filling the next batch.
- One oversized value is accepted when `queued` is otherwise empty and becomes immediately eligible
  to flush.
- Memory-only buffers never wait; they reject an incoming value that would cross their limit.

`writing_rows` may exist for status reporting. It must not affect coalescing or backpressure.

## Mutation model

The hot mutation API is:

```julia
append!(buffer, key, row)
edit!(buffer, key, row)
delete!(buffer, key)
```

These operations change `queued` only and normally return without database I/O.

Coalescing preserves the cheapest final intent:

| Mutations in the same `queued` batch | Final queued intent |
|---|---|
| append → append/edit | append the latest row |
| append → delete | nothing |
| edit → edit | edit the latest row |
| edit → delete | delete |
| delete → append/edit | edit/replace |
| delete → delete | delete |

An edit submitted while an earlier append for the same key is already in `writing` remains an edit in
the new `queued` batch. The buffer does not inspect `writing` while coalescing.

### One flush operation

The generic buffer has one database mutation extension point:

```julia
_flush_to_db!(buffer, writing)
```

It is dispatched on the buffer row type and defined by ProjectCache. It uses the buffer's write
connection and applies the complete batch in exactly one explicit `BEGIN TRANSACTION`/`COMMIT`
transaction. Failure rolls that transaction back and is surfaced clearly.

For an ordinary table, one mixed transaction performs:

1. batched deletes;
2. bulk Appender writes for rows whose final intent is append;
3. one batched update/upsert/merge for rows whose final intent is edit;
4. commit.

Appends always use the Appender even when the same batch also contains edits. Appended rows are never
included in the edit query.

For the payload store, the one transaction may:

1. find the old physical locations for deletes and edits;
2. delete those old physical rows and directory entries;
3. register previously unseen DataFrame shapes;
4. append new physical rows, grouped by shape;
5. append replacement directory entries;
6. commit;
7. update the buffer's small in-memory location index before `writing` is cleared.

The payload buffer is the sole owner of every table touched by that transaction.

### Flush task

Keep only helpers with separate, documented jobs:

- one predicate may decide whether nonempty `queued` work is due because of its deadline, capacity,
  or closure;
- one background loop waits, transfers `queued` to `writing`, calls `_flush_to_db!`, reports failure,
  and clears `writing` after the transaction finishes.

Do not preserve `_flush_once!`, `_flush_and_wait!`, `_run_buffer`, forwarding aliases, or renamed
equivalents. Every retained helper must state whether its caller holds the condition lock and whether
it performs I/O.

The condition exists only to:

- protect the very short queue transfer and queue mutation sections;
- wake the flush task when new work arrives;
- wake disk-backed writers when transferring `queued` frees capacity;
- allow cold `clear!` or permanent `close!` to wait for `writing` to finish.

The read and write connections require no shared connection lock because they are different
connections.

## Reads

The only public buffer reads are:

```julia
read(buffer, key)
read(buffer)
```

There is no `read_from_db` field, callback supplied at construction, `_read_disk` API, or alternate
ProjectCache read path.

### Ordinary one-table stores

Each ordinary disk-backed buffer receives its table name and key columns:

```julia
CacheBuffer{ItemRow}(db, "items", ("id",))
CacheBuffer{MetaRow}(db, "metadata", ("scope", "entity", "key"))
```

`CacheBuffer` validates and quotes that configuration, constructs its two generic reads, and prepares
the keyed read on its own connection:

```sql
SELECT * FROM <table>
SELECT * FROM <table> WHERE <key1> = ? AND <key2> = ?
```

Rows returned by DuckDB are converted through exactly one pure constructor:

```julia
R(database_row)
```

That constructor converts values already returned by DuckDB. It performs no SQL, connection access,
or other I/O.

The configured key columns determine dictionary keys:

- one key column produces a scalar key;
- multiple key columns produce a tuple in configured order.

`CacheBuffer.jl` contains the generic query construction but no concrete table or column names.

### Payload store

Do not add a payload row or pointer wrapper. Delete `PayloadEntry`. The disk-backed payload buffer
queues the existing `DataItem` type:

```julia
CacheBuffer{DataItem}(
    db,
    "item_data",
    ("item_id",);
    row_limit=CACHE_BUFFER_ROW_LIMIT,
    disk_index=Dict{String,Tuple{UInt16,UInt32}}(),
)
```

`store_processed!` creates `DataItem(record, item_data(item))` and queues it by `record.id`.
`item_data` stores only processed payload pointers, so remove its redundant `stage` column.

The payload buffer owns `item_data`, `dataframe_schemas`, and every physical DataFrame table. It uses
the same queueing, connections, flushing, clearing, and closure as every other buffer. Only its
storage mapping differs.

`disk_index` is an in-memory copy of the `item_data` routing columns:

```julia
disk_index[item_id] = (storage_id::UInt16, seq::UInt32)
```

Persist `storage_id` as `USMALLINT`; `dataframe_schemas` maps that number to the DataFrame shape and
the physical table is named from the same number. No second in-memory shape catalog is needed.
`UInt16` permits 65,536 shapes. `UInt32` avoids the 65,536-version limit that replacement would
quickly make unsafe for sequences.

The complete payload read loads `disk_index` and the data-less processing-ready state from
`item_data`, then overlays `writing` and `queued`. It never loads payload bodies. The item pointer no
longer duplicates source and item fingerprints: validity belongs to the current cached `ItemRecord`
and work revision, and source replacement deletes the item's payload. Crash-atomic cross-buffer
cascades remain deferred.

A keyed payload read resolves `queued`, then `writing`, then `disk_index`. A queued or writing
`DataItem` supplies its resident `item_data`. A committed payload issues exactly one query against
its physical table:

```sql
SELECT * EXCLUDE (<internal sequence column>, <internal row-order column>)
FROM <physical table>
WHERE <internal sequence column> = ?
ORDER BY <internal row-order column>
```

The DuckDB result itself supplies the user column names and types, including for an empty result.
`dataframe_schemas` is not joined for every point read.

ProjectCache defines the payload-specific `read(::CacheBuffer{DataItem}, key)` method because it owns
the physical layout. This is the same public buffer operation, uses only that buffer's read
connection, and introduces no callback or second read API. `read_item_data` wraps a committed body as
`DataItem(record, body)`. Each successful payload flush updates `disk_index` before `writing` is
cleared, so point reads do not query `item_data` again.

### Result-state keys

Result-state buffer keys are the SQL-shaped tuple `(Int8(kind), entity)`. ProjectCache converts
between that tuple and `CacheResultKey` at its semantic boundary. Do not add a row-specific key rule
to the generic buffer.

### Future global queries

The query engine is not implemented in this phase, but this storage correction must preserve its
direct path:

1. use item parameters and statistics to select item IDs;
2. group their in-memory locations by compatible physical payload shape;
3. issue one DuckDB query per compatible shape;
4. join the selected item directory rows to the physical table and perform requested filtering,
   projection, or aggregation inside DuckDB.

It must never issue one directory query and one body query per selected item. Runtime should be
dominated by the requested physical rows, not by routing.

Vectors, matrices, custom cache representations, and the project-facing representation API remain
deferred. They must receive their own queryable physical layouts under the same single-owner rule;
this rewrite must not add a speculative codec framework.

## Cold clear and permanent close

Whole-store clearing is rare and synchronous. It does not participate in hot mutation state.

`clear!(buffer)`:

1. discards `queued`;
2. waits for any current `writing` transaction to finish;
3. holds the condition lock so no operation can enter;
4. clears the complete owned store in one transaction;
5. clears `disk_index` when present;
6. releases the lock.

The coordinator must stop producers before rebuilding and clearing buffers. Do not add `reset!`,
`queued_reset`, `writing_reset`, cancellation, transaction interruption, or a `clearing` lifecycle
state.

`clear_cache_index!` calls synchronous `clear!` on every buffer, then deletes `meta` through its
existing explicit short-lived ProjectCache connection. `meta` remains the sole unbuffered exception;
do not move it into a buffer or a general maintenance helper.

`close!` rejects new operations, wakes the background task, waits for its final flush, closes both
connections, and permanently rejects later use. Repeated closure is an explicit error.

## Remove connection and scheduling bypasses

Delete:

- `with_reader`;
- `with_reader_snapshot`;
- `with_maintenance`;
- `_with_maintenance`;
- `wait_cache_flushed!`;
- `_flush_and_wait!`;
- `cached_item_data_ids`;
- `_fill_item_data_request!`;
- `_item_data_cacheable`;
- `retain_source_items!`;
- the immediate maintenance deletion cascade.

`meta` reads/writes, schema setup, rebuild setup, and the closing checkpoint use explicit short-lived
ProjectCache connections only for those named operations. There is no general helper that hands a
connection to arbitrary work.

Rename `write_scan_identity!` to `write_meta_header!`.

Changed-source replacement is represented by buffered deletes followed by buffered appends/edits.
ProjectCache exposes semantic deletion operations; Workspace never names buffers or storage rows.

## Required edits

### `src/Cache/CacheBuffer.jl`

- remove `read_from_db`, copied snapshots, revisions, reset fields, and snapshot reconciliation;
- add explicit table/key configuration for ordinary disk-backed reads;
- use pure `R(database_row)` construction;
- replace the single pending dictionary with `queued` and `writing`;
- make queue transfer the only meaningful locked section during flush;
- coalesce only within `queued`;
- preserve append/edit/delete intent;
- count only `queued` for backpressure;
- retain one background loop and one `_flush_to_db!` call;
- implement synchronous cold `clear!` and permanent `close!`;
- retain no concrete ProjectCache names or SQL identifiers.

### `src/Cache/ProjectCache.jl`

- delete all six ordinary `_read_*` callbacks and the callback field wiring;
- construct ordinary buffers with their table and key columns;
- provide pure database-row constructors for every ordinary row type;
- keep the payload table family under the sole payload-buffer owner;
- delete `PayloadEntry` without replacing it and queue existing `DataItem` values;
- build and retain `disk_index::Dict{String,Tuple{UInt16,UInt32}}`;
- persist numeric `storage_id` values and remove the redundant `item_data.stage` and fingerprint
  columns;
- define the payload buffer's pointer-only complete read and one-query keyed body read;
- make the payload point read use `SELECT * EXCLUDE` rather than joining
  `dataframe_schemas`;
- use `(Int8(kind), entity)` tuple keys for the result-state buffer;
- split every flush by final intent so the Appender always receives append rows only and the edit
  statement receives edit rows only;
- make payload replacement and deletion one optimized transaction;
- move whole-store deletion into synchronous `clear!`;
- remove every general connection escape hatch and legacy readiness/maintenance API;
- keep semantic Cache operations as the only Workspace-facing surface.

### `src/Cache.jl`

- expose only the rewritten semantic Cache API needed by the later integration phase;
- remove exports and forwarding for deleted legacy operations;
- do not add compatibility aliases.

### `docs/cache.md`

- preserve the existing motivation and operational context;
- document the single-owner table-family rule;
- document queued/writing transfer, intent-aware flushes, and synchronous cold clear;
- document the in-memory payload locations and one-query point-read path;
- explain why native shape-grouped physical tables support future global aggregation;
- remove callback, snapshot/revision, reset, and maintenance descriptions.

## Implementation protocol

The cache correction is performed by the existing cache agent before work-graph implementation or
integration.

The current cache-buffer implementation is a counterexample, not a design to preserve. The agent
traces it only to identify behavior and code that must be removed. It follows this plan strictly,
prefers deletion and direct rewrites, and does not retain helpers, lifecycle states, connection paths,
compatibility aliases, or fallback behavior merely because current callers use them.

The cache agent owns only:

- `src/Cache/CacheBuffer.jl`;
- `src/Cache/ProjectCache.jl`;
- `src/Cache.jl` when required by the rewritten Cache module.

The primary agent owns this plan and `docs/cache.md`.

The cache agent must not:

- run any test, package load, benchmark, invariant grep, formatter, or validation command;
- edit tests or documentation;
- repair Workspace, ItemIndex, GUI, Browser, or other callers;
- add temporary compatibility methods to make the application load;
- weaken or reinterpret a requirement to accommodate existing code;
- continue after finding a real conflict between requirements.

If requirements conflict or an essential decision is absent, the agent stops immediately and reports
the exact conflict, affected symbols/files, and required decision. It does not guess.

This phase intentionally leaves the codebase broken. The application is not expected to load because
callers outside Cache still target removed APIs. The agent reports the rewritten Cache API, removed
API, remaining integration breakage, and every unresolved issue, then stops without committing.

The user and primary agent review the output together. Corrections are made only from explicit review
instructions. The work-dependency agent starts only after that review and after its separate plan has
been reviewed. That later agent completes integration, reviews the cache implementation in the
context of the complete refactor, and performs final validation as part of finishing the PR.
