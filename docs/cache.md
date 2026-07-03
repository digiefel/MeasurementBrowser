# Source Cache

Scanning a source is inherently slow. The package must discover source items, interpret each one into
logical items, process their data, and compute statistics. The source cache stores reusable results
from that work so a reopened workspace can publish its previous index immediately and avoid repeating
work whose inputs have not changed.

The cache is generated package data. Source files and project code remain the source of truth.
Project code supplies interpretation, processing, and statistics; the package decides where and when
their reusable outputs are stored. The cache may always be rebuilt, and its on-disk layout is an
internal implementation detail.

Opening a workspace does not discard an old cache by default. If the cache schema is out of date,
the workspace stays usable with a memory-only cache for that session. Script callers can opt in with
`open_workspace(...; rebuild=true)` to discard and rebuild the generated cache immediately. The
browser shows the same stale-cache condition as a prompt and only discards the generated cache when
the user chooses to rebuild.

The cache does not schedule work and does not own the published workspace index. It reports which
persisted results already exist. The work dependency graph decides what remains to be computed, and
`WorkspaceIndex` owns the records, hierarchy, parameters, statistics, and failures shown by the
application.

Cached item records contain item-local parameters only. Source-owned collection parameters are read
from the open source and applied to restored records exactly as they are in a memory-only workspace.

## Vocabulary

These terms have one meaning throughout the cache code:

- **Result store** — one independently buffered set of keyed values, such as source-item rows,
  logical-item rows, metadata, result states, or processed payloads. Most result stores map to one
  table; processed payloads span their pointer table, schema registry, and physical DataFrame tables.
- **Buffer** — bounded memory in front of one result store. It owns pending mutation intent, reads,
  capacity, its connections when disk-backed, and background flushing.
- **Append, edit, or delete** — a mutation submitted to a buffer. Submission changes memory only and
  normally returns without touching DuckDB.
- **Flush** — one buffer snapshot being applied to its result store in one DuckDB transaction.
- **Backpressure** — a disk-backed payload store waiting because its configured row limit is full.
  Bookkeeping stores are unlimited, and memory-only stores drop an incoming value instead of waiting.
- **Published data** — completed data in `WorkspaceIndex`. Buffered or persisted cache rows are not
  themselves the published workspace.

There is no separate producer/consumer framework. Package code submits semantic cache operations;
`CacheDB` translates them into typed buffer mutations.

## Ownership

`ProjectCache.jl` owns everything specific to MeasurementBrowser's project cache:

- `CacheDB` and the set of result stores it contains;
- row types, keys, and cache result kinds;
- the DuckDB schema, table names, and columns;
- source-item deletion cascades;
- cache identity and reconstruction of `ProjectCacheIndex`;
- semantic operations used by Workspace.

`CacheBuffer.jl` owns the fully generic `CacheBuffer{R}` mechanism:

- keyed pending append, edit, and delete intent;
- reset intent;
- coalescing repeated mutations;
- payload-row accounting and backpressure;
- persistent read and write connections for a disk-backed buffer;
- one background flush loop;
- memory-first reads;
- permanent closure.

The generic buffer names no concrete row type, table, column, cache stage, or ProjectCache type, but it can generate the generic SQL for any payload.
`CacheDB` owns the buffers and the DuckDB database handle, but it does not own a general-purpose
connection.

Workspace code uses semantic Cache operations. It never runs SQL, opens a DuckDB connection, names
cache storage, accesses an individual buffer, or decides how a result is written.

## One mechanism, two backing modes

### Disk-backed buffers

A disk-backed buffer opens one persistent read connection and one persistent write connection when it
is constructed, then immediately launches its background flush loop. Reads and flushes may proceed
concurrently because they never share a connection across tasks. There is no connection exposed to
ProjectCache callers and no maintenance connection hidden behind a callback.

Source items, logical items, metadata, failures, result states, and cacheable processed payloads are
disk-backed.

### Memory-only buffers

Some item data is useful only as a short-lived input to downstream work:

- interpreted item data waiting to be processed;
- processed data that the project declares non-cacheable.

These values use the same keyed buffer mechanism without a connection or flush task. They are never
written to DuckDB. If an incoming value would exceed the configured row limit, it is not retained.
A later miss causes the required upstream work to be performed again.

This prevents interpretation from being throttled by processing and prevents non-cacheable processed
data from becoming an unbounded Julia object cache.

## Capacity and backpressure

Capacity is counted directly in retained payload rows. ProjectCache supplies the row limit when it
constructs each payload buffer and defines how that payload type reports its DataFrame row count.
There is no estimated byte weight.

Bookkeeping buffers have no row limit. They still report the number of pending keys for workspace
status, but they never call the payload row-count operation and never backpressure.

A disk-backed payload submission waits only when accepting it would exceed the row limit and another
pending value can be flushed to make room. One oversized value is accepted when the buffer is
otherwise empty and is made immediately eligible for flushing, so it cannot deadlock.

A memory-only payload submission never waits. If replacing or adding the value would exceed its row
limit, the incoming value is not retained.

## Mutation intent

The buffer exposes four mutation operations:

```julia
append!(buffer, key, row)
edit!(buffer, key, row)
delete!(buffer, key)
reset!(buffer)
```

They update pending memory only. Repeated operations for one key coalesce to the final useful intent:

| Pending sequence | Final pending operation |
|---|---|
| append followed by append or edit | append the latest value |
| append followed by delete | no operation |
| edit followed by edit | edit to the latest value |
| edit followed by delete | delete |
| delete followed by append or edit | edit or replace |
| delete followed by delete | delete |

Preserving this distinction matters. A pure append batch can use DuckDB's bulk Appender. An edit to
an ordinary keyed table can use a batched update, upsert, or merge without first deleting every row.
A processed-payload replacement is different: its physical shape and row count may change, so its
transaction removes the old physical rows before appending the replacement and its pointer.

`reset!` discards earlier pending operations, makes reads ignore committed rows, and marks the next
flush to clear the complete result store before applying mutations submitted after the reset.

## Flushing

Each disk-backed buffer has one background loop. Pending work becomes flushable when:

- the existing flush deadline is reached;
- the payload row limit is reached;
- the buffer is closing.

Continued submissions do not move the existing deadline. A flush snapshots the already-coalesced
intent and calls one ProjectCache extension point:

```julia
_flush_to_db!(buffer, snapshot)
```

The ProjectCache method inspects the whole snapshot and chooses the fewest suitable batched statements
for that result store. It uses `buffer.write_connection` and performs exactly one transaction
containing the snapshot's reset, appends, edits, and deletes. There are no separate append, edit, or
delete flush hooks.

Submissions may continue while database I/O runs. After a successful commit, the buffer removes only
the exact mutation revisions present in the snapshot. A newer mutation for the same key remains
pending.

Cache outputs are reproducible. If a background flush fails, the error is reported and that snapshot
is discarded so the flush task remains alive and a blocked submission cannot wait forever. Missing
outputs are recomputed when demanded. Cross-buffer crash atomicity and automated cache repair are
separate concerns.

## Reads

The only generic read operations are:

```julia
read(buffer, key)
read(buffer)
```

ProjectCache supplies a result store's disk-reading behavior when it constructs the buffer. That
behavior remains private to the buffer and cannot be called as a side-door disk read.

A keyed read observes pending reset and mutation state before consulting DuckDB, then rechecks pending
state after the database read. A complete read obtains committed rows and overlays pending
appends and edits, removes pending deletes, or starts from an empty committed view while reset is
pending.

The full lookup order for disk-backed processed data is:

```text
pending buffer value → DuckDB value → required upstream work
```

For memory-only interpreted or processed data it is:

```text
pending buffer value → required upstream work
```

Loading `ProjectCacheIndex` reads data-less rows, processed-payload pointers, and result states. It
does not materialize processed payload bodies merely to determine which work is already satisfied.

## Lifecycle and connections

A buffer is live when constructed. It is never started, stopped, paused, or restarted.

`close!` is the only lifecycle transition. It rejects new operations, wakes the flush loop, drains the
final snapshot, waits for the loop to finish, closes the owned read and write connections, and
permanently closes the buffer. Repeated closure is an error.

The tiny `meta` header is the only unbuffered data. ProjectCache uses explicit, short-lived
connections for:

- initial schema creation and validation before buffer construction;
- reading, writing, or clearing the meta header;
- the final checkpoint after every buffer has closed.

These are named operations, not a general connection callback. There is no `with_reader`,
`with_maintenance`, or equivalent escape hatch.

## Stored state

The fixed DuckDB stores are:

- `meta` — schema version plus project and source identity;
- `source_items` — source-item ID, fingerprint, path, and timestamp;
- `items` — data-less logical-item records and their source ownership;
- `metadata` — typed item parameters, item statistics, collection parameters, and collection
  statistics;
- `item_failures` — source interpretation and logical-item failures published with the index;
- `result_states` — independent processing, item-statistics, and collection-statistics outcomes;
- `item_data` — pointers to cacheable processed payloads;
- `dataframe_schemas` — the ordered user-column names belonging to each payload shape;
- per-shape DataFrame tables — native rows of processed tabular payloads.

Queued and running work is not persisted. It belongs to the work dependency graph.

### Persisted result state

Successful empty statistics produce no metadata rows, so absence of metadata cannot mean "not
computed." `result_states` therefore records:

- processing failure;
- processing success for a persisted payload;
- item-statistics success, including an empty result;
- item-statistics failure;
- collection-statistics success, including an empty result;
- collection-statistics failure.

Each state records the effective-input fingerprint that produced it. Non-cacheable processed success
remains memory-only.

These outcomes are independent. Processing failure prevents item statistics for that item revision,
but an item-statistics failure does not erase a successfully processed payload. Collection-statistics
failure does not alter member item statistics.

### Processed DataFrames

Cacheable processed `AbstractDataFrame` values are stored in native DuckDB tables. One physical table
is used for each ordered combination of user column names and element types. The shape's stable digest
becomes an internal storage ID and table name; user text is never used directly as a table name.

Each payload receives a compact integer sequence. Its `item_data` pointer and all physical rows share
that sequence instead of repeating the logical item ID in every row. Two reserved internal columns
store the sequence and the row order. User-derived identifiers are quoted, and the reserved internal
prefix cannot collide with user columns.

The pointer, schema registration when needed, and physical payload rows are committed in one
transaction. A reader cannot observe a committed pointer without its payload. A zero-row DataFrame
still has a pointer and schema, so it reconstructs as an empty value with the correct columns and
types.

Stored input fingerprints are validated once, when the cache index is loaded: an entry whose
fingerprint matches the current source-item identity, item identity, and effective parameters seeds
a completed node in the workspace work graph; a mismatch seeds nothing and the work is re-enqueued
(the interpreted item record remains reusable). At runtime the work graph is the single freshness
authority — readers consult it and then read payloads raw; `result_states` is only written, never
re-read, until the next load.

## Source changes

ProjectCache exposes one semantic cascade for removing a source item:

```julia
delete_source_item!(cache, source_item_id, old_records)
```

It derives item IDs and item-scoped metadata, failures, result states, disk payloads, and memory-only
values from the old published records. It does not decide which collection statistics are affected.

The work dependency layer separately identifies affected collection paths and calls the semantic
`delete_collection_stats!(cache, collection_keys)` operation. A changed source is represented as
buffered deletion of its old subtree followed by buffered stores for the replacement. Per-key
coalescing turns retained keys into replacement edits before any database transaction is chosen.

No maintenance cascade writes directly to DuckDB.

## Cache identity

Each project has one predictable cache file:

```text
DEPOT_PATH[1]/measurementbrowser/<project-name>/cache.duckdb
```

The project name must be one safe path component and is used unchanged. The meta header binds the file
to the project name and `source_id(source)`. A cache belonging to a different source fails clearly
instead of silently replacing valid data.

Old schema migration is not supported by this refactor. An incompatible generated cache fails with a
specific error and can be explicitly rebuilt from its source.

DuckDB's buffer pool is configured when a cache is opened. `CACHE_MEMORY_LIMIT_MIB` controls the
default for subsequently opened caches.

## Concurrency assumptions

DuckDB permits concurrent connections inside one process and append transactions do not conflict.
Updates or deletes of the same rows can conflict.

The cache avoids ambiguous ownership structurally:

- each disk-backed result store has one owning buffer with separate persistent read and write
  connections;
- only that buffer mutates the store;
- one flush contains the complete coalesced snapshot transaction;
- payload pointer, schema, and physical-row mutations are coordinated by the single payload buffer;
- fixed schema creation finishes before buffer connections are opened;
- meta and checkpoint operations are explicit and do not run as alternate writers to buffered stores.

The design does not depend on a hidden global writer, a maintenance connection, or callers waiting for
buffers to flush before they can read the logical cache state.
