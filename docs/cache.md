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

Cached item records carry the entries layer of their metadata. Collection nodes separately persist
the package-owned label, metadata, and registration-name projections needed to reconstruct the
hierarchy and inherited collection metadata after reopen. Concrete project collection values do not
enter the cache.

## Vocabulary

These terms have one meaning throughout the cache code:

- **Result store** — one independently buffered set of keyed values, such as source-item rows,
  logical-item rows, metadata, result states, or processed payloads. Most result stores map to one
  table; processed payloads span their pointer table, schema registry, and physical payload tables.
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

The cache is its own package, `lib/DataBrowserCache/`, depending only on `DataBrowserAPI` (for the
item contracts and `ItemIndex` types it reconstructs on reopen), `DataBrowserProfiling`, and its
storage backend (DuckDB/DBInterface). `DataBrowserCore`'s `Workspace` consumes it. The cache never
requires a specific table container: payloads come in as anything implementing Tables.jl (the
`cacheable_data` default is `Tables.istable`) and come back out as the container type they were
stored with. A non-tabular type can opt in by dispatching `cacheable_data` and implementing the
interface.

`project_cache_domain.jl` owns everything specific to DataBrowser's project cache:

- `CacheDB` and the set of result stores it contains;
- row types, keys, and cache result kinds;
- the DuckDB schema, table names, and columns;
- source-item deletion cascades;
- cache identity and reconstruction of `ProjectCacheIndex`;
- semantic operations used by Workspace.

`cache_buffer.jl` owns the fully generic `CacheBuffer{R}` mechanism:

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
constructs each payload buffer and defines how that payload type reports its row count.
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

- `meta` — schema version, project and source identity;
- `source_items` — source-item ID, fingerprint, path, and timestamp;
- `items` — data-less logical-item records, their source ownership, and each item's integer `item_key`
  surrogate plus the private key of its leaf collection;
- `collections` — one row per collection record, with its compact integer key, parent key, durable
  occurrence ID, resolved label and metadata, and registration name when applicable; arbitrary
  user-defined collection values and canonical ID inputs are not stored;
- `source_item_metadata` — the entries layer per item, keyed by `item_key`; reload restores
  `ItemRecord.metadata`;
- `analyzed_item_metadata` — the delivered metadata dict per item (inherited ⊕ entries ⊕ computed
  layers), keyed by `item_key`; this is the query surface;
- `analyzed_collection_metadata` — collection `analyze` output per collection, keyed by collection key;
- `wide_columns` — the registry naming each wide table's columns and their logical type;
- `item_failures` — source interpretation and logical-item failures published with the index;
- `result_states` — independent item-processing and item-analysis outcomes keyed by item identity;
- `collection_result_states` — collection-process and collection-analysis outcomes keyed by the
  private integer collection-record key;
- `item_data` — pointers to cacheable processed payloads, keyed by `(item_key, stage)`;
- `dataframe_schemas` — the ordered user-column names belonging to each payload shape;
- per-shape payload tables — native rows of processed tabular payloads.

Queued and running work is not persisted. It belongs to the work dependency graph.

### Wide metadata tables

`source_item_metadata`, `analyzed_item_metadata`, and `analyzed_collection_metadata` are wide,
dynamically-widened
typed tables: one row per entity, one bare column per metadata name (`max_current DOUBLE`,
`polarity VARCHAR`). A new name is registered on first write — the flush adds the column with
`ALTER TABLE ADD COLUMN` and a `wide_columns` row in the same transaction. The first type registered
for a name wins; a later mismatched write drops that one value and surfaces a graceful failure naming
the key and types, while everything else in the write proceeds. `missing` is not stored — the key is
simply absent after reload. Symbol/String and Date/DateTime share a SQL type, so the `wide_columns`
registry is the only correct decoder.

### The item-key surrogate

Each logical item has one integer `item_key`, minted at interpretation and shared across the wide
tables and the payload store. `open_cache_db` loads the committed `item_key`s from `items`; a rebuild
resets them. Payloads are keyed by `(item_key, stage)`, where the stage is either the item's own
`process` output or the collection `process` rewrite.

Collection records use a separate package-owned integer key. The index assigns it after resolving
the complete path by deterministic collection ID. That ID combines the parent ID, concrete
collection type, and a canonical encoding of `id(collection)`; it never uses process-dependent
`Base.hash`. The cache persists the key, parent edge, final ID, resolved label and metadata, and
optional registration name. It stores neither user-defined collection values nor canonical ID
inputs. The compact key may change after a clean rebuild, while saved selection and annotation state
reconnect through the deterministic ID.

### Persisted result state

Successful empty analysis produces no metadata rows, so absence of metadata cannot mean "not
computed." `result_states` records item processing and analysis, while
`collection_result_states` records collection processing and analysis using integer collection
keys. Both store success — including empty — or failure. Source-backed typed-item success remains
memory-only.

`result_states` carries no per-result input fingerprint: a result's validity is derived from its
position downstream of unchanged sources, not from a stored claim. At reopen the cache index
restores hierarchy and metadata; missing `RESULT_READY` rows are gap-filled by enqueueing live work,
and the source-metadata diff invalidates whatever changed while the workspace was closed.

These outcomes are independent. Processing failure prevents item analysis for that item, but an
item-analysis failure does not erase a successfully processed payload. Collection-analysis failure
does not alter member item metadata.

### Query surface

`query_items(cache, predicate)` returns item ids whose delivered metadata satisfies a SQL predicate,
over a view `items LEFT JOIN analyzed_item_metadata USING (item_key)`. The view is recreated only
when the effective-metadata column signature changes. Queries read committed DB state, so a value
written within the buffer flush window may not appear yet.

### Processed tabular payloads

Cacheable processed values implementing the Tables.jl interface are stored in native DuckDB tables,
provided every column's element type maps to a DuckDB column type. One physical table is used for
each ordered combination of user column names and element types. The shape's stable digest becomes
an internal storage ID and table name; user text is never used directly as a table name. The
payload's container type is stored with its pointer, and a read rebuilds that container via
`Tables.materializer` — the caller gets back what it stored, whether that was a `DataFrame`, a
`NamedTuple` of vectors, or any other Tables.jl sink. A container type that no longer deserializes
is a cache miss and the work reruns.

Each payload receives a compact integer sequence. Its `item_data` pointer and all physical rows share
that sequence instead of repeating the logical item ID in every row. Two reserved internal columns
store the sequence and the row order. User-derived identifiers are quoted, and the reserved internal
prefix cannot collide with user columns.

The pointer, schema registration when needed, and physical payload rows are committed in one
transaction. A reader cannot observe a committed pointer without its payload. A zero-row payload
still has a pointer and schema, so it reconstructs as an empty value with the correct columns and
types.

Freshness is derived from source data, not from stored per-result claims. When the cache index is
loaded, the result-state tables restore which steps already finished; the work graph holds only live
jobs and gap-fills any missing step. `DirectorySource` refreshes package-owned collection values
from the current `metadata.txt` contents after loading the cached hierarchy. The workspace compares
old and new effective metadata, invalidates only affected item and collection work, and leaves source
files themselves cached. The same reconciliation runs for live metadata-file changes. At runtime
delivery consults live nodes, then cached result state, then enqueues work; payloads are read only
after a row is current.

## Source changes

ProjectCache exposes one semantic cascade for removing a source item:

```julia
delete_source_item!(cache, source_item_id, old_records)
```

It deletes each old record's item rows, wide-metadata rows, failures, result states, disk payloads,
and memory-only values by `item_key`. It does not decide which collections are affected.

The work dependency layer separately identifies affected collection paths and calls the semantic
`delete_collection_metadata!(cache, collection_keys)` operation. A changed source is represented as
buffered deletion of its old subtree followed by buffered stores for the replacement. Per-key
coalescing turns retained keys into replacement edits before any database transaction is chosen.

No maintenance cascade writes directly to DuckDB.

## Cache identity

Each project has one predictable cache file:

```text
DEPOT_PATH[1]/databrowser/<project-name>/cache.duckdb
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
