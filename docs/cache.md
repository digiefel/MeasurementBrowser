# Source Cache

Scanning a source is inherently slow: the engine must enumerate every source item, interpret each one
into data items, and then compute analysis values. The tree is known before analysis finishes. The
**cache** stores the result of that work so the previous tree can appear immediately while a new scan
checks the source.
The cache is generated package data, fully rebuildable from the source files. Project code supplies
source interpretation, processing, and statistics; the package decides where and when their reusable
results are stored. Nothing in the cache is a source of truth: its only jobs are to let a reopened
workspace show results quickly and to let processing skip recomputation. Because it is disposable, its
on-disk shape is an implementation detail we are free to change.

## Vocabulary

These words have one meaning each. The rest of this document and the code use them precisely.

- **Buffer** — the high-level API to the database. Code never runs SQL or opens a database connection
  directly; it talks to the buffer. The buffer holds recently touched rows in memory, decides how to
  apply them to the database, serves reads, and writes to disk on its own schedule.
- **Store into the buffer** — hand an object to the buffer (`append`, `edit`, or `delete`). This is the
  only way data enters the cache.
- **Flush** (also **write**) — the buffer draining its held rows to the database. This is the actual
  disk write, and it happens in the background.
- **Backpressure** — a disk-backed buffer making a *store* wait until a flush frees room. The only
  thing that ever blocks. (Memory-only buffers never backpressure — they drop instead.)

There is no "producer" and no "consumer". There is code that sends things to the buffer, and there is
the buffer.

## The buffer model

One buffer sits in front of **each** database table. Each buffer owns its own persistent write
connection and its own persistent read connection, and runs its own background flush task. Because no
two buffers (and no two tasks) ever write the same table, the database can never hit a write–write
conflict — this is a structural guarantee, not something we police at runtime.

The buffer is the single door to the cache. Every read and every change in the rest of the codebase
goes through it.

### Two backing modes (disk-backed and memory-only)

A buffer is *bounded memory in front of something*. That something is one of two things:

- **Disk-backed.** The common case: the buffer sits in front of a DuckDB table, flushes its held rows
  to it on a schedule, and serves reads from memory or that table. This is everything described below
  unless noted.
- **Memory-only.** Some data has **no disk backing at all** — most importantly the
  interpreted-but-unprocessed payloads. There is no point writing them: they are a short-lived
  intermediate that processing consumes to produce the processed data we actually cache, and writing
  them was historically the dominant cache-write cost. A memory-only buffer is a bounded in-memory store
  with **no connection, no table, and no flush task**. It is drained not by a flush but by
  **consumption** (processing taking an item) and by eviction under the memory ceiling. Everything else
  is the same machinery — the memory ceiling, the single read door — but it is **best-effort** and
  **never throttles its producer**: on overflow it simply drops a payload (the item's record is written
  to its own disk-backed buffer regardless), and a later read that misses recomputes from source. So the
  tree fills at interpretation speed no matter how slow processing is. Dropping is safe at any cost ratio
  because the buffer only overflows when interpretation is outrunning processing — exactly when
  re-interpreting is the cheap stage.

This is one mechanism in two configurations. It also unifies the awkward case: a **non-cacheable**
processed item behaves exactly like an interpreted payload — memory-only, never written, recomputed from
source when asked for again. The memory ceiling bounds both, so a stream of either can never leak.

### Storing (potentially blocking)

Storing an object normally returns immediately: it goes into the buffer's in-memory rows and is instantly
visible to reads. A store into a **disk-backed** buffer blocks in exactly one situation — the buffer is
over its memory ceiling — and only until a background flush drains it below the ceiling. The flush never
blocks anything; backpressure is just the buffer holding the store back to keep memory bounded. A store
into a **memory-only** buffer never blocks: on overflow it drops instead (see "Two backing modes").
Reads are never blocked by either, and a flush in progress never blocks a store.

The memory ceiling counts the rows that actually cost memory — the payload data rows. The bookkeeping
buffers (records, metadata, failures) hold tiny rows and in practice never reach their ceiling, so they
effectively never backpressure; the ceiling exists to bound the payload buffers.

### Flushing (background, periodic, non-blocking)

Each disk-backed buffer flushes its held rows to its table when a short interval has elapsed (a couple
of seconds) or when its memory ceiling is hit, whichever comes first — one database transaction per flush, using
DuckDB's bulk Appender. A flush that fails is logged and dropped; the affected data is simply
recomputed from source on the next read, and the flush task keeps running (a dead flush task would
leave a blocked sender stuck forever).

### Reading (memory first, then disk)

A read asks the buffer. If the buffer still holds the rows in memory, it answers from memory with no
disk trip. Otherwise it reads them from the database over its persistent read connection. The caller
cannot tell which happened — that is the point of a single door.

If the buffer has neither the rows in memory nor on disk, the read is a miss. The layer above the
buffer then falls back to the source (re-reading and re-interpreting the original file). So the full
lookup order for a disk-backed item is: **buffer memory → database → source**. For a memory-only buffer
there is no database step, so the order is just **buffer memory → source**. The buffer owns the
in-memory (and, when disk-backed, the database) step; the workspace owns the source fallback.

### The buffer API (append / edit / delete / read)

The buffer exposes the smallest possible verb set (illustrative signatures; the verb names are what
matter):

- `append(buffer, key, object)` — add new rows.
- `edit(buffer, key, object)` — update existing rows.
- `delete(buffer, key)` — remove rows.
- `read(buffer, key)` — read rows back (memory hit, otherwise disk).

Multiple dispatch on the object's type selects the right table and the right way to append/edit/delete
it. The generic buffer machinery never names a table or a column; each table is described by a small
row type plus the dispatch methods for it. Adding a table later is a struct plus a couple of methods.

> **Scope limitation (first implementation).** Only `append` is implemented. The streaming write path
> is **append-only**: during a build we only ever add rows. `edit` and `delete` are the intended API
> but are not built yet, and the operations that genuinely need them are out of scope for now (see
> "Mutations, deferred" below). This is a deliberate scope cut to get a correct, simple buffer landed;
> the code and the API carry comments marking it.

## Stored state

DuckDB stores source-item fingerprints, data-less item records, typed parameters and statistics,
failures, and only **one** item-data stage:

- `processed`: data returned by `process` and consumed by views.

The other stage, `interpreted` (data returned by `data_items` before `process`), is **never written to
DuckDB** — it lives in the memory-only interpreted buffer and is recomputed from source on a miss. Whether
an item has been processed is **derived** from a processed entry existing — it is not stored as a flag.

Cacheable `AbstractDataFrame` payloads live in native DuckDB tables, one physical table per distinct
column shape (column names **and** element types — same names with different types is a different
shape). The first time a shape is seen, its table is created with typed columns; after that, payload
rows are appended row by row with the Appender. Each stored payload gets a compact integer surrogate
(`seq`) from a catalog sequence; the payload rows and the pointer row that locates them share that
`seq` instead of repeating the long item id. The payload rows and their pointer are flushed in the
**same transaction**, so a reader can never see a pointer whose data is not yet on disk. A zero-row
payload still gets a pointer and is reconstructed as an empty, correctly-typed result on read.

A disk-backed payload buffer is also a **read-through cache**: a freshly produced payload is held in
memory so it reads back instantly, and is evicted only **after** its rows are durably committed — and
only if the in-memory copy is still the exact object that was written (a later store may have replaced
it). That ordering is why a reader never reaches a pointer whose payload isn't on disk yet.

Two internal columns ride alongside the user's columns in every measurement table: the `seq` surrogate
and a within-payload row index. They use a reserved `__mb_` prefix that user column names may not take,
so they can never collide. Identifiers derived from data are made SQL-safe before use — a shape's table
name is its content hash rendered as hex (never user text), and user column names are quoted in the
generated `CREATE TABLE`. On read, a cached payload counts as a hit only if its stored fingerprints match
the requested record's; a mismatch is a stale miss, not a hit.

`cacheable(item)` decides whether the **processed** payload is written to disk; interpreted payloads are
never written regardless (they are always memory-only). The built-in `DataItem` DataFrame path is
cacheable; the low-level default is false. A **non-cacheable processed item is therefore handled exactly
like an interpreted one** — held in a memory-only buffer, never written, dropped on overflow or once it
has been consumed, and recomputed from source on a later request. Keeping such items resident
indefinitely would be an unbounded memory leak, so the memory ceiling bounds them like everything else.

### Tables

The fixed base tables, created once at startup before any per-table buffer opens its connection:

- `source_items` — one row per discovered source item: id, fingerprint (hex text), path, timestamp.
- `items` — one data-less record per logical item: id, source-item id, label, kind, collection.
- `metadata` — typed key/value rows (a compact EAV encoding). Item parameters, item statistics, node
  parameters, and node statistics are all `metadata` rows, distinguished by a scope tag.
- `item_failures` — one row per failed item: item id, source-item id, message.
- `item_data` — the **pointer** table: `(item_id, seq, …)` locating a stored payload by its `seq`.

Each base table has its own buffer. The `item_data` pointer table and the per-shape **measurement
tables** (created on first sight of each payload shape) share one owning buffer, so a payload and its
pointer commit together. Interpreted payloads have no table — they are memory-only.

## Reads under load

DuckDB is the only package-level shared cache; there is no Julia object cache layered on top. An active
plot or inspector owns the item objects it currently displays; a later selection reads them back
through the buffer (memory or disk) or repeats the required upstream work. Reads use scalar `seq`
predicates and reconstruct DataFrames as views over the query result columns rather than copying them
again. DuckDB's buffer pool is limited to 1 GiB by default through `CACHE_MEMORY_LIMIT_MIB`;
`set_cache_memory_limit!` changes the default or an open workspace.

## Mutations, deferred

The cold rebuild — the common path, and the one we benchmark — is naturally append-only. The
operations that genuinely change existing rows are:

- **Reprocess** an item that already had processed data (replace its processed payload and statistics).
- **Rescan** a source item whose fingerprint changed (replace its records and dependent results).
- **Clear a failure** when an item that previously failed later succeeds.

These need per-key `edit`/`delete` on the streaming path, which the first implementation does not
provide. Until they exist, **reprocess**, **rescan-of-changed-source**, and **clear-a-failure** are not
incremental yet. They are surfaced as an explicit *work-in-progress* limitation (a quick, honest "not
implemented" marker rather than a silent half-update). The supported way to pick up such a change is a
full **Rebuild Cache** — a wipe-and-append, which is pure append and needs no per-key `edit`/`delete`.

The few bulk/DDL mutations that the cache *does* perform — stamping the identity rows, wiping for a
Rebuild, deleting the source items a scan no longer contains, and the closing checkpoint — still go
through the buffer, never a second writer. Each runs while the per-table buffers are quiescent (the
identity stamp and the deletes are opening moves before any append; the checkpoint runs after the
buffers stop), so the buffer performs them on a transient maintenance connection it opens and closes
for the operation. The single-door guarantee holds: there is no long-lived writer parallel to the
buffers.

Designing per-key reprocess and rescan into the buffer's `edit`/`delete` verbs — so they become
incremental yet fast and batched — is the next piece of work after this buffer lands.

## Location and identity

Each project has one predictable DuckDB file:

```text
DEPOT_PATH[1]/measurementbrowser/<project-name>/cache.duckdb
```

The project name must be one safe path component and is used unchanged. The cache records the project
name and `source_id(source)`. Opening the same project name for a different source fails rather than
creating another cache. Old hash-named caches are ignored.

## Recovery

Corrupt or schema-incompatible generated caches are rebuilt automatically, because no useful user
choice exists. A project/source identity conflict is different: it fails clearly, because silently
replacing another source's valid cache would lose useful state. Failures remain attached to the
smallest source item, logical item, or collection that failed, and do not invalidate unrelated
committed work.

## What we rely on from DuckDB (and what we don't)

DuckDB is mostly safe for concurrent use across connections; we do not assume otherwise. We rely on
only two specific, observed facts:

- **Concurrent `CREATE TABLE` from different connections can conflict on the catalog.** We avoid this
  by giving the measurement-payload tables a single owning buffer that creates them, and by creating
  the fixed base tables once at startup before the per-table connections open.
- **`append_blob` is broken in our DuckDB.jl version.** Fingerprints are therefore stored as hex text
  so the bulk Appender writes them inline.

Everything else (multiple connections reading and writing different tables at once, snapshot reads on a
persistent read connection) is used as normal supported behaviour.

DuckDB's optimistic concurrency model backs this up: **appends never conflict**, even on the same table,
and only two transactions editing the *same row* conflict (the second fails). That is why append-only
streaming is conflict-free regardless of how many sources fan in, and why the later `edit`/`delete`
verbs route every mutation of a table through that table's single owning buffer.

### Handling Concurrency (https://duckdb.org/docs/current/connect/concurrency)

#### Single Process
In in-process mode, DuckDB has two configurable options for concurrency:
Read-write mode: one process can both read and write to the database.
Read-only mode: multiple processes can read from the database, but no processes can write (access_mode = 'READ_ONLY').
When using read-write mode, DuckDB supports multiple writer threads using a combination of MVCC (Multi-Version Concurrency Control) and optimistic concurrency control (see Concurrency within a Single Process), but all within that single writer process. The reason for this concurrency model is to allow for the caching of data in RAM for faster analytical queries, rather than going back and forth to disk during each query. It also allows the caching of function pointers, the database catalog, and other items so that subsequent queries on the same connection are faster.
#### Concurrency Model within a Single Process
DuckDB supports concurrency within a single process according to the following rules. As long as there are no write conflicts, multiple concurrent writes will succeed. Appends will never conflict, even on the same table. Multiple threads can also simultaneously update separate tables or separate subsets of the same table. Optimistic concurrency control comes into play when two threads attempt to edit (update or delete) the same row at the same time. In that situation, the second thread to attempt the edit will fail with a conflict error.
#### Multiple Processes
Writing to DuckDB's native database format from multiple processes is supported through the Quack remote protocol, which turns DuckDB into a client-server database. Quack in beta stage as of DuckDB v1.5.2, and is expected to become mature by DuckDB v2.0 in fall 2026.
For a stable solution, consider using the DuckLake format with PostgreSQL as the catalog database. By coordinating through a central PostgreSQL catalog, DuckDB instances can achieve concurrent read-writes on the same database. The DuckLake v1.0 specification and its DuckDB implementation, both intended for production use, were published in April 2026.

### Optimistic Concurrency Control
DuckDB uses optimistic concurrency control, an approach generally considered to be the best fit for read-intensive analytical database systems as it speeds up read query processing. As a result any transactions that modify the same rows at the same time will cause a transaction conflict error:
```
Transaction conflict: cannot update a table that has been altered!
```
