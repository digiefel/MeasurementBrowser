# Sharded cache writer

Status: design agreed, not yet built. Branch `feature/buffered-cache-scheduler`.
Supersedes the single-locked-writer design currently on that branch. Forward-looking; see
`docs/ARCHITECTURE.md` for current state.

## Why (the problem this fixes)

A profiling trace from the real RuO2 project (`bench/trace.json`, ~5852 source items / ~18185 items,
9 Julia threads) showed processing throughput collapse to **<30–40 items/s** while the scan ran, then
recover after the scan ended. The user has seen ~400 items/s before and expects processing pinned at
100% throughout. Diagnosis from the trace:

- The **single writer connection is saturated**: `writer_busy` climbs to ~30.7s of 41s wall (~70%),
  and **95% of that is interpreted DataFrame payload writes** (`write_item_payloads` = 22.7s of the
  23.8s `reconcile` total; index/metadata/source rows are ~1s combined).
- The flusher barely matches the scan (both ~1300 items in 41s), so the buffer **pins at its
  2,000,000-row ceiling**.
- At the ceiling the backpressure (`_await_capacity!`) parked **every** deposit — including
  processing's background stat deposits, which add **zero** rows. Processing is the *consumer* that
  drains the buffer; parking it there throttled the very stage that relieves the pressure. Each
  `processing/item` span was median 16.7ms / mean 80ms / p95 618ms while the actual `process()` work
  was **1.8ms median, 8s total** — i.e. ~125 thread-seconds spent *parked*, not computing.
- Consequences: the processing queue grew unbounded (455 → 15,660), progress advanced in jumps (the
  scan parks on the same ceiling and drains in flush-sized bursts), and throughput recovered only once
  the scan stopped producing interpreted payloads and the writer drained.

Root cause: **one writer + one shared buffer/ceiling couples four independent streams** (interpreted
writes, processed writes, stats writes, reads) onto one serialized resource, and the shared ceiling
lets the producer (scan) throttle the consumer (processing).

Why the in-repo benchmark missed it: `bench/realistic_browse.jl` payloads are light enough that the
writer never saturates (~35–39% busy), so the coupling never manifests. The benchmark needs a
**write-heavy payload mode** to reproduce this and guard the regression.

## How the sibling branches relate

- **`codex/simplify-cache-reads` == `a513eaa`** (the pre-buffer "fast base", where processing pinned at
  100%): it had **two independent write streams** — interpreted writes ran synchronously in the scan
  task (`flush_cache_batch!` → `reconcile_source_items!`), while processed writes went through a
  *separate* async `processed_write_worker!` with its **own** backpressure
  (`PROCESSED_WRITE_BACKPRESSURE_ROWS = 2_000_000`, coalesce at `PROCESSED_WRITE_BATCH_ROWS =
  1_000_000`). Processing was never gated by the interpreted-write rate. That separation is exactly
  what we lost by unifying everything into one buffer.
- **`integration/best-of-both` (4f22ef3)**: fought this same coupling and never escaped it — commits
  like "Decouple scan interpretation from DuckDB cache writes", "bound scan payload memory", "wait for
  staged interpreted data instead of reconciling during scan", "Fix write-buffer deadlock while the
  cache writer lock is held", and the tip "Fix scan freeze when cache writer backpressure blocks
  cancel." It patched symptoms (deadlocks, freezes, cancel hangs) instead of removing the coupling, got
  spaghettified, and broke at an unimplemented `Base.tryput!`.

The backpressure coupling is the recurring hazard across every integration attempt. This plan removes
it at the root.

## Stopgap already landed

Commit **`23cc9ca`** ("Never block processing on the cache buffer's row ceiling") gates the existing
single buffer's ceiling on the scan producer only; `stage_processed!` no longer calls
`_await_capacity!`. That restores a513eaa's "consumer never blocks" behaviour as a half-measure. It is
**not** sufficient on its own: if analysis (processing) is slower than interpretation, "consumer never
blocks" just moves the failure — the processing queue grows unbounded and the only way to bound it is
to stall the scan, which couples them from the other side. A single shared buffer/ceiling cannot
express independent bounds. This plan does, and supersedes the stopgap.

## The two rules

1. **One writer connection per table.** No table ever has two writers, so DuckDB write-write conflicts
   and retries are impossible *by construction* — not handled, absent. Multiple connections scale here
   because writes to different tables are independent.
2. **Nobody ever waits, except a producer being throttled to keep its own stream from exhausting
   memory.** Producers never write inline — they drop rows into in-memory staging and move on. Readers
   use their own snapshots. The only blocking call in the system is a per-stream memory ceiling, and it
   only ever parks the producer filling that stream.

## The table → writer map

A writer = {in-memory staging queue, one DuckDB connection, one flush loop task}. Because each writer
is the sole task on its connection, there is **no DB write lock** — only a tiny lock around its
in-memory staging queue.

| Writer | Owns | Fed by |
|---|---|---|
| **Payload pool** (N writers) | the big `dataframe_<sha1>` tables, each pinned to one pool member by `hash(storage_id) % N` | scan (interpreted payloads) + processing (processed payloads) |
| **Index** (1) | `item_data` (the pointers) + `dataframe_schemas` | scan + processing, *after* payload commit |
| **Metadata** (1) | `metadata` (item params, item stats, node params/stats) | scan (params) + processing (stats) + analysis (node stats) |
| **Catalog** (1) | `items`, `source_items`, `meta`, `item_failures` | scan + failures |

Every table appears exactly once → no shared tables → no conflicts. The heavy 22.7s lives entirely in
the payload pool and fans out across its connections; the interpreted stream and processed stream land
in different `dataframe_*` tables (processing adds columns → different schema → different table → no
clash), so they write **simultaneously on different connections** instead of fighting one lock.

**Honest limit on fan-out:** parallelism granularity is *per schema* (per payload table), bounded by
how many distinct data shapes exist, not by N. For a few-kind project that's a handful of payload
writers — enough, because the win we need is separating interpreted-payload writing from
processed/stats/reads, and those are always different tables. Parallelizing *within* one schema would
require partitioning a single payload table; not in scope.

## The one cross-writer rule: pointer after target

The only ordering constraint: an `item_data` pointer is written **only after** the payload it points
to is durable.

- Mechanism: a pool writer commits payloads, then hands the resulting
  `(item_id, stage, storage_id, seq, row_count, sif_hash, if_hash, column_names)` tuples to the
  **index** writer, which writes `item_data` + `dataframe_schemas` together. A reader that sees no
  pointer treats it as a cache miss and recomputes — so a torn write is at worst wasted work, never a
  dangling pointer into a half-written payload.
- `seq` (the surrogate tying `item_data` to the payload's `__mb_seq`) is minted by an **atomic
  counter** at deposit time, so payload rows and the pointer carry the same value. This deletes the
  `item_data_seq` DB sequence *and* the `_mint_item_seqs!` round-trip.
- `metadata` and `catalog` carry no pointers into payloads, so their writers flush freely, no ordering.

## Producers (deposit, never wait — except OOM)

- **Scan** interprets a source item → publishes its items to the interpreted read-through store and
  routes rows to: catalog (`items`/`source_items`), metadata (item params), and the payload pool
  (interpreted payloads → then index). Parks **only** when interpreted staging is over its memory
  ceiling.
- **Processing** reads interpreted data (from the read-through store, or payload tables on a miss) →
  produces the processed item → publishes it to the processed read-through store and deposits processed
  payload (pool → index) + stats (metadata). Parks **only** when processed staging is over its own,
  independent ceiling. **Selected/visible items keep priority**: processed first, flushed first, and
  exempt from the *soft* ceiling — they only ever wait at the hard memory limit, the one case the rule
  allows.
- **Failures** → catalog (`item_failures`), never wait.

A slow interpreted-payload writer now throttles **only** the scan (its ceiling); it cannot touch
processing, stats, or reads — different connections. Symmetric: a slow processed writer throttles only
processing. This is the stopgap's intent made structural and two-sided.

## Read-through store (never waits)

Two in-memory stores — interpreted (feeds processing) and processed (feeds plots) — keyed by item id. A
deposit is published immediately, so a just-produced item is readable before it is durable.

**Eviction is uniform and memory-bounded; nothing is pinned indefinitely:**

- **Cacheable item**: evicted once its payload + pointer are committed (durable). A later read falls
  through to DuckDB.
- **Non-cacheable item**: never written to disk, so it cannot fall through — but it is **not** kept
  resident either (that would risk OOM on a stream of uncacheable things). It transits the store like
  everything else and is evicted to bound memory. If an evicted non-cacheable item is needed again, it
  is **re-read from source** via `source_fallback` (re-interpret the source item). Recovery for
  cacheable = disk; for non-cacheable = re-read. Either way the store stays bounded.
- Eviction is tied to consumption / the memory bound, not to a disk commit, so non-cacheable items
  leave memory once they are no longer needed (consumed) rather than lingering. Worst case under
  pressure: an item is evicted before its consumer reads it and gets re-read — acceptable; OOM
  prevention outranks avoiding a re-read.

The plot reader keeps its own multi-threaded ephemeral snapshot connections — already independent,
MVCC, never blocked by any writer.

## Backpressure = memory only

Each producer stream (interpreted, processed) has its own byte/row ceiling. Crossing it parks **only**
that producer, never a consumer, never a reader. Independent ceilings → a slow interpreted writer
throttles only the scan; a slow processed writer throttles only processing. Priority (selected/visible)
items bypass the soft ceiling and only wait at the hard OOM limit.

## What disappears

- The single `writer` connection and `writer_lock` (→ per-table connections, no global write lock).
- Write-conflict retries (impossible by construction).
- The `item_data_seq` DB sequence and `_mint_item_seqs!` probe (→ atomic counter).
- The shared-ceiling coupling (→ independent per-stream ceilings; the `23cc9ca` stopgap becomes
  structural).

## Where the code changes

- **`src/Cache/ProjectCache.jl` (`CacheDB`)**: hold the writer set (payload pool + index + metadata +
  catalog connections) instead of one `writer`/`writer_lock`. Each writer single-tasked → no global
  lock. Schema unchanged except dropping `CREATE SEQUENCE item_data_seq`.
- **`src/Workspace/CacheBuffer.jl`**: the single buffer/`_flush_once!` becomes the two read-through
  stores + the four table-writers, each with its own staging queue, flush loop, and (for producers)
  memory ceiling. The pointer-after-payload hand-off lives here.
- **`src/Cache/ProjectCache.jl` / `src/Cache/ItemDataCache.jl`**: split the write functions by target
  table; `reconcile_source_items!` / `write_cached_item_data!` / `persist_stats!` become "route rows to
  the owning writer"; seq via the atomic counter; payload→index hand-off implements the ordering rule.
- **Read path** (`_read_cached_item_data`, `read_cached_item_data`): essentially unchanged — still joins
  `item_data → dataframe_schemas → dataframe_*` on a committed snapshot.
- **`bench/realistic_browse.jl`**: add a **write-heavy payload mode** (big interpreted DataFrames) that
  saturates the writer the way the real project does, so the stall reproduces and the before/after is
  measurable. This is the verification artifact; the current bench is blind to this failure class.

## Build order

1. **Payload pool + index writer** — the split that ends the 70% saturation and proves the model
   (interpreted vs processed payloads on different connections; pointer-after-payload ordering; seq
   counter). Verify against the new write-heavy bench mode.
2. **Metadata + catalog writers** — fold in stats/params/identity once the model holds.
3. Delete the single `writer_lock` path and the stopgap's now-redundant pieces.

## Open verification

Cannot run the user's real project here. Verify via (a) the new write-heavy bench mode (reproduce the
stall, measure before/after), and (b) the user re-running their project on the branch — processing
should pin and the queue should stay flat with the scan.
