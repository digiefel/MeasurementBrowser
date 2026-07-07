# Live job graph + cache authority

**Goal:** shrink the in-memory work graph to the *live frontier* and let the smarter graph — not
a side index and not the cache — answer every scheduling question. Fewer lines, less cognitive
load, O(1) scheduling, and per-item memory bounded by work in flight, not work ever done.

---

## 1. Problem

`workspace.work.nodes` keeps every processing step for the session (~2 per item), almost all
`:ready`. Three costs follow:

1. **Memory** — a fat `WorkNode` (dependency vector + waiters + reverse-edge entries) survives per
   finished step, plus a global `dependents` dict.
2. **Status latency** — the status bar recomputes "what's running?" by scanning the whole map
   (`work_counts`), so cost grows with items *scanned*, not items *processing*.
3. **Scheduling cost** — `wake_ready_dependents!` re-evaluates a dependent's *entire* dependency
   set on every dependency completion. An N-member collection drains in O(N²) readiness checks.

**Fix:** `nodes` holds only `:waiting` / `:queued` / `:running` jobs. Finished work lives in the
cache. Each node tracks *how many* of its dependencies remain (`pending`) and *who* is waiting on
it (`dependents`); it does **not** remember which dependencies those were. Completion decrements;
readiness is `pending == 0`. No global dict, no counter side-index, no cache reads on the hot path.

---

## 2. Glossary

| Term | Meaning |
|------|---------|
| **Job / node** | One unit of background work + its queue entry. `WorkNode`, state `:waiting`, `:queued`, or `:running` only. |
| **Job key** | `WorkKey(kind, entity)`. Kinds: `SOURCE_INTERPRET`, `ITEM_PROCESS`, `ITEM_ANALYZE`, `COLLECTION_PROCESS`, `COLLECTION_ANALYZE`. Stable across revisions. |
| **Dependencies** | What a node is waiting on (member analyses, an item process, …). Tracked only as a **count** (`pending`). |
| **Dependents** | Who is waiting on a node. Stored on the node as `dependents :: Set{WorkKey}`. |
| **Job graph** | `workspace.work.nodes` — only waiting/queued/running jobs. |
| **Cache** | DuckDB or memory store: payloads + `result_states` (`RESULT_READY` / `RESULT_FAILED`). Authoritative for finished work. |
| **Delivery** | Serving processed data to a plot/table (`request_processed_items`). |
| **Invalidation** | Inputs changed → affected work reruns (`invalidate_records_work!`). |
| **Publish** | A worker finished → write cache, update index, remove node, decrement dependents. |

---

## 3. Design rules (non-negotiable)

1. **Successfully finished** ⇒ `result_states` row `RESULT_READY` (+ payload). **Not** in the graph.
2. **Failed** ⇒ `result_states` row `RESULT_FAILED` (+ `analysis_errors`). **Not** in the graph.
3. **Live (any of waiting/queued/running)** ⇒ in the graph, and *only* these are in the graph.
4. **Unmet-dependency invariant:** if a job has an unmet dependency, that dependency has a live
   node. Equivalently: *no node for a key ⇔ that key is finished (ready or failed)*. This is what
   the scheduler is allowed to assume, and what makes `pending` seedable from `haskey(nodes, dep)`
   with zero cache reads. It holds because dependencies are always enqueued **before** the
   dependents that need them (creation order), and because invalidation (rule 6) re-enqueues the
   whole dependent subtree.
5. **Delivery** ⇒ one path: node present → wait on it; else cache `READY` → read; `FAILED` →
   fail; absent → enqueue and wait.
6. **Invalidation** ⇒ enqueue the entire affected **dependent** subtree immediately (item process,
   item analyze, and every ancestor collection process/analyze). Jobs whose dependencies aren't
   ready sit `:waiting`. Never wait for a parent to publish before scheduling children; never
   delete cache payloads as the mechanism.
7. **Status** ⇒ active count is an O(live) scan of `nodes` (small by construction). No side counter.

---

## 4. The scheduling core

Two facts per node, nothing else:

```julia
mutable struct WorkNode
    key        :: WorkKey
    revision   :: UInt16            # (no realistic need for more than 16 bits)
    state      :: Symbol            # :waiting | :queued | :running
    priority   :: Int
    dependents :: Set{WorkKey}      # who is waiting on me
    pending    :: UInt64            # how many of my dependencies are still unmet
    waiters    :: Vector{Channel{Any}}
    queued_ns  :: UInt64
end
```

**Removed:** `dependencies :: Vector{WorkKey}` (a node no longer remembers what it waits on),
and the graph's `dependents :: Dict{WorkKey,Set{WorkKey}}` (each node owns its own).

**Create a node** (under `graph.lock`): for each requested dependency `dep`, if
`haskey(nodes, dep)` it is live and unmet — push my key into `dep.dependents` and increment my
`pending`; otherwise it is finished — ignore it. Then `pending == 0` ⇒ queue immediately.

**Publish A** (under `graph.lock`): walk `A.dependents`, decrement each live dependent's `pending`,
queue any that hit zero and are still `:waiting`, then remove A's node.

```
publish A:
    for D_key in A.dependents
        D = get(nodes, D_key, nothing)
        D === nothing && continue            # removed since; skip
        D.pending -= 1
        D.pending == 0 && D.state === :waiting && queue!(D)
    delete!(nodes, A.key)                     # A.dependents goes with it
```

**Re-enqueue reseeds; it never patches.** Invalidating a key replaces its node with a fresh one
(revision bumped, `pending` reset to 0) and re-seeds `pending`/`dependents` from the dependency list
the caller supplies *again* — the node never stored it. Replacing the value in `nodes` gives correct
*identity resolution*, but a fresh node has an **empty** `dependents` set, so continuity does not
come from the swap: it comes from rule 6 re-enqueueing the *whole* dependent subtree in dependency
order. A replaced dependency `A'` and every dependent `B` that waited on it are rebuilt in the same
pass — `B` re-seeds against the live `A'` and re-registers into `A'.dependents`; the old node's set
is discarded with it. Because `dependents` is keyed by `WorkKey`, re-registering an already-listed
dependent is idempotent, so a still-live dependency decrements each dependent exactly once. This is
why re-enqueue must **reset** `pending`, not add to it, and why no node needs its predecessor's edges
carried forward.

**Failure needs no special scheduler branch.** Publish decrements dependents whether A succeeded or
failed; the node vanishes either way. A dependent that later runs on a failed input reads the
missing/`FAILED` payload from cache and fails itself — same net result the per-kind
`dependency_ready_for` used to produce, minus the branching. (Collections already tolerate member
failure by folding the error; that behavior is unchanged and now falls out for free.) **This
deletes `dependency_ready_for` and the readiness re-scan in `dependencies_ready`.**

**Why no ready-index and no cache read here:** the only scheduling question is "is dependency
`dep` still unmet?", and rule 4 makes that exactly `haskey(nodes, dep)` — an in-memory O(1) check.
The cache is consulted only at *delivery*, which is user-paced.

---

## 5. Delivery: `ensure_uptodate!`

Single entry for every delivery path (`Processing.jl`):

```
ensure_uptodate!(workspace, key; dependencies, priority):
    node present for key?  → attach waiter, wait, done
    else cache_work_status(workspace, key):
        READY   → done (caller reads payload)
        FAILED  → surface failure
        ABSENT  → enqueue_work!(key; dependencies, priority); attach waiter; wait
```

`request_processed_items` becomes: per record, `delivery_gate` → `ensure_uptodate!` → read payload
via existing `_delivered_payload` / `read_item_data`. The `work_ready`-node branch is gone. The
collection gate is not a separate function: the collection branch promotes its member `ITEM_PROCESS`
jobs and calls `ensure_uptodate!` on the `COLLECTION_PROCESS` key with `dependencies` = the member
`ITEM_ANALYZE` keys — the same `enqueue_work!` path everything else uses. `enqueue_collection_delivery!`
is deleted, not wrapped.

---

## 6. Structural changes

### `WorkDependencyGraph` (`Workspace.jl`)
- **Remove** `dependents` dict. **Keep** `nodes`, `queue`, `workers`, locks, `total`, `completed`,
  `source_*`, `closed`.
- `nodes` holds live jobs only.

### `WorkNode` (`Workspace.jl`)
- **Add** `dependents :: Set{WorkKey}`, `pending :: Int`.
- **Remove** `dependencies :: Vector{WorkKey}`.
- **Keep** `revision` on the live node only (stale-publish rejection). Revision is not persisted
  across completions: when no live node exists, the next enqueue starts at 1. Stale-worker guards
  (`work_node_current`, `current_completion_node`) make this safe.

### `Processing.jl`
| Action | Function | Note |
|--------|----------|------|
| Add | `ensure_uptodate!` | §5 |
| Add | `cache_work_status(workspace, key)` → `:ready` / `:failed` / `:absent` | maps `WorkKey.kind` → cache kind (§9) |
| Rewrite | `enqueue_work!` | Create-node path seeds `pending`/`dependents` per §4; re-enqueue resets `pending` to 0; no `:ready`/`:missing` states; terminal short-circuit answers from cache, not from a node |
| Rewrite | `wake_ready_dependents!` → decrement loop | §4; no readiness re-scan |
| Rewrite | `dependencies_ready` | `node.pending == 0` |
| Rewrite | `work_counts` | `(completed, total, count(queued/running))` over the now-small `nodes` |
| Rewrite | `work_kind_running` | scan of `nodes` (small) |
| Rewrite | `finish_work_node!` → publish path | `delete!(nodes, key)` + decrement dependents; no `state = :ready` |
| Rewrite | `cancel_waiting_work!` | delete the waiting/queued nodes and fail their waiters; no `:missing`, no `pending` fixup (lazy skip handles orphaned edges; orphaned `pending` on survivors is acceptable — profile restart calls `reset_work_graph!` on idle) |
| Delete | `dependency_ready_for` | subsumed by "always decrement" |
| Delete | `work_ready` | callers use `ensure_uptodate!` / `cache_work_status` |
| Delete | `bump_work_missing!` | invalidation enqueues instead |
| Delete | `seed_work_node!`, `set_dependencies!`, `remove_dependency_edges!`, `replace_work_node!`'s edge bookkeeping | no forward list, no eager teardown (lazy skip on publish) |

### `Operations.jl`
| Action | Function | Note |
|--------|----------|------|
| Add | `enqueue_dependent_subtree!` | The one fan-out helper: given records/collection keys, enqueue `ITEM_PROCESS`, `ITEM_ANALYZE`, and every affected `COLLECTION_PROCESS` / `COLLECTION_ANALYZE` in dependency order. Every former `mark_collection_work_missing!` / `bump_work_missing!` caller routes here. |
| Delete | `seed_cached_work_nodes!` | reopen no longer hydrates finished nodes |
| Rewrite | `apply_cache_index!` | gap-fill: for each item lacking `RESULT_READY` process/analyze, `enqueue_processing!` / `enqueue_item_analysis!` (no node seeding); failure restore reads `cache_work_status(key) === :failed` → set `analysis_errors`, not a `:failed` node |
| Rewrite | `invalidate_records_work!` | call `enqueue_dependent_subtree!` for the changed records; keep metadata/index cleanup; drop `bump_work_missing!` |
| Delete | `mark_collection_work_missing!` | its callers (metadata invalidation, source removal in `ingest_source_changes!`, `SOURCE_INTERPRET` failure) call `enqueue_dependent_subtree!` |
| Rewrite | `publish_work_success!` / `publish_work_completion!` | chain the steady-state next step (`enqueue_item_analysis!` / collection work) **while the node is still present**, then the publish path deletes the node — the delete's decrement wakes the chained job; keep `work_node_current` + revision guards so a superseded worker neither writes cache nor decrements |
| Rewrite | `current_completion_node` | node still `:running`, same object |

### `Cache/ProjectCache.jl`
- **No new query APIs.** `cache_work_status` lives in `Processing.jl` (~25 lines): map
  `WorkKey.kind` → ledger row, read `result_states` (or interpret tables for `SOURCE_INTERPRET`).
- **MemoryCacheDB parity only:** track `result_states`, `source_items`, and `failures` on
  store/delete so the in-memory cache matches `CacheDB` ledger writes.
- **Keep** `store_processed!`, `store_item_metadata!`, `store_result_failure!`.
- **Do not** add payload deletes to invalidation.

### `Status.jl`
- Unchanged; benefits automatically from the fixed `work_counts`.

---

## 7. Phasing

The node-deletion, cache-delivery, and counter-scheduling changes are **one atomic landing** — you
cannot delete `:ready` nodes without cache-authoritative delivery and counter dependencies at the
same time, and tests won't be green in between. Sequence the edits, but treat 1–3 as a single PR.

1. **Scheduling core.** Swap `WorkNode` to `dependents`/`pending`; rewrite create + publish +
   `dependencies_ready`; delete `dependency_ready_for`, `set_dependencies!`,
   `remove_dependency_edges!`, forward list, global dict.
2. **Delivery.** Add `cache_work_status` + `ensure_uptodate!`; rewrite `request_processed_items`
   (collection branch inline, `enqueue_collection_delivery!` deleted) and `enqueue_work!`'s terminal
   short-circuit to the cache; delete `work_ready`.
3. **Publish deletes.** `finish_work_node!` removes the node; chain-then-delete ordering in the
   publish path; fix `current_completion_node`; rewrite `cancel_waiting_work!`; `work_counts` →
   O(live) scan.
4. **Invalidation.** Add `enqueue_dependent_subtree!`; route `invalidate_records_work!`, source
   removal, and `SOURCE_INTERPRET` failure through it; delete `bump_work_missing!` and
   `mark_collection_work_missing!`.
5. **Cache load.** Delete `seed_cached_work_nodes!` / `seed_work_node!`; rewrite `apply_cache_index!`
   gap-fill and failure restore.
6. **MemoryCache.** Confirm `result_states` + payload reads work identically without a persistent DB.
7. **Cleanup + review.** Remove any `:ready`/`:missing`/`:failed` node-state remnants and
   `bump_revision!` if unused; drop `dependents` dict cleanup from `reset_work_graph!`; code
   review; update `docs/ARCHITECTURE.md` work-graph section and fold this into it.

---

## 8. Tests

- Update `test/test_work_graph.jl`: dependencies via `pending`/`dependents`, not persistent nodes.
- Invalidate one item's metadata → full dependent subtree enqueued → jobs wait until `pending`
  clears → correct data after drain.
- **Re-invalidation while waiting:** invalidate `A` while a dependent `B` is `:waiting` → both
  re-enqueued → `B` re-registers into the new `A'` → `B.pending` reaches 0 → runs. (Guards point 2.)
- Collection analyze runs only after member analyses complete (counter reaches zero).
- **Failure propagation:** `process` fails on one member → dependents still decrement → collection
  proceeds with the fold error; item-kind dependent surfaces the failure. Net metadata unchanged
  from today.
- **Failed → fixed → rerun:** re-register recipe / re-edit invalidates the `FAILED` key →
  re-enqueue → success clears the `FAILED` row.
- Reopen: no node hydration; only missing steps enqueued; no O(items) allocation.
- Stale worker: invalidate while a job runs (new revision, same key) → old worker's completion
  neither writes cache nor decrements dependents.
- **Scaling guard:** `bench/scaling.jl` already times `status_refresh` (`workspace_status`) and fits
  its exponent — confirm it stays flat (~0) with the live-only graph. No new test file.

---

## 9. `WorkKey` → cache mapping

| `WorkKey.kind` | `result_states` kind | Payload stage |
|----------------|----------------------|---------------|
| `ITEM_PROCESS` | `PROCESSING_RESULT` | `:processed` |
| `ITEM_ANALYZE` | `ITEM_ANALYSIS_RESULT` | metadata tables |
| `COLLECTION_PROCESS` | `COLLECTION_PROCESS_RESULT` | `:collection_processed` / `:processed` |
| `COLLECTION_ANALYZE` | `COLLECTION_ANALYSIS_RESULT` | collection metadata |
| `SOURCE_INTERPRET` | *(interpreted-items table; not in `result_states`)* | `:interpreted` |

`SOURCE_INTERPRET` is the one non-uniform status: `cache_work_status` special-cases it
(interpreted-items presence; interpret failure surfaced through the same failure record path).

---

## 10. Done criteria

- [ ] `work_counts` / status refresh O(live), not O(items scanned).
- [ ] `WorkNode` has no `dependencies` field; graph has no `dependents` dict; no ready-index.
- [ ] Scheduling touches only in-memory `nodes`; cache is read only at delivery.
- [ ] Plot on processed selection: reads cache, no reprocess.
- [ ] Metadata edit: full dependent subtree enqueued immediately; plot gets new data after drain.
- [ ] Collection analyze runs after member analyses complete.
- [ ] Reopen: no `seed_cached_work_nodes!`; only missing steps enqueued; failure restore via `cache_work_status`.
- [ ] Re-invalidating a dependency while a dependent waits leaves `pending` correct (reaches 0).
- [ ] Superseded worker completion updates neither cache nor graph.
- [ ] One fan-out helper (`enqueue_dependent_subtree!`) serves metadata, removal, and interpret-failure paths.
- [ ] `dependency_ready_for`, `work_ready`, `bump_work_missing!`, `mark_collection_work_missing!`, `seed_*_work_node*`, `enqueue_collection_delivery!` deleted.
- [ ] Net smaller codebase, zero dead branches.
- [ ] `bench/scaling.jl` `status_refresh` exponent stays flat.
- [ ] `Pkg.test()` green.
