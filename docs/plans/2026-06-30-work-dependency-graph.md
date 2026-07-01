# Work dependencies publish into WorkspaceIndex

## Boundary

`WorkspaceIndex` remains the sole owner of completed, published records, hierarchy, parameters,
statistics, and failures. `WorkDependencyGraph` replaces the scheduling responsibilities currently
split across `ProcessingQueue`, the readiness probe, `Workspace.analysis`, `start_analysis!`, and
scan-time processing enqueue logic.

```text
source item
    └─ interpretation
         ├─ item records ─────────────────────────────→ WorkspaceIndex
         └─ interpreted item data (memory-only)
                    ↓
               processing ────────────────────────────→ ProjectCache
                    ↓
               item statistics ───────────────────────→ WorkspaceIndex + ProjectCache
                    ↓ much later
               collection statistics ────────────────→ WorkspaceIndex + ProjectCache
```

The graph owns work identities, dependencies, revisions, state, priority, running jobs, and selected
waiters. It references source, item, and collection identities but does not duplicate completed data
from `WorkspaceIndex`.

## Work nodes

Use four explicit work kinds, not a generic graph framework:

- source interpretation;
- item processing;
- item statistics;
- collection statistics.

Each node has an entity identity, input revision, state (`missing`, `queued`, `running`, `ready`, or
`failed`), priority, dependencies, and any selected waiters. Item jobs capture the current
`ItemRecord` revision. Collection membership or member-stat changes increment the collection
revision. Completion for a revision that is no longer current is discarded.

Persisted ProjectCache state seeds ready and failed work when the workspace opens. A cached downstream
result removes any need to reconstruct its upstream inputs. Queued and running state remains
memory-only.

## Source changes

```julia
struct SourceChanges{S<:AbstractDataSourceItem}
    upserts::Vector{S}
    removals::Vector{String}
end
```

`apply_source_changes!` is the only source-mutation entry point. A rescan discovers a snapshot,
compares IDs and fingerprints with the current source snapshot, and submits one `SourceChanges`
batch. A future watcher submits singleton or batched changes through the same function.

- An unchanged fingerprint creates no work and no index mutation.
- A new or changed source item creates one interpretation node.
- A removed source item immediately removes its published outputs and schedules ProjectCache deletes.
- A source-wide collection-parameter change upserts every affected source item.
- Interpretation builds and validates its replacement privately.
- While interpretation runs, the prior published result may remain visible.
- Successful interpretation atomically replaces that source item's published records.
- Failed interpretation removes the stale records and publishes the source failure.

## WorkspaceIndex publication

Add per-source mutation operations to `WorkspaceIndex`:

- validate replacement IDs against unaffected records before mutation;
- remove or atomically replace all records owned by one source item;
- update source fingerprints and failures;
- remove stale item statistics and errors;
- insert replacement records;
- prune empty hierarchy nodes;
- invalidate affected collection statistics and ancestors;
- preserve unrelated records, hierarchy branches, statistics, and selections.

A global collection-parameter change may reapply parameters across the whole hierarchy.

Workers never mutate `WorkspaceIndex`. They emit typed completion events carrying work kind, entity,
revision, and result or failure. `poll_workspace!` rejects stale completion, updates the work graph,
publishes into the index, sends persistence operations to ProjectCache, unlocks dependents, and
recomputes `WorkspaceStatus`.

## Scheduling

Priority order:

1. processing required by the current selection;
2. source interpretation caused by source changes;
3. background processing and item statistics;
4. collection statistics.

Permanent background demand requests processing and item statistics for every current item.

- One result revision has at most one job.
- Persisted ready results satisfy demand without jobs.
- Selection promotes and joins existing processing work.
- Fresh interpreted data feeds processing directly.
- When processing is missing after reopen, one source-interpretation job reconstructs interpreted
  values for the required items from that source item.
- Processing publishes and stores its result independently.
- Fresh processing passes its resident output directly into item statistics.
- Cached processed data with missing item statistics creates a stats-only job.
- Processing failure blocks item statistics for that revision.
- Item-stat failure is independent and terminal.
- Collection work becomes eligible only after the current source-change batch settles and all member
  item-stat nodes are terminal.
- Only affected collections and ancestors are invalidated.
- Failed revisions are retried only after a source change or explicit rebuild.
- An evicted memory-only processed result becomes missing when requested and is recomputed.

`WorkspaceStatus` remains the only watcher-facing summary and derives progress from the work graph
rather than separate processing and analysis schedulers.

## Cache contract

Workspace uses only semantic Cache operations for lifecycle, cache loading, storing interpreted and
processed data, storing independent statistics/results/failures, reading processed data, deleting a
source subtree, invalidating collection results, and reporting pending cache work. It never accesses
`CacheBuffer`, `CacheDB` buffer fields, SQL, connections, tables, or columns.

ProjectCache persists successful empty item/collection statistics and independent result failures so
the work graph never guesses after reopen. Processed success remains represented by its payload
pointer; non-cacheable processed success remains session-only. Buffer mechanics and connection
ownership belong exclusively to the companion cache-buffer plan.

## Verification

- cold build computes every result exactly once;
- unchanged reopen performs no project callback and reads no payload body;
- selection joins and promotes existing background work;
- add/change/remove affects only the corresponding descendants;
- source change during running work rejects stale completion;
- empty statistics remain complete after reopen;
- processing, item-stat, and collection-stat failures are independent;
- collection work waits for terminal item-stat dependencies;
- per-source replacement preserves unrelated index output and selections;
- rescan with no changes creates no work.

The current-state architecture, source/cache/processing, and streaming-status documentation must be
updated with the implementation.
