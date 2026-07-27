# Project Cache Refactor

This design grew out of two apparently small lifecycle problems: accepting a stale-cache prompt
repeats work already completed by the memory fallback, while rebuilding a cache blocks the caller and
recreates more of the session than the operation requires. Fixing either problem locally would add
another transfer path between representations that already overlap and sometimes disagree. The
lifecycle bugs are therefore symptoms of an unresolved ownership boundary rather than isolated GUI
or DuckDB problems.

The goal is a top-down cache model that works identically for GUI, REPL, script, and future CLI
callers. Disk persistence must remain optional; a memory-only cache is a supported storage policy,
not an error state or temporary approximation. Changing policies must preserve completed work when
the destination can represent it, while rebuilding must intentionally discard completed work and
restart asynchronously. Neither operation should reopen the workspace, rewarm project extensions,
or expose backend-specific sequencing to callers.

This plan assumes the typed-API consolidation is complete. `read`, `entries`, item and collection
`process`, and item and collection `analyze` are already the engine's canonical stages;
`cacheable(stage, value)` and `construct` already define persistence and rehydration at those
boundaries. This refactor consumes that vocabulary rather than designing an interim cache around the
current monolithic interpretation path. In that sense it is also a conformance test: the cache must
be implementable from the persistence-neutral API contracts without depending on recipe carriers,
Core work-node types, or concrete workspace internals.

## What we actually have

`AbstractCacheDB` is barely an abstraction:

- Core imports roughly forty cache operations.
- `cache_work_status` branches on `isa CacheDB` and directly accesses `.lock`, `.failures`,
  `.source_items`, and result-state dictionaries.
- Other Core code directly reads `.stage_ledger`.
- `MemoryCacheDB` has several deliberately incomplete or no-op implementations.
- `CacheDB` is already hybrid: DuckDB stores plus three bounded `MemoryStore`s.

Completed state currently exists in several overlapping forms:

- `WorkspaceIndex`
- `ProjectCacheIndex`
- `CacheStageLedger`
- DuckDB row stores
- partial dictionaries in `MemoryCacheDB`

Core then mutates the workspace and cache separately. Sometimes payloads are stored before
publication; sometimes the index is changed before cache bookkeeping. The eviction/`RESULT_READY`
bug is one consequence of that divided ownership.

The roadmap already names this unresolved problem: clarify ownership between workspace, index,
project cache, database, and write buffers; remove duplicated state and layer-skipping call paths.

## The better seam

A cache should mean the complete collection of finished pipeline state, independent of its physical
storage:

- source fingerprints and entry records;
- hierarchy and metadata;
- result status and failures;
- available stage payloads.

The workspace should own:

- the open source;
- the live work graph;
- scheduling and cancellation;
- selection;
- the active cache.

Its visible index should be the cache's lightweight in-memory index, not another separately
maintained authority. `workspace.index` can initially remain as an accessor or reference to that
object so existing headless and GUI consumers do not need to change simultaneously.

Queued and running work remains exclusively in Core. Completed work moves atomically into the cache.
This is materially different from today even though the package dependency diagram remains similar.

## The abstraction should be compositional

The current "disk cache versus memory cache" distinction is misleading. `CacheDB` already mixes
disk-backed and bounded-memory stores.

A better logical shape is:

```text
ProjectCache
├── identity
├── lightweight index
└── results by pipeline stage
    ├── read
    ├── entries
    ├── item process/analyze
    └── collection process/analyze
```

Each result store can use an appropriate policy: durable DuckDB, bounded memory, or a future
specialized backend. The consolidated typed pipeline gives every expensive stage its own retention
and persistence boundary; this refactor gives those boundaries one coherent owner.

Memory mode means:

- a complete lightweight index and ledger in RAM;
- bounded payload stores;
- no durable backing.

Disk mode means:

- the same logical cache;
- durable backing where the destination supports it;
- bounded memory where it does not.

Different cache shapes then become natural compositions instead of implementations forced around
today's three resident queues.

## Construction and conversion

Backend conversion should use ordinary concrete destination construction:

```julia
DuckDBCache(old::AbstractProjectCache)
MemoryCache(old::AbstractProjectCache)
```

The destination owns conversion policy. It may use a generic logical copy, specialize for a concrete
source, adopt compatible stores, stream between representations, or omit values it cannot retain.
There is no public operation tied specifically to today's concept of "resident results." Custom
cache implementations define their own destination constructors and remain free to use different
physical shapes.

Every current cache already carries `ProjectCacheIdentity`, including `MemoryCacheDB`, so a
conversion can inherit logical ownership from its source. In the final model, logical identity and
backend configuration should be distinct: identity binds completed work to a project and source,
while a DuckDB path, capacity limits, and similar choices configure one physical implementation.
That separation need not block the ownership refactor.

Conversion is a paused handover:

```text
pause producers
    → construct destination from active cache
    → atomically replace workspace.cache
    → resume producers
```

The source cache remains intact until construction succeeds. A failed conversion closes and removes
the incomplete destination, keeps the old cache active, and resumes work. A stale-cache acceptance
uses this conversion path. It does not clear the tree or scan the source.

Rebuild has intentionally different semantics: pause producers, replace the active cache with an
empty cache using the same configured policy, clear the visible completed index, and immediately
start a fresh asynchronous scan. The operation begins promptly; repopulation remains real background
work.

## Package responsibilities

The forward-looking split is:

- **DataBrowserAPI:** the logical cache contract and persistence-neutral index/result types, without
  DuckDB.
- **DataBrowserCache:** built-in memory and DuckDB implementations, result stores, buffering, and
  conversion constructors.
- **DataBrowserCore:** workspace execution, the work graph, publication, pause/resume, and atomic
  cache replacement.
- **GUI, REPL, scripts, and future CLI:** invoke the same Core lifecycle operations and never inspect
  a concrete cache implementation.

This differs from today in the seam, not the package names. `AbstractCacheDB` currently lives beside
the DuckDB implementation, Core knows both concrete cache layouts, and GUI imports cache-specific
types for lifecycle decisions. Moving the logical contract to API and removing concrete field access
makes cache implementations substitutable in practice rather than only in their declared subtype.

## Scope and migration

The consolidated typed API is a prerequisite, not part of this refactor. This plan starts from its
final stage, cacheability, metadata-layer, and rehydration contracts and changes ownership and
storage beneath them:

1. Consolidate `WorkspaceIndex` and `ProjectCacheIndex` into one lightweight cache-owned index.
2. Make memory mode retain all lightweight cache state faithfully.
3. Replace Core's concrete cache field access with the logical contract.
4. Add destination constructors between the built-in cache implementations.
5. Implement pause, destination construction, atomic swap, and resume.
6. Define rebuild as replacement with a fresh empty cache of the same configuration.
7. Implement physical result stores against the established typed-stage contracts without exposing
   the current queue layout as public API.
8. Move stale-cache and rebuild UI behavior onto the corresponding Core lifecycle operations so the
   render task never performs cache work.

This is larger than fixing two buttons, but it removes existing duplication and establishes the seam
required to apply the typed API's per-stage persistence consistently. A workspace-specific replay
function would preserve the current divided ownership and duplicate stage semantics that already
belong to the API.

## Target vocabulary and architecture

The intended end state uses these names and responsibilities:

- **`AbstractProjectCache`** — the persistence-neutral contract for completed project work.
- **`ProjectIndex`** — the single lightweight, in-memory index of completed records, hierarchy,
  metadata, failures, fingerprints, and result availability. It replaces the overlapping
  `WorkspaceIndex`, `ProjectCacheIndex`, and externally inspected portions of `CacheStageLedger`.
- **`MemoryCache`** — a volatile `AbstractProjectCache` with a complete `ProjectIndex` and bounded
  payload stores.
- **`DuckDBCache`** — a durable `AbstractProjectCache` with the same logical index and stage results,
  backed by DuckDB where supported and bounded memory where necessary.
- **`CacheIdentity`** — the project/source identity of completed work.
- **cache configuration** — backend-specific location and resource policy, separate from identity.
- **`Workspace`** — the headless execution owner: source, selection, live work graph, cancellation,
  and one active `AbstractProjectCache`. It does not own a second completed-state index.

The runtime ownership diagram is:

```text
source changes
      │
      ▼
DataBrowserCore.Workspace
├── WorkDependencyGraph       queued and running work only
├── selection/cancellation    session control
└── cache::AbstractProjectCache
    ├── ProjectIndex          all lightweight completed state
    └── stage result stores   available payloads under backend policy
```

The package diagram is:

```text
DataBrowserAPI
  AbstractProjectCache, ProjectIndex, CacheIdentity, result contracts
        ▲
        │ implements
DataBrowserCache
  MemoryCache, DuckDBCache, buffers and physical stores
        ▲
        │ consumed through the API contract
DataBrowserCore
  Workspace, WorkDependencyGraph, lifecycle and atomic handover
        ▲
        │ observed/controlled
GUI / REPL / scripts / CLI
```

The final architecture has one invariant from which the lifecycle behavior follows: every completed
result has exactly one logical owner, the active cache. A result reported as ready must be retrievable
from that cache; eviction removes readiness atomically. The work graph contains only unfinished work,
observers read only the cache-owned `ProjectIndex`, constructors convert complete cache values
without involving the source, and rebuild is the explicit operation that replaces that value with an
empty one. Implementation planning must therefore trace each current field and mutation in
`WorkspaceIndex`, `ProjectCacheIndex`, `CacheStageLedger`, `MemoryCacheDB`, and `CacheDB` into this
single ownership model. The resulting cache uses the consolidated API's stage identities, result
values, `cacheable` policy, metadata layers, and `construct` contract as its complete domain context;
source discovery, live work, selection, and frontend state remain outside it.
