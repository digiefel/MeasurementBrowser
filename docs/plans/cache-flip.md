# Cache Edge Flip — Execution Plan

**Status (2026-07-11): executed.** Kept for reference; current behavior is documented in
[../cache.md](../cache.md) and [../data-model.md](../data-model.md).

Flip the residual dependency edge from the package split: today `DataBrowserCache` is a re-export
shell over `DataBrowserCore.Cache`, and Core's [Cache.jl](../../lib/DataBrowserCore/src/Cache.jl)
reaches across packages with a `joinpath(@__DIR__, "..", "..", "DataBrowserCache", "src")` include
hack. The target is `DataBrowserCore -> DataBrowserCache`, with Cache a real package depending only
on `DataBrowserAPI`, `DataBrowserProfiling`, and its storage backend (DuckDB/DBInterface/DataFrames).

Written to be executed phase by phase; **the full test suite must pass and the app must be
committable after every phase**. Run tests with
`julia --project --threads=4 -e 'using Pkg; Pkg.test()'`.

Out of scope (do not touch): the diagnostics/bench overhaul, any generic-plotter work, changing what
the cache stores or how (the DuckDB schema, buffers, and write paths move verbatim), and the
`MB_BENCH_ENGINE_ONLY` cleanup (that lands when `bench/` moves to Core alone, per the roadmap).

## Current state (verified 2026-07-11)

- The cache implementation already lives in `lib/DataBrowserCache/src/` (`build_metrics.jl`,
  `cache_buffer.jl`, `project_cache_domain.jl`, ~3.2k lines) but is `include`d by
  `lib/DataBrowserCore/src/Cache.jl` (47 lines: module header + imports + cross-package includes).
  `lib/DataBrowserCache/src/DataBrowserCache.jl` is a shell that does `using DataBrowserCore` and
  re-exports `DataBrowserCore.Cache` names.
- `build_metrics.jl` and `cache_buffer.jl` have **zero** Core coupling.
- `project_cache_domain.jl` is the only file importing Core symbols — 11 names from `..ItemIndex`:
  `DataItem`, `Hierarchy`, `ItemFailure`, `ItemRecord`, `MetadataValue`, `SourceScan`, `all_items`,
  `collection_path_tuple`, `emit_progress`, `insert_item!`, `metadata_dict`. The cache rebuilds the
  item index on warm reopen, so it needs the real types.
- `ItemIndex` (`lib/DataBrowserCore/src/ItemIndex.jl`, 713 lines, plus `ItemIndex/Scanning.jl`,
  93 lines) depends only on `DataBrowserAPI`, `DataBrowserProfiling`, `Dates`, and one line of
  `DataFrames` (`cacheable(item::DataItem) = item.data isa AbstractDataFrame`, line 279). No other
  Core code.
- `DataBrowserAPI` deps today: `CancellationTokens`, `Dates`, `InteractiveUtils` — no Profiling, no
  DataFrames. `DataBrowserProfiling` is lightweight (JSON, Profile, Statistics).
- `ItemIndex` consumers outside Core: `lib/DataBrowserGUI/src/Browser/{Tags,Operations,Persistence}.jl`,
  `lib/DataBrowserGUI/src/Gui/{TreePanel,Layout}.jl`, `lib/DataBrowserPlots/src/{DataBrowserPlots,PlotPanel}.jl` —
  all via `DataBrowserCore.ItemIndex`.
- `Cache` consumers outside Core: `lib/DataBrowserGUI/src/Browser/Operations.jl`
  (`ProjectCacheSchemaError`), `lib/DataBrowserGUI/src/Gui/Layout.jl` (`ProjectCacheIdentity`),
  `src/DataBrowser.jl:95` (`import DataBrowserCore.Cache as Cache`), tests
  `test_wide_cache.jl`, `test_work_graph.jl`, `test_open_options.jl`, `test_metadata_pipeline.jl`.
- Inside Core, `Workspace.jl` imports ~50 names via `using ..Cache:` (including the unexported
  `ProjectCacheIndex` and `_load_source_item_fingerprints`) and `project_engine.jl:33` imports from
  `.ItemIndex`.

## Decisions already made (do not relitigate)

- **ItemIndex moves into `DataBrowserAPI`** as a submodule (`DataBrowserAPI.ItemIndex`), both files.
  Item records, hierarchy, and source scans are shared vocabulary of the whole family, next to
  `MetadataDict`/`MetadataValue` which already live in API. (User decision, 2026-07-11.)
- **API does not gain a DataFrames dependency.** The one line needing it becomes a dispatch hook:
  API defines `cacheable_data(::Any)::Bool = false` and
  `cacheable(item::DataItem) = cacheable_data(item.data)`; `DataBrowserCache` — which owns the
  columnar storage that makes DataFrames natively cacheable — defines
  `DataBrowserAPI.cacheable_data(::AbstractDataFrame) = true`. Core always loads Cache, so
  behavior is identical in the app. No wrapper helpers beyond this hook; dispatch is the contract.
- **No re-export shims.** Core does not keep a `const ItemIndex = ...` alias or re-export Cache
  names it doesn't use. Every call site updates to the new home (`DataBrowserAPI.ItemIndex`,
  `DataBrowserCache`). GUI adds a direct `DataBrowserCache` dep for the two error/identity types it
  already imports.
- **No migrations and no compatibility.** Nothing persisted changes format; if anything incidental
  does, it just changes.

## Phase 0 — ItemIndex moves to API

1. `git mv lib/DataBrowserCore/src/ItemIndex.jl lib/DataBrowserAPI/src/ItemIndex.jl` and
   `git mv lib/DataBrowserCore/src/ItemIndex lib/DataBrowserAPI/src/ItemIndex` (the `Scanning.jl`
   subdirectory). Include `ItemIndex.jl` from `DataBrowserAPI.jl` **after** the existing includes
   (it uses `Project`, `cacheable`, the metadata types, and the scan-profile hooks).
2. Fix the moved module's imports: `import DataBrowserAPI` / `import DataBrowserAPI: ...` become
   relative (`import ..DataBrowserAPI`, `import ..DataBrowserAPI: ...`) — a package cannot depend
   on itself. `import DataBrowserProfiling as Profiling` stays absolute; add `DataBrowserProfiling`
   to `lib/DataBrowserAPI/Project.toml` `[deps]` (path source not needed — Profiling has no path
   deps; mirror how other lib packages declare it, with a `[sources]` entry like Core's).
3. Replace `cacheable(item::DataItem)::Bool = item.data isa AbstractDataFrame` (ItemIndex.jl:279)
   with the `cacheable_data` hook per the decisions above: the generic
   `cacheable_data(::Any)::Bool = false` lives in `lib/DataBrowserAPI/src/item_contract.jl`
   next to `cacheable`, with a docstring saying storage backends extend it for payload types they
   can store natively. Delete `using DataFrames: AbstractDataFrame` from ItemIndex.jl. Until
   Phase 1 adds the DataFrame method in DataBrowserCache, add it temporarily at the top of
   Core's `Cache.jl` (it moves in Phase 1) so behavior never changes.
4. Update every `ItemIndex` reference:
   - `lib/DataBrowserCore/src/DataBrowserCore.jl:21` — drop `include("ItemIndex.jl")`.
   - `lib/DataBrowserCore/src/project_engine.jl:33`, `lib/DataBrowserCore/src/Workspace.jl:87`,
     `lib/DataBrowserCore/src/Cache.jl:28` — `.ItemIndex`/`..ItemIndex` becomes
     `DataBrowserAPI.ItemIndex`.
   - GUI (`Browser/Tags.jl`, `Browser/Operations.jl`, `Browser/Persistence.jl`, `Gui/TreePanel.jl`,
     `Gui/Layout.jl`) and Plots (`DataBrowserPlots.jl`, `PlotPanel.jl`) —
     `DataBrowserCore.ItemIndex` becomes `DataBrowserAPI.ItemIndex`. Grep the whole repo
     (including `test/` and `bench/`) for `ItemIndex` to catch stragglers.
5. Docs, same commit: CLAUDE.md table row "Scanning, item records, hierarchy" now points to
   `lib/DataBrowserAPI/src/ItemIndex/`; sweep [../data-model.md](../data-model.md) and
   [../api.md](../api.md) for statements that the item index lives in Core.
6. Run the suite.

## Phase 1 — DataBrowserCache becomes a real package

1. Rewrite `lib/DataBrowserCache/src/DataBrowserCache.jl`: delete the `using DataBrowserCore`
   shell. The module body is Core's `Cache.jl` header moved in — the `DataBrowserAPI` imports, the
   Profiling import, `DuckDB`/`DBInterface`/`DataFrames`/`SHA`/`Serialization`/`Dates`, and the
   ItemIndex import (now `import DataBrowserAPI.ItemIndex: ...`) — followed by plain local
   `include("build_metrics.jl")` etc. for the three files, followed by the existing `export` list.
   Add `DataBrowserAPI.cacheable_data(::AbstractDataFrame)::Bool = true` here (with the
   temporary copy in Core's Cache.jl deleted along with that file). Keep the module docstring but
   drop its "re-exports DataBrowserCore.Cache" claim.
2. `lib/DataBrowserCache/Project.toml`: remove `DataBrowserCore`; add `DataBrowserAPI`,
   `DataBrowserProfiling`, `DuckDB`, `DBInterface`, `DataFrames`, `SHA`, `Serialization`, `Dates`
   with `[sources]` path entries for the two DataBrowser packages and compat bounds copied from
   Core's.
3. Delete `lib/DataBrowserCore/src/Cache.jl` and its include in `DataBrowserCore.jl:24`. Core's
   `Project.toml` gains `DataBrowserCache` (dep + `[sources]` path). Then grep Core's `src/` for
   `DuckDB`, `DBInterface`, `SHA`, `Serialization` — drop each from Core's deps if the cache was
   the only user.
4. Update Cache call sites:
   - `lib/DataBrowserCore/src/Workspace.jl` — `using ..Cache:` → `using DataBrowserCache:` and
     `import ..Cache:` → `import DataBrowserCache:` (same name lists; unexported names like
     `ProjectCacheIndex` still resolve through the explicit list). Any remaining `Cache.` qualified
     use gets `import DataBrowserCache as Cache` in that module.
   - `src/DataBrowser.jl:95` — `import DataBrowserCache as Cache`.
   - GUI `Browser/Operations.jl` and `Gui/Layout.jl` — `using DataBrowserCache: ...`; add
     `DataBrowserCache` to `lib/DataBrowserGUI/Project.toml` (dep + source).
   - Tests (`test_wide_cache.jl`, `test_work_graph.jl`, `test_open_options.jl`,
     `test_metadata_pipeline.jl`) — the `DataBrowserCore.Cache` consts become `DataBrowserCache`.
     Grep repo-wide for `DataBrowserCore.Cache` and `..Cache` to catch stragglers.
5. Manifest check, not load assertions: resolve the standalone package
   (`julia --project=lib/DataBrowserCache -e 'using Pkg; Pkg.resolve()'`) and confirm
   `DataBrowserCore` does not appear in `lib/DataBrowserCache/Manifest.toml`. Same check that
   `DataBrowserAPI`'s manifest stays free of DataFrames/DuckDB.
6. Run the suite.

## Phase 2 — Docs and residuals

1. [../cache.md](../cache.md): the cache is owned by `lib/DataBrowserCache/`; Core consumes it
   through `Workspace`. Fix any path references to `lib/DataBrowserCore/src/Cache.jl`.
2. CLAUDE.md table: the DuckDB cache row points at `lib/DataBrowserCache/` only.
3. [roadmap.md](roadmap.md): mark the Cache-edge paragraph executed, pointing here.
4. [gui-extension-architecture.md](gui-extension-architecture.md) target graph already shows
   `DataBrowserCache -> API, Profiling, DataFrames/DuckDB`; add the now-true `Core -> Cache` note if
   the text still hedges.

## Checks

- `lib/DataBrowserCache/Manifest.toml` contains no `DataBrowserCore`; `lib/DataBrowserAPI`'s
  manifest contains no DataFrames, DuckDB, or GUI packages.
- No `joinpath(@__DIR__, "..", "..", ...)` cross-package includes remain anywhere in `lib/`.
- Warm reopen still works: `test_wide_cache.jl` / `test_open_options.jl` pass unchanged in what
  they assert (only their import lines change).
- A DataFrame payload still reports `cacheable(item) == true` once DataBrowserCache is loaded
  (covered by existing cache tests; add one direct `cacheable_data` test only if none fails
  without the hook wired).
