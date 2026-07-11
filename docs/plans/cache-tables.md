# Tables-Native Cache — Execution Plan

Make the cache speak the Tables.jl interface instead of requiring DataFrames: any tabular payload
goes in, the *same container type* comes back out, and the user never notices the database in
between. This also drops DataFrames from `DataBrowserCache`'s dependencies — DuckDB.jl itself
depends only on `Tables`, so nothing else in the cache's environment pulls DataFrames in.

Written to be executed phase by phase; **the full test suite must pass after every phase.**
Run tests with `julia --project --threads=4 -e 'using Pkg; Pkg.test()'`.

## Background (verified 2026-07-11)

Tables.jl defines no container. It is an interface — functions like `Tables.columns`,
`Tables.getcolumn`, `Tables.schema`, `Tables.rowcount` — that many containers implement:
`DataFrame`, `NamedTuple` of vectors, DuckDB query results, TypedTables, StructArrays. There is no
`Table` supertype and no `DataFrame <: Table`. "Moving to Tables" therefore means two separate
things, both in scope here:

1. **Internally**, the cache stops holding columns in `DataFrame`s and uses Base's NamedTuple of
   vectors (`Tables.columntable`) where it needs a container of its own.
2. **At the payload boundary**, the cache accepts any object implementing the interface and
   reconstructs the caller's original container on read via `Tables.materializer(T)` — the one
   function Tables.jl provides for exactly this.

DataFrames use in the cache today, all replaceable:

- `DataFrame(DBInterface.execute(...))` (cache_buffer.jl:1028) → `Tables.columntable(...)`.
- `names`/`nrow`/`eltype` (`_dataframe_shape` at cache_buffer.jl:90, rows accounting at :577, :891;
  project_cache_domain.jl:574) → `Tables.schema`, `Tables.rowcount`.
- `body[!, name]` column access in the Appender write loop (cache_buffer.jl:1300) →
  `Tables.getcolumn`.
- `groupby(rows, MB_SEQ_COLUMN)` + `DataFrame(group[:, data_columns])` (cache_buffer.jl:1041-1043),
  splitting one bulk payload read into per-item tables → a hand-rolled split over the seq column
  (~20 lines; the query already orders by seq).
- The payload gate `payload isa AbstractDataFrame` (project_cache_domain.jl:726).

## Decisions

- **`cacheable_data` default becomes `Tables.istable`.** The API fallback
  `cacheable_data(::Any) = false` changes to `cacheable_data(data) = Tables.istable(data)`;
  `DataBrowserAPI` gains the (interface-only, tiny) `Tables` dep. Core's
  `cacheable_data(::AbstractDataFrame)` method is subsumed and deleted. Tables are first-class:
  anything tabular is cacheable by default, and a type can still opt out by dispatch.
- **Column eltypes must map to DuckDB SQL types** (`_duckdb_sql_type`). The gate checks this;
  a tabular payload with an unmappable column eltype is treated as non-cacheable (memory-only),
  exactly like non-tabular data today. No serialized-blob fallback in this pass.
- **Round-trip returns the original container type.** The stored schema row gains the payload's
  container type, serialized like fingerprints (hex text via `Serialization`); read rebuilds with
  `Tables.materializer(T)(columntable)`. If the type cannot be deserialized (project code changed),
  the entry is a cache miss and the work reruns — the cache may always be rebuilt.
- **Schema version bumps once** for the type-tag column and any renamed tables; old caches rebuild.
  No migration. Internal names saying "dataframe" (`_dataframe_shape`, `_dataframe_table_name`,
  `dataframe_schemas`, `DataFrameShape`) are renamed to say what they now are (tabular payloads).

## Phase 0 — Internal container swap

1. In `cache_buffer.jl`: replace `_dataframe_shape` with a `Tables.schema`-based helper; loosen
   `TabularBodyBatch.body` and the read-path types from `AbstractDataFrame` to `Any` (the contract
   is the interface, checked at the gate); rewrite the bulk-read split (columntable + seq-column
   ranges) and the Appender loop (`Tables.getcolumn`). Delete `using DataFrames` from
   `DataBrowserCache.jl`; add `Tables`.
2. In `project_cache_domain.jl`: rows accounting via `Tables.rowcount`; rename the "dataframe"
   identifiers; bump `PROJECT_CACHE_SCHEMA_VERSION`.
3. `lib/DataBrowserCache/Project.toml`: drop `DataFrames`, add `Tables`. Reads still return
   NamedTuple columntables internally; `Workspace.materialize_items` behavior is unchanged for
   DataFrame payloads only because Phase 1 restores the container — so Phases 0 and 1 land
   together if the suite catches the intermediate container change, separately if not.

## Phase 1 — Payload boundary

1. Gate: `cacheable(item) && Tables.istable(payload) && eltypes map to SQL types`.
2. Store the container type with the schema row; rebuild via `Tables.materializer` on read;
   unresolvable type = miss.
3. API: `cacheable_data` default delegates to `Tables.istable`; Core's DataFrame method deleted;
   API `Project.toml` gains `Tables`.

## Phase 2 — Tests and docs

1. Round-trip tests: DataFrame in → DataFrame out; NamedTuple of vectors in → NamedTuple out
   (including through warm reopen); non-tabular and unmappable-eltype payloads stay memory-only.
2. `docs/cache.md` ownership/vocabulary sections describe the Tables contract; AGENTS.md row
   unchanged.

## Checks

- `lib/DataBrowserCache/Manifest.toml` standalone-resolves without `DataFrames`.
- A payload cached before the change is not silently misread (schema bump forces rebuild).
- `read_item_data` of a NamedTuple payload returns a `NamedTuple`, not a DataFrame.
