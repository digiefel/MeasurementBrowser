# Profiling & benchmarking

Performance is the reason this app exists, so the engine ships cheap summary metrics, explicit
debug timing summaries, and a headless benchmark harness.

For profiling the whole application on the real RuO2 project (compilation, startup, scan
throughput, overhead, GUI cost), use the `databrowser-profiling` skill
(`.agents/skills/databrowser-profiling/`). The harnesses below run synthetic workloads instead:
use them when you need controlled item counts or metrics the skill does not cover (plot latency,
warm reopen, scaling exponents).

The `bench/` environment is separate from the package (`julia --project=bench`) so dev tools
(BenchmarkTools, CairoMakie, …) do not become runtime dependencies. Instantiate once:

```bash
julia --project=bench -e 'using Pkg; Pkg.instantiate()'
```

Every run writes a timestamped directory under `bench/results/` (gitignored). Synthetic data and
DuckDB caches live in temp dirs and are deleted on exit; only the result files are kept. Compare
runs by diffing the summary artifacts below against `benchmark.log` (realistic) or the printed
git context in the terminal (scaling).

## Scaling sweep

```bash
julia --project=bench bench/scaling.jl [n1,n2,...]
```

Default sizes: `500,1000,2000,4000`. Pass a comma-separated list while iterating — the largest
size builds an O(N²) scan and a wide sweep takes minutes.

Times GUI-hot operations that must **not** grow with item count, plus cold scan-build time normalized
per item. The `scan_build_per_item` exponent is the throughput guard: exponent ~0 means stable
throughput, while exponent ~1 means total scan time is quadratic. The sweep must cross the cache
buffer row ceiling (~1000 items); below it metadata reads look artificially flat.

**Persistent output:** `bench/results/<yyyymmdd-HHMMSS>-scaling/scaling.csv` — one row per operation
with `exponent`, `r2`, and `ms_n<size>` columns for each sweep point.

**Compare:** diff `scaling.csv` between runs. `status_refresh` and `items_panel` should stay near
exponent 0; `metadata_publish` per-call cost grows once past the buffer ceiling (known O(N) per
publish × N publishes during scan).

## Debug timings

Benchmark code can wrap a selected phase in `with_debug_timings(timings) do ... end` and save the
result with `write_debug_timings(outdir, timings)`. The output is a readable timing table and CSV;
it aggregates task-owned timers after the benchmark reaches its intended idle point. Use Julia
sampling profiles and pprof for call-path attribution.

## Realistic browsing

```bash
julia --project=bench --threads=auto bench/realistic_browse.jl [scale]
```

The single benchmark harness. It builds a deterministic three-kind project with compact file count
but beyond-RuO2 fatigue pressure: each synthetic big file expands into more cycles and rows per cycle
than the typical real fatigue file. It does what a user does — selects items and **renders real
plots** (GLMakie figures built from the cached data) while the cache is still being built, and again
after it settles. The default scale models 636 files and about 2,156 items, but still crosses the
cache buffer row ceiling; smaller scales are useful for iteration. The `bench/` environment is
separate from the package so its dev tools (BenchmarkTools, CairoMakie, …) don't become runtime
dependencies.

It measures, on real functions and real data:

1. **build throughput** — parallel scan, item processing, and collection analysis (items/s each);
2. **interactive plot probe latency** — the full select, materialize, setup, and draw round trip,
   sampled during the build and after it settles, per item kind;
3. **cache writes** — interpreted/processed/stats call counts, mean latency, writer occupancy
   (busy vs queued-wait), and an explicit processed-payload saturation pass;
4. **normalized averages** — rows per file/item, milliseconds per file/item, write nanoseconds per
   payload row, and memory per item so scale sweeps are comparable;
5. **process memory** — RSS, GC-live, and index/queue counts sampled across the build;
6. **warm reopen** — closing and reopening on the same cache, timing the incremental rescan and the
   first plot, with the allocation it costs (the cached-index handling path);
7. **database aggregation query latency** — one `sum/avg` query per processed payload schema.

Outputs land under `bench/results/realistic-<timestamp>/`: `benchmark.log` (git branch, commit,
Julia version, threads, env), `scorecard.csv` (the one-line summary),
`responsiveness.csv` (every interactive sample), `memory_samples.csv`, `saturation.csv`,
`database_aggregation_queries.csv`, `reopen.csv`, `benchmark.log`, and explicit timing summaries
in `debug_timings.txt` and `debug_timings.csv`.

**Compare:** `scorecard.csv` is the primary before/after line — build throughput, normalized
ms/file and ms/item, plot latencies, memory, and warm reopen. Use `debug_timings.csv` to compare
instrumented operation totals. `benchmark.log` records the exact git commit and tunables for each run.
