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

Source data is three static CSVs under `bench/templates/` (`kind1.csv`, `kind2.csv`, `kind3.csv`).
`bench/custom_data_source.jl` presents each file as many source items (same bytes, distinct ids).
`scale` / the item-count list only changes how many aliases are advertised.

Every run writes a timestamped directory under `bench/results/` (gitignored). The DuckDB cache lives
in a temp depot and is deleted on exit; only the result files are kept. Compare runs by diffing the
summary artifacts below against `benchmark.log` (realistic) or the printed git context in the
terminal (scaling).

Peak RAM is not sampled inside Julia. Wrap a harness so `ps` records the high-water mark (kilobytes
on macOS and Linux):

```bash
bench/run.sh realistic_browse.jl [scale]
bench/run.sh scaling.jl [n1,n2,...]
```

That writes `peak_rss_kb.txt` in the result directory. If the number jumps by a large factor between
similar runs, RAM usage exploded. Running the `.jl` files directly skips that file.

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

Loading `DataBrowserProfiling` turns on `@timed_dbg` in the engine. A benchmark clears the
accumulator with `DataBrowserProfiling.reset_debug_timings!()` before the measured work, then
calls `DataBrowserProfiling.take_debug_timings!()` to get a `TimerOutputs.TimerOutput` for that interval.
The realistic harness writes that tree as `debug_timings.txt` and `debug_timings.csv`. Use Julia
sampling profiles and pprof for call-path attribution.

## Realistic browsing

```bash
julia --project=bench --threads=auto bench/realistic_browse.jl [scale]
```

The single benchmark harness. It builds a three-kind project against the static templates: compact
alias counts, fatigue-style row volume in `kind3.csv`. It does what a user does — selects items and
**renders real plots** (GLMakie figures built from the cached data) while the cache is still being
built, and again after it settles. The default scale models 636 aliases and about 2,156 items, and
crosses the cache buffer row ceiling; smaller scales are useful for iteration.

It measures, on real functions and real data:

1. **build throughput** — parallel scan, item processing, and collection analysis (items/s each);
2. **interactive plot probe latency** — the full select, materialize, setup, and draw round trip,
   sampled during the build and after it settles, per item kind;
3. **cache writes** — interpreted/processed/stats call counts, mean latency, writer occupancy
   (busy vs queued-wait), and an explicit processed-payload saturation pass;
4. **normalized averages** — rows per alias/item, milliseconds per alias/item, write nanoseconds per
   payload row so scale sweeps are comparable;
5. **warm reopen** — closing and reopening on the same cache, timing the incremental rescan and the
   first plot, with the allocation it costs (the cached-index handling path).

Outputs land under `bench/results/realistic-<timestamp>/`: `benchmark.log` (git branch, commit,
Julia version, threads, env), `scorecard.csv` (the one-line summary),
`responsiveness.csv` (every interactive sample), `saturation.csv`,
`reopen.csv`, and explicit timing summaries in `debug_timings.txt` and `debug_timings.csv`.
`bench/run.sh realistic_browse.jl` also writes `peak_rss_kb.txt`.

**Compare:** `scorecard.csv` is the primary before/after line — build throughput, normalized
ms/file and ms/item, plot latencies, and warm reopen. Use `debug_timings.csv` to compare
instrumented operation totals. `benchmark.log` records the exact git commit and tunables for each run.
