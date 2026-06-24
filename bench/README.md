# Profiling & benchmarking

Performance is the reason this app exists, so the engine ships permanent, zero-cost-when-off
instrumentation plus a headless benchmark harness.

## Engine timers (`MeasurementBrowser.Profiling`)

Hot **sequential** regions are wrapped with `@timeit_debug TIMER "label" …`. They compile to nothing
until you enable them, so normal runs pay nothing.

```julia
using MeasurementBrowser
MeasurementBrowser.Profiling.enable!()   # turn on (do this at top level — see note)
# … exercise a workload …
MeasurementBrowser.Profiling.report()    # print per-region time / calls / allocations
MeasurementBrowser.Profiling.disable!()
```

Note: `enable!()` redefines methods, so calls made *in the same function* that called `enable!()`
won't see it yet (Julia world age). Enable at top level (or in a separate REPL step) before running
the workload. The parallel scan is **not** instrumented this way (a shared timer isn't thread-safe);
its per-file parse cost is in `MeasurementBrowser.scan_profile_summary(project)` instead.

## Benchmark harness

```bash
julia --project=bench --threads=auto bench/cache_pipeline.jl [n_files] [rows_per_file]
```

Creates a deterministic synthetic CSV dataset once under `bench/data/`, reuses it on later runs,
and saves each report under `bench/results/`. The harness reports real functions on real data:

1. **per-item micro-benchmarks** — CSV parse vs serialize vs deserialize vs DuckDB blob round-trip
   vs DuckDB columnar round-trip (this settles the fastest way to cache item data);
2. **whole-dataset macro timings** — the real parallel source scan;
3. **source-item expansion** — copied child tables versus views;
4. **workspace cold builds** — full scan, item stats, native DuckDB writes, and cache-only reads;
5. **incremental reopen** — no-change and one-changed-source read counts and elapsed time;
6. **item-analysis throughput** — one source item expanded into many independently processed items;

The `bench/` environment is separate from the package so its dev tools (BenchmarkTools, …) don't
become runtime dependencies.

## Realistic browsing

```bash
julia --project=bench --threads=auto bench/realistic_browse.jl [scale]
```

Builds a deterministic three-kind project with thousands of files and measures real item
materialization while the cache is being built and after it settles. The default scale models 4,500
files and about 10,400 items; smaller scales are useful for iteration. `scorecard.csv` records scan,
processing, and end-to-end throughput; mean write latency and writer occupancy; and aggregate
during-build and after-build median, p90, p99, and maximum read latency.
