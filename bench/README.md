# Profiling & benchmarking

Performance is the reason this app exists, so the engine ships cheap aggregate metrics, opt-in
structured profiling, and a headless benchmark harness.

## Structured profiling

The realistic benchmark can record the same workspace-owned trace as the internal GUI:

```bash
MB_PROFILE_INTERNAL=1 julia --project=bench --threads=auto bench/realistic_browse.jl
MB_PROFILE_INTERNAL=1 MB_PROFILE_CPU=1 julia --project=bench --threads=auto bench/realistic_browse.jl 0.1
MB_PROFILE_INTERNAL=1 MB_BENCH_START_PROFILE=0 julia --project=bench --threads=auto bench/realistic_browse.jl
```

The first command captures the full default workload without CPU sampling. The second runs the 10%
CPU-sampled workload. The third measures enabled-but-idle overhead without recording events.
`profile.json` is Perfetto-compatible, `profile_summary.csv` contains grouped latencies, and
`scorecard.csv` records event and drop counts. See [profiling.md](../docs/profiling.md) for runtime
flags and data policy.

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
during-build and after-build median, p90, p99, and maximum read latency. Diagnostics that add locks or
flushes are disabled in the normal run.
