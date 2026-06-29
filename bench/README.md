# Profiling & benchmarking

Performance is the reason this app exists, so the engine ships cheap summary metrics, opt-in
structured profiling, and a headless benchmark harness.

## Structured profiling

The realistic benchmark records the same workspace-owned trace as the internal GUI by default:

```bash
MB_PROFILE_INTERNAL=1 MB_PROFILE_CPU=1 julia --project=bench --threads=auto bench/realistic_browse.jl 0.1
MB_PROFILE_INTERNAL=1 MB_BENCH_START_PROFILE=0 julia --project=bench --threads=auto bench/realistic_browse.jl
MB_PROFILE_INTERNAL=0 julia --project=bench --threads=auto bench/realistic_browse.jl
```

The first command adds CPU samples to a 10% workload. The second measures enabled-but-idle overhead
without recording events. The third disables structured profiling and therefore skips event
histograms. `profile.json` is Perfetto-compatible, `profile_summary.csv` contains grouped latencies,
and `scorecard.csv` records event and drop counts. See [profiling.md](../docs/profiling.md) for
runtime flags and data policy.

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

Outputs land under `bench/results/<timestamp>/`: `scorecard.csv` (the one-line summary),
`responsiveness.csv` (every interactive sample), `memory_samples.csv`, `saturation.csv`,
`database_aggregation_queries.csv`, `pipeline_event_times.csv`, `pipeline_event_summary.csv`,
`pipeline_timeseries.csv`, `reopen.csv`, `benchmark.log`, and — when CairoMakie is present —
`pipeline.png`. Diagnostics that add locks or flushes are disabled in the normal run.
