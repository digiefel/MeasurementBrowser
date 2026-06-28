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

## Realistic browsing

```bash
julia --project=bench --threads=auto bench/realistic_browse.jl [scale]
```

The single end-to-end harness. It builds a deterministic three-kind project with thousands of files
and does what a user does — selects items and **renders real plots** (GLMakie figures built from the
cached data) while the cache is still being built, and again after it settles. The default scale
models 4,500 files and about 10,400 items; smaller scales are useful for iteration. The `bench/`
environment is separate from the package so its dev tools (BenchmarkTools, CairoMakie, …) don't
become runtime dependencies.

It measures, on real functions and real data:

1. **build throughput** — parallel scan, item processing, and collection analysis (items/s each);
2. **interactive plot latency** — the full select → materialize → setup → draw round trip, split into
   data-load vs render time, sampled during the build and after it settles, per item kind;
3. **cache writes** — interpreted/processed/stats call counts, mean latency, and writer occupancy
   (busy vs queued-wait);
4. **process memory** — RSS, GC-live, and index/queue counts sampled across the build;
5. **warm reopen** — closing and reopening on the same cache, timing the incremental rescan and the
   first plot, with the allocation it costs (the cached-index handling path);
6. **whole-table aggregates** — one `sum/avg` per processed payload schema, the analysis-query shape.

Outputs land under `bench/results/<timestamp>/`: `scorecard.csv` (the one-line summary),
`responsiveness.csv` (every interactive sample), `memory_samples.csv`, `global_queries.csv`,
`reopen.csv`, and — when CairoMakie is present — `responsiveness.png`. Diagnostics that add locks or
flushes are disabled in the normal run.
