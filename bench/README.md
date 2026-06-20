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

Generates a synthetic CSV dataset and reports, from real functions on real data:

1. **per-item micro-benchmarks** — CSV parse vs serialize vs deserialize vs DuckDB blob round-trip
   vs DuckDB columnar round-trip (this settles the fastest way to cache item data);
2. **whole-dataset macro timings** — the real parallel scan and a single-threaded analysis loop;
3. **engine per-region timings** — the `Profiling` report over a representative cache exercise.

The `bench/` environment is separate from the package so its dev tools (BenchmarkTools, …) don't
become runtime dependencies.
