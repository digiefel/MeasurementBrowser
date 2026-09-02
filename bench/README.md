# Profiling & benchmarking

Performance is the reason this app exists. The numbers that should stay stable or improve live in
`bench/status.txt` (committed). If a run fails, that file is left empty so a bad baseline is not
committed.

Unit tests and benchmarks share the same `bench/` environment, i.e. this directory.
Note that `Pkg.test()` is not the test command, it would compile a second environment.

For a benchmark-only run:
```bash
bench/run.sh
```

To run tests and benchmark:
```
julia --project=bench --threads=auto test/runtests.jl
```

`test/runtests.jl` runs the unit tests, then `run_performance()` which writes `status.txt`.
`bench/run.sh` is the same command plus peak RSS via `ps`, appended as `peak_rss_kb`.

We try to keep the performance testing under two minutes to allow fast iteration.


Every line in `status.txt` has a `#` comment that defines the field. Those comments are written by
`bench/run.jl`. Sampled plot times are median then variance; a header row sits above them. Compare
with `git diff bench/status.txt`.

Default workload: realistic scale `0.05` (32 aliases, ~127 items) and scaling sizes `250,500,1000`.
Override with `MB_BENCH_SCALE` and `MB_BENCH_SCALING_SIZES`. Scale `1.0` triggers an optional
heavier run and is the only scale that guarantees processed-writer saturation past the
cache row ceiling.

Source data is three static CSVs under `bench/templates/`. `bench/custom_data_source.jl` presents
each file as many source items. The DuckDB cache lives in a temp depot and is deleted on exit.

Standalone scripts (print to stdout, do not write `status.txt`):

```bash
julia --project=bench --threads=auto bench/realistic_browse.jl [scale]
julia --project=bench bench/scaling.jl [n1,n2,...]
```
