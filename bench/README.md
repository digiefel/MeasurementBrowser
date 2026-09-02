# Profiling & benchmarking

Performance is the reason this app exists. The numbers that should stay stable or improve live in
`bench/status.txt` (committed). If a run fails, that file is left empty so a bad baseline is not
committed.

Correctness tests stay `Pkg.test()` on the package. This directory is the performance run:

```bash
julia --project=bench -e 'using Pkg; Pkg.instantiate()'
bench/run.sh
```

`bench/run.sh` starts one Julia process (`bench/run.jl`), records peak RSS with `ps` (kilobytes on
macOS and Linux), and appends `peak_rss_kb` after a successful write. Wall target is under two
minutes on a warm cache.

Every line in `status.txt` has a `#` comment that defines the field. Those comments are written by
`bench/run.jl` (that file is the catalog). Sampled plot times are median then variance; a header row
sits above them. Compare with `git diff bench/status.txt`.

Default workload: realistic scale `0.05` (32 aliases, ~127 items) and scaling sizes `250,500,1000`.
Override with `MB_BENCH_SCALE` and `MB_BENCH_SCALING_SIZES`. Scale `1.0` is a heavier optional run
and is the only scale that requires processed-writer saturation past the cache row ceiling.

Source data is three static CSVs under `bench/templates/`. `bench/custom_data_source.jl` presents
each file as many source items. The DuckDB cache lives in a temp depot and is deleted on exit.

Standalone scripts (print to stdout, do not write `status.txt`):

```bash
julia --project=bench --threads=auto bench/realistic_browse.jl [scale]
julia --project=bench bench/scaling.jl [n1,n2,...]
```

For profiling the real RuO2 project, use the `databrowser-profiling` skill
(`.agents/skills/databrowser-profiling/`).
