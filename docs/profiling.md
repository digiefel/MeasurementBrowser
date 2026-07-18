# Profiling

Normal workspace diagnostics keep lightweight build metrics, callback timing rows, workspace health,
and plot-phase samples. They support the browser and do not retain engine event traces.

## Explicit debug timings

Benchmarks and diagnostics can time selected engine calls with one explicit object:

```julia
timings = DebugTimings(start_ns=time_ns())

with_debug_timings(timings) do
    workspace = open_workspace(project, root)
    wait_workspace_idle!(workspace)
end

write_debug_timings(outdir, timings)
```

`@time_dbg` is debug-only. A call such as `Profiling.@time_dbg open_cache_db(identity)` uses
`open_cache_db` as its label, executes the original call directly when debug timing is inactive, and stores
results only in the active `DebugTimings` object. `write_debug_timings` writes `debug_timings.txt`
and `debug_timings.csv` with call counts, inclusive total and average time, and process allocation
changes.

Timers are task-owned because engine work runs concurrently and tasks can migrate across Julia
threads. Enter `with_debug_timings` before spawning the work to measure and collect after the
workspace has reached its intended idle point. Row totals can exceed wall time when calls on several
tasks overlap. Allocation changes come from Julia's process allocation counter, so concurrent rows can
overlap; use a controlled single-worker pass or Julia allocation profiling for attributable allocations.

## Sampling profiles

Julia's `Profile` tools remain the attribution path for a hot aggregate timing label. The plot
diagnostic keeps its bounded sampling helper, and benchmark/skill scripts continue to save
ProfileView and pprof artifacts for CPU or wall-time investigations.
