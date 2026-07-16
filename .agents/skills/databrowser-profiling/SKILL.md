---
name: databrowser-profiling
description: Profile and benchmark the DataBrowser application on the real RuO2 v2 project, headless or with the native GUI. Use for DataBrowser startup, compilation, throughput, overhead, scan-timeline, or GUI performance investigations in this repository.
---

# DataBrowser application profiling

Use the script `scripts/timed_profile.jl` to measure DataBrowser's performance. It runs the real
RuO2 project on the real measurement data, times every phase from precompilation to a live scan,
and prints one summary at the end. A full run should take about two minutes. To see the effect of a
code change, run it once before the change and once after, and compare the two summaries.

Run it from the repo root:

```bash
julia --project=bench --threads=auto .agents/skills/databrowser-profiling/scripts/timed_profile.jl headless
julia --project=bench --threads=auto .agents/skills/databrowser-profiling/scripts/timed_profile.jl gui
```

`headless` runs the engine without a window and works in any environment. `gui` additionally
opens the real DataBrowser window, which only works in an unsandboxed process. To find out how
much the GUI slows down the engine, run both modes with the same settings and compare their
throughput numbers.

There are four options:

- `--budget=SECONDS` sets how long the measured scan runs. Default 30.
- `--warmup=SECONDS` sets how long the throwaway warmup scan runs before the measured one, so
  that compilation happens during warmup instead of polluting the measurement. Default 20. Set
  it to 0 if you specifically want to measure cold-start behavior.
- `--no-fresh` skips deleting the DataBrowser compiled caches. The run starts faster, but the
  precompile number then measures an empty rebuild and means nothing.
- `--profile=cpu` or `--profile=allocs` additionally runs Julia's sampling profiler during the
  measured scan and writes a flat report (`cpu_profile.txt` or `allocation_profile.txt`) next to
  `timeline.csv`. Use this to find which functions the scan spends its time or allocations in.
  Profile the two separately: allocation sampling changes runtime cost.

## How to read the summary

The summary has five sections. This is what each line tells you:

**Startup.** How long each phase took, in order: `precompile` rebuilds the nine DataBrowser
packages from scratch (the per-package times are printed above the summary by Pkg), `using
DataBrowser` loads the package, `project include` runs the project file with all its `register!`
calls, and `open_workspace` creates the workspace. In gui mode there is also `open_browser`,
which times how long starting the render loop takes. It does not time how long until the window
actually responds; the app has no hook for that yet.

**Scan window.** The measured scan runs on an empty temporary cache, so it always does the full
work. `runtime JIT in window` is compilation that happened during the measurement; with the
default warmup it should be around a second, and if it is much larger the warmup was too short.
The four milestone lines tell you when discovery, interpretation, processing, and analysis
finished, or `not reached` if the budget ran out first. With a 30-second budget on the RuO2 data
you will normally only see interpretation running; processing is lower priority by design and
starts later.

**Throughput.** Items per second for interpretation and processing, given separately for the
first and second half of the window. If the two halves differ a lot, look at `timeline.csv` to
see whether the rate actually degrades or whether the scan just hit a stretch of heavier files.

**Pipeline callbacks vs capacity.** This answers how much time goes into the actual work versus
the machinery around it. The engine records the time spent inside the project's own callbacks
(detect, read, entries, process, analyze). The summary compares that total against the window
length times the number of worker threads. In a saturated scan, the difference is the overhead
of everything the app does around the callbacks: scheduling, caching, publishing.

**Timeline events.** Wall-clock timestamps of every phase of the run, so you can see where the
two minutes went.

The per-sample counts (every 0.25 s: sources found, interpreted, processed, analyzed, and so on)
are written to `timeline.csv` in `bench/results/<timestamp>-ruo2-v2-timed-<mode>/`.

## What the script uses

The project definitions are a copy of the RuO2 v2 project, checked in at
`project/definitions.jl` together with its `data/`, `analysis/`, and `plots/` files. If the real
project changes, copy the files over again. The measurement data itself is not copied: the
script reads the real `electricaldata` directory (6321 files, about 65k items). Set
`RUO2_DATA_ROOT` to use a different data root.

The script never touches your real caches. Workspace caches go to a temporary directory that is
deleted afterwards, and a fresh run deletes only the `DataBrowser*` entries from the compiled
cache, which Julia rebuilds in the precompile step.

## When this script is not the right tool

- To check that GUI-hot operations stay flat as item counts grow, run `bench/scaling.jl`. It
  sweeps synthetic workspaces of increasing size and fits scaling exponents.
- To benchmark one isolated function properly, use the `julia-bench` skill.
