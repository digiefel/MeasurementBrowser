---
name: databrowser-profiling
description: >-
  Profile DataBrowser startup, GUI, and RuO2 fresh or resumed-cache performance in a
  persistent Revise session, with numeric timeline output, ProfileView, and
  pprof text/web navigation. Use for scan stalls, throughput, GUI render-loop
  noise, CPU or wall-time hotspots, profile filtering, and iterative performance
  work in this repo.
---

# DataBrowser profiling

## Goal

Keep one DataBrowser session open while using the app, profiling it, editing source with Revise,
and checking whether the app actually becomes faster. Save each run as files the user can inspect
directly.

## Measurement contract

User code is the code registered with the project API: file detection, reading, splitting a source
into entries, processing, analysis, and collection callbacks. Everything DataBrowser does around
those callbacks is application overhead, including repeated reads, scheduling, cache work,
coordination, publishing results, and GUI updates.

For a parallel measurement:

```text
useful work = elapsed time summed over all registered callback calls
total capacity = measurement duration × number of pipeline workers
overhead fraction = 1 - useful work / total capacity
```

Only report this percentage when a script has measured every registered callback kind and uses the
configured pipeline-worker count. The current callback summary omits collection callbacks, so it
does not yet provide a complete overhead percentage. Present the measured callback times and leave
the percentage unavailable; do not reconstruct it manually from an incomplete table.

The scripts currently measure these values directly:

- `Pkg.precompile()` wall time in a fresh process;
- `using DataBrowser` wall time;
- loading the copied project script, which includes its helper files and `define_project`/
  `register_*` calls;
- total `open_workspace` wall time;
- marked DataBrowser calls such as cache opening, source discovery, source publication, item
  publication, and cache writes in `debug_timings.txt` and `debug_timings.csv`;
- time from `open_browser` until the first rendered frame;
- progress and throughput by sources, items, and collections during a fixed measurement period;
- runtime compilation during that period;
- elapsed time inside the registered detect/read/entries/process/analyze callbacks;
- sampled call stacks for finding which functions ran during application overhead.

The project-script time is not a pure measurement of the `register_*` calls because the script also
loads and compiles its helper code. Say that plainly when reporting it.

`profile_scan!` owns one explicit `DebugTimings` object around workspace opening and the measured
scan. It closes the workspace before writing the files so calls still finishing at the boundary are
included safely. Read `debug_timings.txt` for human-readable durations and `debug_timings.csv` for
programmatic comparison. The rows are durations of marked function calls, not sampled estimates.
Parent rows include marked child calls, and calls on different tasks can overlap, so row totals may
exceed the measurement wall time. If a requested operation has no marked row, use the sampling
profile to locate it or state that its duration was not recorded.

## Workflow

1. Check that the copied RuO2 project still matches the real project.
2. Start one Revise-backed Julia session and open the browser once.
3. Use the app normally and profile the period that contains the behavior under investigation.
4. Give the user direct links to the summary, timeline, and profile files.
5. Use pprof to find the first meaningful DataBrowser call path and explain it in plain language.
6. Edit one relevant piece of code, profile the same behavior again, and compare the files.
7. Rebuild the cache from the open app when a fresh build is needed. Close the session when done.

Run only one profiling workspace at a time. The scripts use
`$TMPDIR/databrowser-profiling.pid` to prevent overlapping runs.

## Check the workload copy

`project/definitions.jl` and its included files copy the user's RuO2 v2 project:

```bash
USER="/Users/davide/Documents/OneDrive/OneDrive - Lund University/projects/Borg/202501_RuO2test/analysis/v2"
COPY=".agents/skills/databrowser-profiling/project"
DATA="${RUO2_DATA_ROOT:-/Users/davide/Library/CloudStorage/OneDrive-LundUniversity/projects/Borg/202501_RuO2test/electricaldata}"
test -d "$DATA"
diff -qr "$COPY/data" "$USER/data"
diff -qr "$COPY/analysis" "$USER/analysis"
diff -qr "$COPY/plots" "$USER/plots"
```

Also compare the registered item and plot kinds in `COPY/definitions.jl` with `USER/browser.jl`.
If a difference changes the workload, ask whether to update the copy before profiling.

## Start the persistent session

Revise struct revision is enabled in `bench/LocalPreferences.toml`, so method and type edits
normally do not require restarting Julia.

Start the bench session with the thread count used for the comparison:

```bash
tmux -L jmux set-environment -g JULIA_NUM_THREADS auto

jmux --restart --project /Users/davide/code/Julia/DataBrowser/bench \
  'using Revise; @includet "/Users/davide/code/Julia/DataBrowser/.agents/skills/databrowser-profiling/scripts/interactive_profile.jl"'
```

Use the same Julia thread count for every run being compared.

## Cache state

The copied project is named `DataBrowserProfilingRuO2`, so it uses the normal isolated cache:

```text
DEPOT_PATH[1]/databrowser/DataBrowserProfilingRuO2/cache.duckdb
```

- `cache_mode=:fresh` opens with `rebuild=true` and starts this cache again from zero.
- `cache_mode=:resume` opens the cache already saved by an earlier run.

The public `workspace_status` API supplies the cached-source and interpreted-item counts written to
the summaries and timelines. For a partial-cache check, compare the final counts from the fresh run
with the starting counts from the resume run. A resume must not go backwards. If the values differ,
report the values instead of assuming that progress was preserved; work may have completed between
the last sample and workspace shutdown.

For a controlled fresh/resume pair:

```bash
jmux --project /Users/davide/code/Julia/DataBrowser/bench \
  'DataBrowserInteractiveProfile.profile_scan!(mode=:headless, budget=10, profiler=:none, cache_mode=:fresh)'

jmux --project /Users/davide/code/Julia/DataBrowser/bench \
  'DataBrowserInteractiveProfile.profile_scan!(mode=:headless, budget=60, profiler=:cpu, cache_mode=:resume)'
```

The harness uses DataBrowser's real default `background_processing=false`.

## Live profiling loop

Open the real GUI workspace and browser once:

```bash
jmux --project /Users/davide/code/Julia/DataBrowser/bench \
  'DataBrowserInteractiveProfile.start_live!(mode=:gui, cache_mode=:resume)'
```

Use the app, then profile the behavior of interest. Measure at least 60 seconds after the first
frame. The command below records a fixed 60 seconds even if the app is idle for part of that time:

```bash
jmux --project /Users/davide/code/Julia/DataBrowser/bench \
  'DataBrowserInteractiveProfile.profile_live!(budget=60, profiler=:cpu)'
```

Use `profiler=:wall` when the question is blocking or waiting rather than CPU work. Use
`profiler=:none` when comparing throughput without profiler overhead.

`profile_live!` produces the fixed-window progress and sampling files. Use `profile_scan!` when the
question requires the marked cache-open, scan, or publication durations in `debug_timings.txt`.

The user can refresh ProfileView or pprof from the same captured samples:

```bash
jmux --project /Users/davide/code/Julia/DataBrowser/bench \
  'DataBrowserInteractiveProfile.refresh_profileview!()'

jmux --project /Users/davide/code/Julia/DataBrowser/bench \
  'DataBrowserInteractiveProfile.refresh_pprof!()'
```

After a source edit, the next `profile_live!` applies pending Revise changes before recording. To
start the cache again while keeping the same Julia process and browser, click **Rebuild Cache** in
the app or call:

```bash
jmux --project /Users/davide/code/Julia/DataBrowser/bench \
  'DataBrowserInteractiveProfile.rebuild_live_cache!()'
```

Close everything owned by the profiling session when finished:

```bash
jmux --project /Users/davide/code/Julia/DataBrowser/bench \
  'DataBrowserInteractiveProfile.stop_live!()'
```

## Files from each run

Each `profile_live!` results directory contains:

- `summary.txt`: setup time, the profiled duration, progress counters, compilation, and project
  callback time;
- `timeline.csv` and `timeline.png`: cached sources and completed items over time;
- `*_profile.jls`: the most complete saved Julia sample and frame data;
- `*_profile.jlprof`: the saved profile opened by ProfileView;
- `*_profile.pb.gz`: the same profile converted for pprof.

`profile_scan!` additionally writes:

- `debug_timings.txt`: a readable table of marked DataBrowser calls;
- `debug_timings.csv`: the same call counts and durations for comparisons.

The summary uses these terms:

- **measurement period**: the fixed wall-clock interval after the first frame in GUI runs;
- **runtime compilation**: compilation performed during the measurement, summed across Julia
  threads; it can exceed wall time when several threads compile concurrently;
- **project callback time**: elapsed time summed across callback calls. `detect` chooses a recipe,
  `read` loads a source, `entries` splits it into items, `process` transforms an item, and `analyze`
  computes item results. Calls may overlap on different workers, so their sum can exceed wall time.

Run the pprof commands in [references/profile-tools.md](references/profile-tools.md) and save their
text beside these files. When replying, provide direct links to the summary, timeline image, pprof
text, ProfileView file, raw profile, and any debug-timing files. Present the files instead of
retyping their tables into chat. Add a short plain-language explanation grounded in the DataBrowser
operation the user recognizes—scanning files, reading measurements, building the item tree, drawing
the GUI—not just Julia runtime function names.

Do not force every run into a fixed verbal template. Report what the files show clearly, state what
the app was doing during the profile, and say when the samples did not capture the intended work.

Profiled and unprofiled throughput are not directly comparable. Compare runs with matching thread
count, cache state, mode, profiler setting, and user action.

## One-shot and fresh-process checks

`profile_scan!` is the secondary one-shot check. It includes `open_workspace`, optionally opens the
GUI, profiles the scan, writes the sampling, timeline, and debug-timing files, and closes the
workspace.

Use `scripts/timed_profile.jl` only when startup in a fresh Julia process is itself the question:

```bash
julia --project=bench --threads=auto \
  .agents/skills/databrowser-profiling/scripts/timed_profile.jl gui \
  --budget=60 --profile=cpu --cache=fresh
```

Use `--cache=resume` for a fresh-process reopen of the dedicated cache. Record the commit and Julia
thread count with final comparisons.

## Scaling

Use `bench/scaling.jl` for the existing controlled N sweep. It uses one source per item, so S and N
change together; it measures total scan/build time per item and the scaling of status refresh, item
panel gathering, and metadata publication. It does not independently separate source-count scaling
from item-count scaling.

Use `bench/realistic_browse.jl` for a synthetic workload with different numbers of files and items
per file. Its `scorecard.csv` reports source and item throughput, and its `debug_timings.txt` and
`debug_timings.csv` report marked cache, scan, publication, and plotting calls. A proper independent
S/N scaling study requires a set of runs that varies file count while holding items per file fixed,
then varies items per file while holding file count fixed. Do not describe a single scale-factor run
as that study.

## Validity

Do not claim an improvement by doing less work, hiding samples, or comparing different cache states.
Put the progress counts and timeline beside elapsed time. If the profile mostly contains waiting or
Julia runtime frames, say that the chosen interval did not isolate useful application work and use
the navigation steps in `references/profile-tools.md` before choosing code to change.
