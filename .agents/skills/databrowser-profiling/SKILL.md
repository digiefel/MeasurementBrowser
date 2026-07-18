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

Use one Revise-backed Julia session to move quickly from a measured stall to a named hot call path.
Give the user the numbers on every pass. Use the same saved samples for the agent's pprof tables and
the user's ProfileView or pprof GUI check.

## Loop

```
- [ ] Diff the copied RuO2 project against the user's current project
- [ ] Start one Revise-backed bench session and call `start_live!` once
- [ ] Let the user use the real GUI; capture repeated fixed windows with `profile_live!`
- [ ] Human: refresh ProfileView or pprof web; agent: query the same `.pb.gz` with pprof CLI
- [ ] Report the required numbers and artifact paths to the user
- [ ] One hypothesis, one edit, one matching window; rebuild the cache only when needed
- [ ] Call `stop_live!` when the iteration is finished
```

Run only one workspace profile at a time. The scripts use a PID-aware
`$TMPDIR/databrowser-profiling.pid` lock to reject overlapping live runs and recover after a killed
Julia process.

In a skill-evaluation loop, the runner is the only agent that executes this skill. Return its
stdout and artifact paths unchanged for review. The reviewer must judge those outputs without
reading this skill or its code.

## Check the workload copy

`project/definitions.jl` and its included files copy the user's RuO2 v2 project. Check them before
measuring so the harness does not optimize a stale workload.

```bash
USER="/Users/davide/Documents/OneDrive/OneDrive - Lund University/projects/Borg/202501_RuO2test/analysis/v2"
COPY=".agents/skills/databrowser-profiling/project"
DATA="${RUO2_DATA_ROOT:-/Users/davide/Library/CloudStorage/OneDrive-LundUniversity/projects/Borg/202501_RuO2test/electricaldata}"
test -d "$DATA"
diff -qr "$COPY/data" "$USER/data"
diff -qr "$COPY/analysis" "$USER/analysis"
diff -qr "$COPY/plots" "$USER/plots"
```

Also compare the `register_*` kinds in `COPY/definitions.jl` with `USER/browser.jl`. If a difference
could alter the workload, stop and ask whether to sync it or proceed with the known difference.

## Start the persistent session

The benchmark environment enables Revise 3.16 struct revision in `bench/LocalPreferences.toml`.
Keep this block:

```toml
[Revise]
revise_structs = true
```

Restart the existing **bench** jmux session once if it predates that preference; do not restart it
during the iteration loop.

Start the bench session with multiple threads and load the interactive harness. `jmux` launches Julia
through its persistent tmux server, so set the thread environment on that server before restarting
the Julia process:

```bash
tmux -L jmux set-environment -g JULIA_NUM_THREADS auto

jmux --restart --project /Users/davide/code/Julia/DataBrowser/bench \
  'using Revise; Revise.includet("/Users/davide/code/Julia/DataBrowser/.agents/skills/databrowser-profiling/scripts/interactive_profile.jl")'
```

`using DataBrowser` happens after Revise, so package source changes are tracked. There is no hidden
or discarded warmup run: the first run and its runtime compilation are part of the reported result.
`profile_live!` calls `Revise.revise(; throw=true)` before every window. There is no separate Revise
assertion.

Keep the exact thread count matched across comparisons. The harness rejects a one-thread session
instead of emitting misleading throughput or capacity numbers.

## Choose the cache state

The copied workload is deliberately named `DataBrowserProfilingRuO2`, so it uses the ordinary,
isolated cache at:

```text
DEPOT_PATH[1]/databrowser/DataBrowserProfilingRuO2/cache.duckdb
```

Do not redirect `DEPOT_PATH`. `cache_mode=:fresh` passes `rebuild=true`, so DataBrowser deletes only
this project's generated `cache.duckdb` and `.wal` before opening. `cache_mode=:resume` uses
`rebuild=false` and resumes whatever state the preceding run closed with. A short fresh run followed
by a resume run is the direct partial-cache benchmark:

```bash
jmux --project /Users/davide/code/Julia/DataBrowser/bench \
  'DataBrowserInteractiveProfile.profile_scan!(mode=:headless, budget=10, profiler=:none, cache_mode=:fresh)'

jmux --project /Users/davide/code/Julia/DataBrowser/bench \
  'DataBrowserInteractiveProfile.profile_scan!(mode=:headless, budget=45, profiler=:cpu, cache_mode=:resume)'
```

Closing the first workspace drains its cache writes. The resume run's initial stage counts are the
truth about the resulting partial state; report them instead of assuming the ten-second seed stopped
at an exact item boundary. A resume run mutates the cache further, so begin a new comparison pair
with another `cache_mode=:fresh` run.

The real RuO2 `browser.jl` uses DataBrowser's default `background_processing=false`. The profiling
harness matches it. Pass `background_processing=true` only when that alternative is the explicit
question, and never compare it against a default run without calling out the difference.

## Primary live workflow

Open the real GUI workspace/browser once. Use `cache_mode=:resume` to continue the dedicated cache,
or `:fresh` to rebuild it as the session opens:

```bash
jmux --project /Users/davide/code/Julia/DataBrowser/bench \
  'DataBrowserInteractiveProfile.start_live!(mode=:gui, cache_mode=:resume)'
```

The user now uses the app normally. Capture a fixed-duration window while scanning, browsing, or
waiting for an interaction. A live window always lasts for its requested budget, even if the
workspace is idle at the beginning:

```bash
jmux --project /Users/davide/code/Julia/DataBrowser/bench \
  'DataBrowserInteractiveProfile.profile_live!(budget=15, profiler=:cpu)'
```

Use `profiler=:wall` for blocking/idle questions. The human refreshes ProfileView or the pprof web
view from the latest capture in the same session:

```bash
jmux --project /Users/davide/code/Julia/DataBrowser/bench \
  'DataBrowserInteractiveProfile.refresh_profileview!()'

jmux --project /Users/davide/code/Julia/DataBrowser/bench \
  'DataBrowserInteractiveProfile.refresh_pprof!()'
```

The agent does **not** inspect either GUI. It queries the emitted `.pb.gz` directly with the bundled
pprof CLI; see [references/profile-tools.md](references/profile-tools.md). The human and agent are
therefore reading different views of the same captured samples.

Edit one source location, then call `profile_live!` again. It applies pending revisions before the
window. Check the numeric stage deltas and throughput. If a clean rebuild is needed, the user can
click **Rebuild Cache** in the open app, or the agent can call:

```bash
jmux --project /Users/davide/code/Julia/DataBrowser/bench \
  'DataBrowserInteractiveProfile.rebuild_live_cache!()'
```

This keeps the same Julia process, workspace, and browser. At the end:

```bash
jmux --project /Users/davide/code/Julia/DataBrowser/bench \
  'DataBrowserInteractiveProfile.stop_live!()'
```

The live loop is the default for performance work. `profile_scan!` is secondary: use it for an
independent fresh/resume measurement that includes `open_workspace`, then closes its workspace and
browser. Prefer headless CPU profiling for engine-only hotspots; use GUI CPU profiling when GUI work
itself is the question.

## Always report these outputs

Copy the numbers, not a qualitative paraphrase:

- mode, profiler, budget, Julia thread count, `background_processing`, cache mode, and exact cache
  path;
- session `open_workspace` and time to first frame once, then each scan-window length, capture
  overhead, and aggregate runtime compilation;
- start → end counts and rates with their units: found/pending/cached **sources**, interpreted/
  processed/analyzed **items**, and processed/analyzed **collections**; starting counts are mandatory
  for resume runs;
- callback seconds, worker capacity, and the explicitly **unattributed** remainder, which includes
  idle/wait time, runtime work, and uninstrumented work;
- total profile samples;
- raw pprof top table with no filters;
- each filter expression followed by its filtered table;
- `summary.txt`, `timeline.csv`, `timeline.png`, portable raw `.jls`, ProfileView-native `.jlprof`,
  and direct pprof `.pb.gz` paths—all derived from the same capture;
- `WIDEST PROJECT BRANCH:` with that frame's cumulative samples, then `WIDEST SHOWN CHILD:` with
  the child's own cumulative samples; never label a whole call chain with an ancestor's count;
- `NEXT HYPOTHESIS:` followed by one concrete query or change justified by the counts and stack.

Do not hide the raw table. Filtering can improve navigation but cannot be used to make a result
look better. The standard render-wrapper filter uses pprof `hide`, which removes wrapper frames from
the view while preserving descendants. Use `ignore` only when intentionally removing entire
matching samples, and report how that changes the total.

Runtime compilation is part of the result. Report it as aggregate compilation time across Julia
threads; do not discard or automatically repeat a run merely to make that number disappear.

With the real default `background_processing=false`, zero processed/analyzed items is not evidence
of a processing stall: those stages may not have been requested. Use source/interpretation progress
or an explicit app interaction, and state which work was requested.

The raw `.jls` remains because native JLPROF discards some system-specific frame-category
information. `.jlprof` is the ProfileView artifact and `.pb.gz` is the agent's direct pprof input.
There is no `cpu_profile.txt`; a flat text dump loses the call-path structure needed here.

## Revise iteration

After a source edit, issue the next jmux command in the same bench session. `profile_live!` invokes
`Revise.revise(; throw=true)` before capturing and keeps the current workspace/browser open. Do not
restart for method or struct edits while the benchmark environment's struct-revision preference
remains enabled.

Use the profile to form one concrete hypothesis. Make one change, run the same budget/profile, and
refresh the agent's CLI table and the human's GUI views. Once the profile signal improves, run the
same workload with `profiler=:none` to compare useful work without sampling overhead.

Profiled and unprofiled item/s are not comparable. Compare unprofiled runs to unprofiled runs under
matched threads, mode, budget, data, and cache state.

## Fresh-process boundary

The persistent loop is the default. Use `scripts/timed_profile.jl` in a fresh Julia process only
when measuring package startup/precompile behavior or performing the final independent check:

```bash
julia --project=bench --threads=auto \
  .agents/skills/databrowser-profiling/scripts/timed_profile.jl gui \
  --budget=45 --profile=cpu --cache=fresh
```

Use `--cache=resume` for the same dedicated resume path and `--fresh-compile` only when compiled
package-cache startup is the question. Use `--background-processing` only for the explicit
non-default experiment. The fresh-process harness emits the same raw `.jls` profile for the
`.jlprof` and `.pb.gz` outputs. A promising result still ends with the user's real `browser.jl` and
real cache; record the commit and `Threads.nthreads()` with that check.

## Integrity

Do not claim improvement from throttled status updates, dropped work, weakened readiness, hidden
samples, or a warm/partial cache presented as cold. Treat callback-capacity remainder as
unattributed scheduling, GC, I/O, compilation, idle capacity, and unsampled work—not as application
overhead. When a number, GUI view, or source path is unavailable, say exactly what is missing.
