---
name: databrowser-profiling
description: >-
  Time and profile DataBrowser cold-cache scan performance on the RuO2 v2
  workload using scripts/timed_profile.jl (GUI or headless). Use when
  investigating scan slowness, pre-scan stalls, GUI vs engine cost, throughput,
  timeline plateaus, or open_workspace/open_browser performance in this repo.
---

# DataBrowser profiling

## Goal

Shorten the time from “browser is up” until the scan is doing useful work, and
raise how fast that work runs. The harness is the everyday measurement tool.
The user’s real `browser.jl` session is the final check that a change actually
helped.

## Workflow

```
- [ ] Diff harness project copy vs user’s RuO2 analysis tree; confirm data root
- [ ] One gui harness run (--budget=45)
- [ ] Read stdout + timeline.png (+ cpu_flame.svg if you profiled)
- [ ] If you need to know which functions are hot: one headless --profile=cpu run
- [ ] One hypothesis → small patch → one matching gui remeasure
- [ ] If the harness looks better: user verifies with real browser.jl
```

Run only one harness at a time. The script takes
`$TMPDIR/databrowser-profiling.lock` and errors if another run holds it.

## Diff before measuring

`timed_profile.jl` loads
`.agents/skills/databrowser-profiling/project/definitions.jl` — a copy of the
registration and plot code from the user’s RuO2 tree. The user’s day-to-day
launcher is `…/analysis/v2/browser.jl`. If the copy is stale, you optimize the
wrong code.

```bash
USER="/Users/davide/Documents/OneDrive/OneDrive - Lund University/projects/Borg/202501_RuO2test/analysis/v2"
COPY=".agents/skills/databrowser-profiling/project"
DATA="${RUO2_DATA_ROOT:-/Users/davide/Library/CloudStorage/OneDrive-LundUniversity/projects/Borg/202501_RuO2test/electricaldata}"
test -d "$DATA"
diff -qr "$COPY/data" "$USER/data"
diff -qr "$COPY/analysis" "$USER/analysis"
diff -qr "$COPY/plots" "$USER/plots"
```

Also compare `register_*` kinds in `COPY/definitions.jl` to `USER/browser.jl`.

If the diff is non-empty and you cannot tell whether it matters, stop and ask
the user whether to sync the copy, proceed anyway, or skip the harness and only
use the real app. Do not guess.

## Run the harness

From the DataBrowser repo root:

```bash
julia --project=bench --threads=auto \
  .agents/skills/databrowser-profiling/scripts/timed_profile.jl gui \
  --budget=45
```

What this does: load DataBrowser from this repo, include the project copy above,
open a workspace on `$DATA` with a **throwaway** cache (never
`~/.julia/databrowser/RuO2/cache.duckdb`), open the GUI, sample the scan for
`--budget` seconds at ~10 Hz, write artifacts under `bench/results/`, then shut
down the browser and workspace.

| Arg | Effect |
|---|---|
| `gui` | Opens the window. Use this for the normal loop. If you omit the mode, the script defaults to `headless`. |
| `headless` | No window. Use for engine CPU profiles and gui-vs-engine comparisons. |
| `--budget=S` | How long to sample the scan after `open_workspace` (default 30). Use 30–60 s of scan time. |
| `--warmup=S` | Throwaway scan first so JIT is not counted in the sample window (default 20). |
| `--fresh` | Wipe compiled `DataBrowser*` package caches before precompile. Only when the question is precompile time. |
| `--profile=cpu` | Sample CPU during the scan window and write `cpu_flame.svg` (flame graph). Prefer **`headless --profile=cpu`**. |
| `--profile=wall` / `--profile=allocs` | Wall-time tree or allocation profile instead. One of these flags per run — not combined. |

Do not compare item/s from a `--profile=*` run to an unprofiled run; sampling changes cost.

If the summary says `budget hit` and you still need a later pipeline stage for
your question, run **once** more with a larger `--budget` before concluding.

## Read the outputs

Stdout ends with `Artifacts: <dir>`.

**Startup block**

- `open_workspace` — workspace construct + cache open + scan kickoff.
- `open_browser` — returned after the render task is started.
- `time to first frame` — seconds from `open_browser` call until the first
  non-blank frame is submitted (startup surface or full UI). **This is the GUI
  startup number to minimize**, not `open_browser` alone.

**Scan window**

- `discovery done` — time in the sample window when `sources_found` first equals
  the final source count for that run (the full file list is known).
- `interpretation done` / `processing done` / `analysis done` — when that stage
  has finished every item, or `not reached` if the budget ended first.
- Throughput — items completed in the window divided by window length.
- Callback vs capacity — time inside project callbacks vs (window × worker
  threads). The remainder is everything else (engine, scheduling, cache, GUI,
  idle). Treat a huge remainder as “need a flame graph,” not as a finished
  diagnosis.

**`timeline.png` / `timeline.csv`**

Line plots of `sources_found`, `cached_sources`, `interpreted`, `processed`,
`analyzed` vs time (~10 Hz). Use this to see dead time (flat) and how fast each
stage advances (slope). This is the main view of scan progress. It does not
name functions — that is what `--profile=cpu` is for.

**`cpu_flame.svg`**

Open in a browser. Wide stacks are where CPU time went. Produce this with
`headless --profile=cpu` when the timeline shows a problem but you need a code
location to patch.

## When you cannot verify something

If a path is missing, a diff is ambiguous, a flame graph will not open here, or
the numbers and the user’s feel disagree — do **not** invent an answer. Ask
whether you should escalate that specific question to the user. Wait for a yes
before pinging them. Phrase the question so they can answer in one or two
sentences (what you saw, what you need them to confirm).

## Profile when the timeline is not enough

The timeline tells you *when* work stalls. It never names *which function* is
hot. When you need that:

```bash
julia --project=bench --threads=auto \
  .agents/skills/databrowser-profiling/scripts/timed_profile.jl headless \
  --budget=45 --profile=cpu
```

Open `cpu_flame.svg` from the artifacts dir.

If you think the **GUI** is slowing the scan (not just the engine), run the
same budget once as `gui` and once as `headless` **without** `--profile`, and
compare `timeline.png` slopes and throughput. Do not use a GUI CPU profile to
judge scan workers — it is dominated by the render loop.

## Patch → remeasure

1. State one hypothesis tied to a concrete signal (e.g. “~11 s after first frame
   with sources known and `interpreted` still 0” or a named frame in the flame
   graph).
2. Make the smallest code change that tests that hypothesis.
3. Run **one** `gui --budget=45` harness again (same threads; change budget only
   if the hypothesis requires it, and say so).
4. Compare the new `timeline.png` and summary to the previous run: shorter dead
   time? steeper `interpreted` slope? earlier processing?

Do not stack several unmeasured patches. A clear harness improvement is the
gate to user verification — not the end of the story.

## User verification (every promising harness win)

The harness uses a temp cache and the project copy. The user cares about the
real app and the real cache path.

1. Fresh Julia process; DataBrowser from this repo.
2. Delete the real cache (only here — not during harness iteration):

   ```julia
   rm(joinpath(first(DEPOT_PATH), "databrowser", "RuO2", "cache.duckdb"); force=true)
   ```

3. `include`  
   `/Users/davide/Documents/OneDrive/OneDrive - Lund University/projects/Borg/202501_RuO2test/analysis/v2/browser.jl`
4. Ask the user whether the pause after the window appears feels shorter and
   whether browsing while the scan runs feels better. Their answer is the gate.
   Record git commit and `Threads.nthreads()` with their answer.

Warm / partially filled cache is later: same steps, skip deleting `cache.duckdb`.

## Validity

Do not make the chip look “Fresh” by throttling updates, dropping work, or
weakening readiness. Project callbacks are the workload; speed up the engine
around them. Claim improvement only under matched conditions, and report useful
work (timeline shape, items done) next to wall time.
