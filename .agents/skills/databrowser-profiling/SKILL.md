---
name: databrowser-profiling
description: Agent-executed diagnose → attribute → patch → remeasure loop for DataBrowser RuO2 GUI cold-cache scan performance. Harness is the hot loop; user feel gates ultimate verification.
---

# DataBrowser profiling

**Canonical path:** `.agents/skills/databrowser-profiling/` (`.claude/...` is the same inode — edit once).

**What this skill is:** an agent playbook. The **Runner** executes it. The **parent** edits this file and the harness. The **user** is the authority on feel. Numbers drive the hot loop; feel gates ship.

**Primary product question:** wall time from browser up until the scan is usefully done, and how fast useful work moves — including wasted time *before* scan work starts.

**What the harness actually measures (do not conflate):**

| Block | Contents |
|---|---|
| **Startup** (reported separately) | precompile / `using` / project include / `open_workspace` / `open_browser` |
| **Scan window** (hot-loop regression target) | milestones + throughput for `--budget` seconds *after* `open_workspace` returns |

User-perceived “GUI opened then sat there” includes **both**. Optimize the scan window and pre-scan waste inside it; treat startup as its own hypothesis when those lines dominate.

## Roles

| Role | Does |
|---|---|
| Parent (skill author) | Writes this skill + tiny harness changes. Never confuses itself by “being the Runner.” |
| Runner | Executes phases. One harness run at a time. |
| Critic (optional) | Judges Runner productivity / whether a skill edit helped. **Must not** invoke the harness while the Runner might. |
| User | Feel, ultimate GUI check, answers when escalated. |

## Hard rules

1. **Never run multiple benchmarks in parallel.** Lock: `$TMPDIR/databrowser-profiling.lock`. If present, wait or inspect processes — do not delete casually.
2. **One harness run per hypothesis.** Patch → one remeasure. No stacked unmeasured changes.
3. **Hot-loop cache is disposable** (harness temp depot). Never ask the user to delete `cache.duckdb` during iteration.
4. **Ultimate verification** (real `browser.jl` + real cache + user feel) after **every promising change** — not the hot loop.
5. If you cannot verify something: ask the parent **“should I ask the user?”** — do not invent, do not silently skip.
6. Do not ask the user to stopwatch. Timings come from the harness (or profile artifacts).
7. Prefer ordinary Julia profiling over heavy internal/Perfetto tooling. Add lightweight markers only if external tooling cannot answer. Perfetto is last resort.

## Stop / escalate

- **Continue hot loop** while there is a concrete next hypothesis tied to a harness signal.
- **Ultimate verify** when a remeasure shows a clear win on the targeted milestone/throughput (matched budget/threads/mode) — then ask the user to confirm feel on the real project.
- **Escalate to user** when feel and numbers disagree, data root/project drift is ambiguous, GUI cannot open, or ProfileView cannot display — via “should I ask the user?”
- **Do not declare victory** on harness numbers alone.

---

## Phase 1 — Checked-in project fidelity

Harness loads `project/definitions.jl`, **not** the external `browser.jl`. Wrong project ⇒ wrong optimization.

| | Path |
|---|---|
| Checked-in project | `.agents/skills/databrowser-profiling/project/` |
| Real launch | `/Users/davide/Documents/OneDrive/OneDrive - Lund University/projects/Borg/202501_RuO2test/analysis/v2/browser.jl` |
| Real data root | `/Users/davide/Library/CloudStorage/OneDrive-LundUniversity/projects/Borg/202501_RuO2test/electricaldata` (`RUO2_DATA_ROOT` overrides) |

**Runner must run (or equivalent):**

```bash
REAL="/Users/davide/Documents/OneDrive/OneDrive - Lund University/projects/Borg/202501_RuO2test/analysis/v2"
CHECK=".agents/skills/databrowser-profiling/project"

test -d "$RUO2_DATA_ROOT" -o -d "/Users/davide/Library/CloudStorage/OneDrive-LundUniversity/projects/Borg/202501_RuO2test/electricaldata"
diff -qr "$CHECK/data" "$REAL/data"
diff -qr "$CHECK/analysis" "$REAL/analysis"
diff -qr "$CHECK/plots" "$REAL/plots"
```

Also spot-check that `definitions.jl` still wires the same `register_*` surface the real tree uses (`measurements.jl` / `browser.jl` includes). Count `register_item!` / `register_collection_analysis!` / kind names; record the delta.

**Pass:** data root exists; `data`/`analysis`/`plots` diffs are empty or only explained noise; registration surface matches.

**Fail / unclear:** stop. Ask **“should I ask the user?”** — sync checked-in copy, proceed with documented delta, or ultimate-only.

---

## Phase 2 — Hot-loop baseline (harness, GUI)

**Always pass `gui` explicitly.** Script default if omitted is `headless` — that is a footgun, not the hot-loop default.

```bash
julia --project=bench --threads=auto \
  .agents/skills/databrowser-profiling/scripts/timed_profile.jl gui \
  --budget=45
```

| Flag | Role |
|---|---|
| `gui` | Harness-owned `open_browser`. **Required** for the hot loop. |
| `headless` | Engine-only. Attribution tool when isolating GUI/render. |
| `--budget=S` | Scan sampling window after `open_workspace` (script default 30). Target **~30–60 s of actual scan work**. |
| `--warmup=S` | In-process JIT throwaway (default 20). Not part of the measured scan window. |
| `--fresh` | Wipe compiled `DataBrowser*` — **only** for precompile questions. |
| `--profile=cpu\|wall\|allocs` | One magnifier on the scan window. Never compare profiled throughput to a clean baseline. |

**Budget rule:** if stdout says `budget hit`, or `analysis done` / `processing done` is `not reached` and you still need that stage for the hypothesis → **one** re-run with higher `--budget` (e.g. 60–90) **before** attributing. Do not attribute on an unfinished window.

**Cold cache:** harness `mktempdir()` on `DEPOT_PATH` — never touches `~/.julia/databrowser/RuO2/cache.duckdb`.

**Read:** Startup block, Scan window milestones, throughput halves, callback vs capacity, `Timeline events`, `Artifacts:` dir (`timeline.csv`, optional `*_profile.txt`).

Expected wall clock for one run is long (precompile line + using + warmup + budget + GUI). Plan for that; do not start a second run.

---

## Phase 3 — Attribute (one magnifier at a time)

Use the baseline report first:

1. Re-read milestones + `timeline.csv` — where does work start? which stage stalls? does `busy` clear?
2. `--profile=cpu` → `cpu_profile.txt` (patchable frames). Prefer this over internal traces.
3. `--profile=wall` if CPU sampling misses waits.
4. `--profile=allocs` if allocation pressure is suspected.
5. Same-commit `headless` vs `gui` when GUI/render might steal time.
6. Interactive flame graph for the user (below) when a visual pass helps.
7. Perfetto / `profile_internal` — **last resort** (`docs/profiling.md`). We are evaluating whether this surface can shrink.

Hypothesis resolution can be a milestone shift, a %, a flame frame, or a file/function — whatever answers *this* question.

### Interactive flame graph (show the user)

Default agent artifact: `cpu_profile.txt`.

When the user should **click** a flame graph, lowest friction:

```julia
using Pkg; Pkg.add("ProfileView")   # once in bench env if missing
using ProfileView
# disposable cache: pushfirst!(DEPOT_PATH, mktempdir()) before open_workspace
@profview begin
    ws = open_workspace(PROJECT, DATA_ROOT; metadata_file="device_info.txt", cache=true)
    while workspace_status(ws).busy; sleep(0.25); end
    close_workspace!(ws)
end
```

- **Engine callback costs:** `headless` + ProfileView / `--profile=cpu` is enough.
- **GUI/render hypothesis:** do **not** pretend a headless flame graph answers it — use `gui` harness + `--profile=wall`, and/or ask **“should I ask the user?”** to look at Performance tabs while a harness `gui` run is up.

If ProfileView cannot open here, ask **“should I ask the user?”** — do not jump to Perfetto for a simple flame graph.

---

## Phase 4 — Patch → one remeasure

1. One hypothesis tied to a Phase 3 signal.
2. Small patch.
3. One `gui` harness run, **same budget/threads** unless the hypothesis requires a change (document why).
4. Compare matched milestones / throughput / callback share.
5. Clear win → Phase 5. Else new hypothesis (back to 3). No stacked patches.

---

## Phase 5 — Ultimate verification (every promising change)

Real project the user opens:

1. Fresh Julia; DataBrowser from this repo.
2. Delete real cache (allowed here):

   ```julia
   rm(joinpath(first(DEPOT_PATH), "databrowser", "RuO2", "cache.duckdb"); force=true)
   ```

3. `include` real `browser.jl`.
4. **User confirms feel** (responsiveness, status → Fresh/Loaded, progress clears). Runner records commit + threads; user is the gate.
5. Optional: user opens Performance → Project / Memory.

`open_browser` return ≠ first frame. Do not infer FPS from scan milestones.

---

## Feel vs numbers

When they diverge: ask specific questions anchored to harness lines; check missing signals (GUI cost, JIT-in-window, unattributed %); run **one** targeted magnifier; do not trust either side blindly.

---

## Later (not this loop)

- Warm reopen / ~10% partial cache (same ultimate path, don’t delete cache).
- `bench/scaling.jl`, `bench/realistic_browse.jl` — not the RuO2 GUI hot loop.

## Validity

- Do not “win” by throttling, coalescing, dropping work, or weaker readiness.
- Project callbacks are workload; optimize scheduling/caching/publishing around them.
- Regressions need matched conditions + useful work reported with time.
- Callback % ≠ theoretical minimum; unattributed includes engine, scheduling, JIT, GC, idle, GUI.
