---
name: databrowser-profiling
description: Profile and benchmark the DataBrowser application on its real RuO2 v2 project, either headlessly or while its native CImGui/GLMakie GUI runs. Use for DataBrowser throughput, cache-build, workspace, plot, responsiveness, memory, Julia CPU/allocation profiling, GLFW, or GUI performance investigations in this repository.
---

# DataBrowser application profiling

Use this skill for whole-application scenarios and the real RuO2 workload. Use `julia-bench` for
isolated kernels, inference, and BenchmarkTools experiments after this workflow identifies a hot
path. Do not reproduce its general microbenchmark guidance here.

Do not use DataBrowser's internal profiler or profiling UI. They are scheduled for removal. Use
Julia's standard `Profile`, `Profile.Allocs`, timing, allocation, and memory tools so this workflow
remains valid after that removal.

## Start every investigation

1. Read `AGENTS.md` and `bench/README.md`.
2. Record the target checkout, commit, Julia version, thread count, and active package path. Never
   assume the RuO2 v2 manifest points at the checkout under test.
3. State the scenario: cold isolated cache, warm reopen, full background processing, or interactive
   GUI work. Compare only identical scenarios.
4. Ask: **“Do you want the native GUI profiling/analysis pass as well, or headless only?”**
   Headless work may proceed while waiting. Do not open a native window without approval.

Never call throughput linear merely because it is much faster or no longer quadratic. Report raw
times and rates, fit total time against item count, and distinguish `O(N)`, mildly superlinear, and
`O(N²)` behavior.

## Real project

Use:

```text
/Users/davide/Library/CloudStorage/OneDrive-LundUniversity/projects/Borg/202501_RuO2test/analysis/v2
```

The data root is `../../electricaldata`; `browser.jl` defines the project. Its manifest may point at
`/Users/davide/code/Julia/DataBrowser`, so the bundled driver launches with the checkout under test
as the Julia project and loads only the v2 definitions. Confirm its printed `DataBrowser source`
before trusting any result.

## Headless runs

Use a fresh process. The driver reads the actual v2 project and electrical data, warms the callbacks
once, and then measures a clean rebuild in an isolated temporary cache. It never deletes the normal
project cache.

The top-level `DataBrowser` package currently imports its GUI dependencies even in a headless run.
On macOS, a sandboxed process can therefore stall while initializing system input services before
the workspace opens. If `load-only` does not promptly print `DataBrowser source`, rerun the same
headless command outside the sandbox with approval. This grants the process access to the required
macOS services; it does not open a window or authorize GUI actions.

CPU profile:

```bash
julia --project=. --threads=auto \
  .agents/skills/databrowser-profiling/scripts/ruo2_v2_profile.jl headless --profile=cpu
```

Allocation profile in a separate run:

```bash
julia --project=. --threads=auto \
  .agents/skills/databrowser-profiling/scripts/ruo2_v2_profile.jl headless --profile=allocs
```

Use `--profile=none` for an uninstrumented baseline. Add `--background-processing` only when the
scenario is intended to process and analyze the complete workspace. `--timeout=SECONDS`,
`--sample-rate=RATE`, and `--output=PATH` control long runs and artifacts.

The driver writes `run.csv` and either `cpu_profile.txt` or `allocation_profile.txt`. Inspect the
flat report first, identify application frames with high self or total samples, then use
`julia-bench` for focused measurement. Do not infer a cause from allocation totals alone.

For scale claims, run `bench/scaling.jl` at several sizes too. The real v2 run is the acceptance
workload; the synthetic harness is controlled attribution, not a substitute.

## GUI runs

After explicit approval, launch the real native app outside the sandbox in a fresh interactive
Julia process with a PTY:

```bash
julia -i --project=. --threads=auto \
  .agents/skills/databrowser-profiling/scripts/ruo2_v2_profile.jl gui
```

Request escalated execution because GLFW needs the macOS WindowServer. Keep the process session so
commands can be sent to its REPL and it can be stopped cleanly. Process startup is not GUI
validation: confirm with the user that the window is visible and responsive.

The driver defines two REPL helpers that profile all Julia tasks while the GUI remains live:

```julia
profile_gui_cpu(30)                 # writes gui_cpu_profile.txt
profile_gui_allocs(15; rate=0.01)  # writes gui_allocation_profile.txt
```

These helpers only start Julia's profilers and wait for the requested sampling period. They do not
drive, automate, or prescribe GUI actions. The user decides what to do in the application during a
capture. Record that scenario in plain language and keep distinct scenarios in separate captures.
Do not require any DataBrowser window, tab, button, or internal profiling API.

For frame-loop questions, CPU-sample a long enough window to collect repeated frames. For latency
of one known action, first identify the callable boundary from the sample profile, then benchmark or
time that boundary directly outside the render loop where possible. Profile allocations separately
because allocation sampling changes runtime cost.

Close the workspace and Julia process cleanly after collecting artifacts. Never claim visible GUI
success from headless checks or from an invisible GLMakie screen.

## GLFW and GLMakie failures

Distinguish a package-precompile failure from an application crash. A stack ending in
`_glfwGetMonitorPosCocoa` while compiling GLMakie means its precompile workload attempted monitor
access without a usable WindowServer; it does not prove that DataBrowser's GUI crashes.

For a GUI run, retry the launch outside the sandbox after approval. For a headless run, or if only
GLMakie's precompile workload fails, bypass package-image generation and warm the real workload
inside the process:

```bash
julia --compiled-modules=no --project=. --threads=auto \
  .agents/skills/databrowser-profiling/scripts/ruo2_v2_profile.jl headless --profile=cpu
```

For an approved GUI fallback, use the same `--compiled-modules=no` flag with the interactive GUI
command. Exclude startup and warm-up from measurements. If an unsandboxed visible launch still
fails, preserve the exact stack, Julia/GLMakie/GLFW versions, and command before diagnosing a real
GUI startup problem.

## Report results

Lead with:

- exact checkout/commit and confirmed `pathof(DataBrowser)`;
- scenario, cache mode, data root, threads, profiler, and whether the GUI was visibly tested;
- raw total time, items/s, allocations or RSS, and interaction latency where relevant;
- fitted scaling exponent with all raw points;
- dominant Julia CPU or allocation frames, separated into application and dependency code;
- limitations, including any GUI pass not approved or not visibly confirmed.

Apply one targeted fix, rerun the identical headless scenario, and repeat an approved GUI capture
when the defect is interactive.
