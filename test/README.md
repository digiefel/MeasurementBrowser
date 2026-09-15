# Tests and performance

Run commands from the repository root. `bench/Manifest.toml` pins dependency versions; Pkg gives
each package test process only that package's dependencies and `Test`.

## Design constraints

- Each package owns its tests and fixtures, using only its declared dependencies and `Test`.
- Tests check public behavior and fundamental invariants. Benchmarks may read internal counters
  for measurements. Missing interfaces belong in concise TODOs, not test-only workarounds.
- Pkg handles dependency resolution and compilation. The runner selects work and records results;
  it must remain simpler than the system being measured.
- Reuse unchanged successful workloads and compiled code during development. Full verification
  measures clean precompilation when package inputs change, then reuses that compilation.
- Measure package overhead with fixed inputs and minimal user callbacks. Include import time,
  process-to-browser readiness, throughput, latency and memory. Faster execution must not fail a test.
- Stop on failure. Publish `status.txt` only after full verification succeeds; test failures do not
  prevent Git commits. This workflow installs no commit hook or GitHub Action.

## Commands

```sh
# Full verification: affected package suites, clean precompilation, engine and real-browser smoke.
julia --project=bench --threads=auto test/runtests.jl

# One package; --force reruns its tests even when the last successful inputs are unchanged.
julia --project=bench --threads=auto test/runtests.jl Core
julia --project=bench --threads=auto test/runtests.jl Core --force

# One test file. Explicit file selections always run and do not certify the entire package.
julia --project=bench --threads=auto test/runtests.jl Core test_workspace.jl

# One benchmark using existing compiled caches; --force reruns it with unchanged inputs.
julia --project=bench --threads=auto test/runtests.jl bench engine
julia --project=bench --threads=auto test/runtests.jl bench browser
julia --project=bench --threads=auto test/runtests.jl bench engine --force

# Measure clean precompilation again, regardless of previous results.
julia --project=bench --threads=auto test/runtests.jl precompile --force
```

`Pkg.test("DataBrowserCore"; julia_args=["--check-bounds=auto"], test_args=["test_workspace.jl"])` in the bench environment is the
standard alternative to the runner. It always runs. The runner uses normal Julia bounds checking
for both tests and benchmarks so Pkg does not compile a second set of package images. Package tests never import the umbrella or
another package's tests. No test-only runtime dependencies are added to make a suite pass.

| Package selector | Test files and behavior |
| --- | --- |
| `API` | `test_identity.jl`: stable collection identity across labels, types and parent paths. |
| `Annotations` | `test_annotations.jl`: coordinates, layouts, tag/note persistence and inheritance. |
| `Sources` | `test_directory.jl`: discovery, exclusions, recursion and change fingerprints. |
| `Recipes` | `test_stages.jl`: registration through typed stages and CSV loading. |
| `Cache` | `test_persistence.jl`: payload replacement/deletion across reopen; cached `nothing`. |
| `Core` | `test_workspace.jl`: reconstruction, metadata precedence, reuse, invalidation and failures. `test_tables.jl`: heterogeneous table values and provenance. |
| `GUI` | `test_clipboard.jl`: clipboard interchange preserves cell boundaries and quoting. |
| `Plots` | `test_view.jl`: live and fixed plot selections survive view persistence. |
| `Profiling` | `test_timings.jl`: concurrent instrumentation, snapshots and lifecycle. |

## What runs again

The runner stores successful input fingerprints in ignored `bench/results/checks.toml`. A package
fingerprint includes its source, tests, fixtures, resolved dependency closure and local dependency
source contents. Julia version, CPU, thread count, project preferences and the execution command
also participate. Runner edits that leave the command unchanged do not invalidate package tests.
Changing a dependency's tests alone does not invalidate its consumers. Failed runs are not recorded.
An unchanged suite is reported as a reused result, not as a newly executed test. `--force` reruns
selected tests and benchmarks; it only forces clean compilation with the `precompile` selector.

Clean compilation runs when package code or its resolved dependencies change. Its dedicated depot
has an empty compiled-package directory; installed package sources and binary artifacts are shared.
Julia's bundled resources remain available. The measured command precompiles all DataBrowser
packages, the umbrella and their dependencies. Subsequent tests and benchmarks use that compiled
cache. Test-only and documentation edits do not trigger clean compilation.

All runner workloads use `bench/results/depot`; Julia builds missing compiled caches there normally.
Package/file tests and benchmark-only runs never request clean compilation. Use those while editing. Run full
verification once when the change is ready. The full command reuses successful unaffected workloads.

## Measurements

A successful full run writes `bench/status.txt`. Partial runs print their results and record them
in ignored `bench/results/checks.toml`; they do not update the snapshot.
The snapshot identifies measurement dates and input fingerprints. Review its diff with the code.
Compare timings on the same machine, Julia version, thread count and workload. Measurements are
regression evidence; noisy timing differences are not automatic test failures.

| Metric | Meaning |
| --- | --- |
| `precompile_s` | Wall time of the clean precompile command, including Julia/Pkg startup. No package compiled images from the user depot are available. |
| `import_s` | `using DataBrowser` in the fresh process for the initial browser opening, with valid compiled caches. |
| `browser_open_s` | Process launch through a presented full browser frame after extension warmup, including import and workspace opening. Preparation screens do not count. |
| `browser_reopen_s` | Same boundary with a saved data cache, again in a fresh process. |
| `frame_ui_ms` | Mean CPU wall time building 60 or more browser frames after the first plot; excludes renderer/vsync time. |
| `index_us_per_item` | Open, index, query IDs and save 10,000 one-row items, including cache close, divided by item count. |
| `write_mib_s` | Process and store 390 payloads, including completed disk writes at close and five concurrent read probes. |
| `read_mib_s` | Materialize all 400 saved payloads after reopen, with the OS file cache warm. |
| `cached_materialize_ms` / `concurrent_materialize_ms` | Materialize ten saved items, idle / submitted alongside bulk materialization. Five fixed read probes run; either reads or writes may finish first. |
| `reopen_ms` | Restore and validate the 400-item workspace index; payload loading is measured separately. |
| `peak_rss_mib` / `browser_peak_rss_mib` | Highest resident-memory sample of engine / browser workers, including Julia and native allocations. Sampling is every 100 ms; shorter peaks may be missed. |

The engine workload uses prepared type-API values: 400 items × 10,000 rows × four Float64 columns.
Throughput uses logical uncompressed payload bytes. User callbacks return prepared values; no
numerical analysis or parsing is profiled. Measurements are medians of three runs after one warmup;
`MB_BENCH_REPEATS` controls repetitions.

The browser smoke project has three items in two collections. It exercises directory discovery,
registered item/collection stages, a real Makie plot, selection, saved-cache reuse, source/metadata
updates and live addition/removal. Expected values and callback counts check correctness. It opens
actual windows and requires a working graphics session; unavailable graphics or renderer errors
fail the workload. Browser startup is measured once per opening and includes script setup; readiness
is observed on the next render-loop iteration after presentation, plus up to 10 ms polling delay.
