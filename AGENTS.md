# Repository Guidelines

For the full architectural model, when needed, see
[docs/ARCHITECTURE.md](docs/ARCHITECTURE.md).

## Commands

```bash
# Run tests (when validation is needed — skip for doc-only / trivial edits)
julia --project --threads=4 -e 'using Pkg; Pkg.test()'

# Precompile without launching UI 
julia --project -e 'using Pkg; Pkg.precompile()'

# Scaling sweep — compare scaling.csv (slow)
julia --project=bench bench/scaling.jl [n1,n2,...]

# Realistic browse — compare scorecard.csv + profile.json (very slow)
julia --project=bench --threads=auto bench/realistic_browse.jl [scale]
```

Benchmark details: [bench/README.md](bench/README.md).

There is no generic launcher script. Project scripts (outside this repo) call `define_project` and
`register_*`, then `open_workspace(project, root)` and `open_browser(ws)`.

## Purpose

The direction is to turn this into the persistent interactive layer around a Julia data project —
not a domain app and not a generic IDE, but the thing that makes “open the data, browse, plot,
annotate, iterate on code, come back tomorrow” work without re-running scripts or holding the
whole pipeline in your head. The package keeps owning everything stateful and expensive (watching
sources, interpreting once, caching at multiple stages, background work with selection priority,
rendering); project code stays thin script logic that shouldn’t know any of that machinery exists.
The bet is on live work: files and project code can change while views stay attached via selections
or rules, plots can be built while the cache is still filling, and the same operations should
eventually be callable from the REPL as from the GUI so a saved workflow is just a replayable session,
not exported figure code. Generic table/column visualizers and composable figures are the main
product expansion; `register_plot!` per kind is the bootstrap, not the end state. Spatial
navigation, tags, and notes are there to make the tree match how people actually think about
devices on a chip, not to replace the tree. Broader “DataBrowser” generalization (non-DataFrame
payloads, tree as a derived view) is a type-system widening on the same workspace model, not a
pivot. Performance isn’t polish — if browse-while-building and warm reopen aren’t fast, none of the
rest matters. What’s explicitly being left behind: figure-script export, browser-owned project
state, bundled experiment projects in the core package, and compatibility layers that slow down
getting to that live workspace. North star: [docs/plans/workspace-vision.md](docs/plans/workspace-vision.md).

## Target applications

Before pre-made domain projects ship with the app, DataBrowser needs to be installable as a
standalone executable (or equivalent distribution) that a user can launch without hand-assembling a
Julia environment.

The long-term plan is to develop bundled projects on top of this engine — each one a complete
workflow for a technique, not just parsers and plot callbacks. The list below is the product
direction: it shows what the platform must eventually support (multi-window layouts, interactive
fitting, responsive updates while parameters change, and so on). Add to this list when a new
application is scoped.

- **Semiconductor / ferroelectric characterization** — IV, CV, PUND, fatigue, and related
  measurements (this kind of work already lives outside this repo at 
  /Users/davide/Documents/OneDrive/OneDrive - Lund University/projects/Borg/202501_RuO2test/analysis/v2).
- **XPS analysis and fitting** — spectrum import, peak models, constraints, and every window needed
  to fit data interactively, responsively, and flexibly (tools like CasaXPS do this poorly today).
- **Ellipsometry analysis and fitting** — layered optical models, maps vs wavelength, live parameter
  exploration; aim toward CompleteEASE-class capability.
- **Further techniques** — this list will grow.

## Architecture

Project scripts describe how to recognize files, parse them into items, and draw plots. The package
handles directory scanning, background processing, DuckDB caching, the item tree, selection, and the
browser UI. Project code should not touch cache files, background jobs, or UI state.

When a workspace opens, the package scans the data root, finds source files, and interprets each one
into logical items using the project's registered callbacks. That work runs through a dependency
graph with five stages: interpret the source file, process each item, analyze each item, then
process and analyze at the collection level. Completed results are published into an index that
the tree and plots read from. Work is event-driven — background workers finish tasks and publish
updates; the GUI does not poll a job queue. If the user selects items that are still processing,
that work gets higher priority. Full detail: [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md).

To register a new measurement type, see [docs/api.md](docs/api.md). For IDs, collection paths, and
the item tree, see [docs/data-model.md](docs/data-model.md). Recipe `detect` callbacks are tried in
registration order; the first match wins, so register specific filename patterns before general ones.

| Editing… | Look in… | Doc |
|---|---|---|
| `register_*`, item callbacks | `lib/DataBrowserAPI/` | [api.md](docs/api.md) |
| Scanning, item records, hierarchy | `lib/DataBrowserCore/src/ItemIndex/` | [data-model.md](docs/data-model.md) |
| Directory traversal, `metadata.txt` | `lib/DataBrowserSources/` | [storage.md](docs/storage.md) |
| DuckDB cache, writes, reopen | `lib/DataBrowserCore/src/Cache.jl` (+ `lib/DataBrowserCache/`) | [cache.md](docs/cache.md) |
| Background work, loading item data | `lib/DataBrowserCore/src/Workspace/` | [ARCHITECTURE.md](docs/ARCHITECTURE.md) |
| Plot rendering, table viewer model | `lib/DataBrowserPlots/` | [gui.md](docs/gui.md) |
| Browser panels, plot embedding | `lib/DataBrowserGUI/` | [gui.md](docs/gui.md) |
| Profiling, benchmark traces | `lib/DataBrowserProfiling/` | [profiling.md](docs/profiling.md) |
| Tags, notes, spatial layout | `lib/DataBrowserAnnotations/` | [annotations.md](docs/annotations.md) |

`docs/*.md` describes current behavior. `docs/plans/` is for designs not yet built. When you change
behavior that affects the model, update the relevant doc in the same commit — do not copy
architecture into this file.

## Working rules

Julia 1.12; 4-space indent; `snake_case` functions, `UpperCamelCase` types. Don't catch errors that
should be fixed. Docstrings on public APIs. Pre-pre-alpha: replace cleanly, no compatibility shims.
When code and docs disagree, fix the doc in the same commit.

## Testing

When a change needs validation, run the full suite once:
`julia --project --threads=4 -e 'using Pkg; Pkg.test()'`. Skip for doc-only, inspection-only, or
harmless local edits. Fixtures in `test/fixtures/`; inline projects in `test/test_project.jl` and
`test/test_scan_profile.jl`. Plot/GUI tests: metadata, labels, figure creation — not pixels.

## Benchmarks

Use `bench/` for performance work (`julia --project=bench`). Results persist under
`bench/results/` (gitignored). See [bench/README.md](bench/README.md).

- **scaling.jl** — times `status_refresh`, `items_panel`, and `metadata_publish` at increasing item
  counts; writes `scaling.csv` with power-law exponents. Pass smaller size lists while iterating.
- **realistic_browse.jl** — synthetic RuO2-shaped workload: scan while plotting, cache saturation,
  warm reopen. Writes `scorecard.csv`, `profile.json` (Perfetto trace, on by default), and
  supporting CSVs. Use `scale=0.1` to iterate.

Compare runs via `scaling.csv` or `scorecard.csv` + `benchmark.log`. Set `MB_PROFILE_INTERNAL=0`
to skip the structured trace.

## Notes

- The UI needs OpenGL (GLMakie + CImGui).
- Use `test/` fixtures or small subfolders during development — full tree scans are slow.
