# DataBrowser Roadmap

The ordered plan to take the codebase from `MeasurementBrowser` today — a single package with the GUI
and engine entangled — to the `DataBrowser` family in [../vision.md](../vision.md): a scriptable,
agent-ready scientific toolkit. Stages run roughly in order; each lists what it delivers and when it
is done. Near-term stages are concrete; later ones are directional and get detailed as they approach.

## Stage 1 — Split into the DataBrowser family

Break the monolith into the functional family in [../vision.md](../vision.md) §6 —
`DataBrowserAPI · Sources · Cache · Core · Plots · GUI`, the leaves `DataBrowserProfiling` and
`DataBrowserAnnotations`, and the `DataBrowser` umbrella. The payoff: the engine builds and tests
without the GLMakie/CImGui stack, and every package edge lands where extensions plug in.

Do it by **extracting one real, finally-named package at a time**, bottom-up (leaves first), each
wired back into the shrinking main package through `[sources]` path deps so the whole thing keeps
compiling and testing at every step. Naming is part of extraction, not a later pass: each package
lands as `DataBrowser*` with its real boundary, so every step is reviewable on its own terms — "does
`DataBrowserCache` hold the right things and depend on the right things?" — instead of as a reshuffle
under the old names. There is no separate rename phase: the residual `MeasurementBrowser`, once
everything else has been lifted out, *is* the thin `DataBrowser` umbrella, and renaming it is the last
step.

The order, leaves first, is Profiling → API → Annotations → Sources → Cache → Core → Plots → GUI →
umbrella-rename, with the cross-cutting rules enforced at the step that owns them (payload-agnostic API
when API is cut, cache-on-project when Cache is cut, workspace-as-context when Core is cut). The
file-by-file mapping, per-package deps, and the keep-it-green procedure live in
[package-split.md](package-split.md).

- *Done when:* `DataBrowserCore` loads and tests headless with no GLMakie/CImGui; the GUI loads on top
  via the umbrella; a domain package can register loaders/plots against `DataBrowserAPI` without
  touching engine internals; a cache built in one session is reused by the next; and `bench/` depends
  on Core alone (the `MB_BENCH_ENGINE_ONLY` branch is deleted).

## Stage 2 — Generic plotter and the window registry

Make the app useful on new data before anyone writes a custom visualizer. Build the default
visualizers in `DataBrowserPlots` (1-D signals, multiple traces, 2-D arrays, images, collections,
scalar metadata, nested structures, tables, Makie figures) and stand up the GUI **window registry** —
the public surface external packages add their own windows through.
- *Done when:* common exploration works with no custom draw functions, and a built-in window is
  registered through the same public surface an external package would use.

## Stage 3 — Command unification

Move important GUI and REPL actions onto the shared commands declared in `DataBrowserAPI`. Add command
discovery, help text, and structured inputs/outputs/errors, so every interface (and later an agent)
drives the same operations.
- *Done when:* a workspace can be driven end to end from the REPL through the same commands the GUI
  issues, and the command set is introspectable.

## Stage 4 — Project persistence and cache discipline

Make projects reopenable as *environment + data-source config + saved state*. Track provenance, cache
entries, package versions, and command history well enough for trust and reproducibility; support both
lightweight recovery and richer archival bundles. The source/cache/project ownership model is in
[project-persistence.md](project-persistence.md).
- *Done when:* closing and reopening a project restores its state and reuses its cache, and a project
  can be shared and reopened elsewhere from its pins (or bundle).

## Stage 5 — Python package

Ship `databrowser` on PyPI over JuliaCall: resolve the Julia family, open/attach a workspace, and
expose Pythonic wrappers over the same commands, plus Python-side loaders and processors.
- *Done when:* a Python user can open a project, load and inspect data, and drive the generic plotter
  without writing Julia.

## Stage 6 — Analysis packages

Define the extension-package interface and build one or two serious domain packages (XPS,
ellipsometry) as downstream proof points — each depending on `DataBrowserAPI` (+ `DataBrowserGUI` for
its windows), pinned by a project's environment.
- *Done when:* a domain package delivers a specialized workflow (loaders, visualizers, windows) with
  no changes to the base packages, validating that the base gives enough for free.

## Stage 7 — Agent support

Add an MCP layer over the same commands. Prioritize inspection, command discovery, data summaries,
plot creation, exports, and reproducible execution — no GUI operation required.
- *Done when:* an agent can discover and run the command set to open data, summarize it, produce a
  plot, and export a result.

## Stage 8 — Packaging

Produce a reliable bundled desktop app for non-technical users; evaluate PackageCompiler and app
bundlers, and selective `juliac` use for smaller headless components (which the split already makes
possible, since they depend only on Core or lower).
- *Done when:* a non-technical user can install and launch the app without assembling a Julia
  environment, and pick a preinstalled domain package.

## Stage 9 — Long-term protocol evaluation

Once the command model is stable, evaluate an internal command-oriented protocol / ABI, a
client-server split, a browser frontend, multiplayer collaboration, or external-language clients — per
[../vision.md](../vision.md) §14. Directional; scoped only when the command model has settled.
