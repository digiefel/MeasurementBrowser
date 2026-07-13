# DataBrowser Roadmap

The ordered plan to take the codebase from a single monolithic package — GUI and engine
entangled — to the `DataBrowser` family in [../vision.md](../vision.md): a scriptable,
agent-ready scientific toolkit. Stages run roughly in order; each lists what it delivers and when it
is done. Near-term stages are concrete; later ones are directional and get detailed as they approach.

## Stage 1 — Public foundation

The functional package family is in place: `DataBrowserAPI · Sources · Cache · Core · Plots · GUI`,
the Profiling and Annotations leaves, and the `DataBrowser` umbrella. Core loads headlessly without
the GUI or Makie stack, Cache is below Core, and Plots extends the GUI host.

The public data foundation has two paths described in [../api.md](../api.md): ordinary-data
registration pipelines and concrete type-based sources/items. Both paths share one workspace engine,
identity model, metadata pipeline, and cache.

- *Done when:* the four projects under [`examples/`](../../examples/README.md) execute through the
  documented public API; custom item types remain intact through processing and visualization;
  registration replacement works under Revise; and no public callback receives internal index,
  cache, scheduler, or browser values.

## Stage 2 — Generic plotter and the window registry

Make the app useful on new data before anyone writes a custom visualizer. Keep lightweight
visualizers for Core-supported data shapes in `DataBrowserGUI`, then build the Makie-backed defaults
in `DataBrowserPlots` (1-D signals, multiple traces, 2-D arrays, images, collections, scalar metadata,
nested structures, enhanced table plots, Makie figures) as a default GUI extension. Stand up the GUI
**window registry** — the public surface external packages add their own windows through.

Add a global pipeline-stage control for inspecting how the selected data changes through the project.
It is a view control: moving it changes what the browser displays without rewinding or invalidating
the workspace. The minimum path is **Source → Loaded → Items**. Optional positions appear when the
project defines them: **Processed → Analyzed → Collection processed → Collection analyzed**. Missing
positions are skipped, giving three to seven meaningful stops rather than exposing internal work
states. The Items position also shows the results of `label`, `collection`, and `id`; those cheap
descriptions are not separate stops.

Moving backwards coalesces sibling items into their shared loaded value and source. Moving into
collection positions groups selected items by collection. Analysis positions show the metadata added
at that point alongside the unchanged data. Intermediate values use generic inspectors when a custom
visualizer only supports the final form. Results come from the existing cache where available;
`read` and `entries` become separately observable internal results, with bounded in-memory retention
for values that are not persisted.

- *Done when:* common exploration works with no custom draw functions, DBPlots is registered through
  the same public surface an external package would use, `DataBrowserGUI` can load without GLMakie,
  and moving the pipeline-stage control gives responsive access to every meaningful result available
  for the current selection without rerunning work that is still cached.

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
ellipsometry) as downstream proof points — each depending on `DataBrowserAPI`, plus
`DataBrowserGUI` for windows and `DataBrowserPlots` for Makie visualizers, pinned by a project's
environment.
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
