# Vision: Agent-Ready Scientific Data Toolkit

This is the product and architecture north star for the `DataBrowser` family.
It defines vocabulary and long-term shape.

## 1. Purpose

DataBrowser is a Julia-native toolkit for scientific data processing, visualization, inspection, and
interactive analysis. It is for domain experts who use software as an instrument, not primarily for
software engineers. It should make messy, heterogeneous scientific data easier to load, inspect,
process, visualize, cache, automate, and package into reusable analysis workflows.

The near-term product is a Julia package family with a GUI, REPL access, GLMakie-based
visualization, CImGui-based interaction, DuckDB-backed caching where appropriate, and a registration
API for user-defined analysis code. The long-term direction is an agent-ready, scriptable, extensible
environment supporting specialized scientific analysis packages, multiple interfaces, and eventually
a more formal internal protocol.

The core design goal is to preserve the freedom of arbitrary scientific code while progressively
adding structure around the parts that benefit from it: workspace state, datasets, visualizers,
pipelines, cache metadata, inspection, provenance, command execution, and packaging.

## 2. Product Vision

The toolkit should sit between a programming environment and an exploratory scientific application,
closer to a powerful scientific instrument than to a data-science dashboard.

A user should start with a directory of files, quickly inspect what is inside, load the data with
minimal ceremony, visualize common structures immediately, and gradually add specialized analysis
logic. A new user should get useful plots without writing custom setup and draw functions. An
advanced user should be able to define custom loaders, processors, visualizers, UI panels, and
complete domain-specific analysis packages.

The system should support both exploratory and repeatable work. A user should be able to click around
interactively, but also replay, script, inspect, export, and automate what happened. This matters for
scientific trust, reproducibility, and agentic use.

The long-term product should make it natural to distribute open, modifiable analysis packages. These
packages may include loaders, processors, visualizers, custom UI components, examples, documentation,
and templates. Commercial value can come from high-quality packaged analysis workflows while
preserving the user's freedom to inspect and modify the underlying code.

## 3. Scope and Strategy

The near-term strategy keeps the system Julia-native, flexible, and usable while avoiding
architectural commitments that would slow down the core product. Julia is the primary implementation
and extension language for now, because that is where the GUI, Makie integration, package system, and
user extension model already work.

This does not rule out broader language-neutral architecture later. Language neutrality should emerge
from a stable workspace model, command model, and data-handling system rather than being imposed
before the product shape is clear.

The system should not force all data into a single canonical representation. Tables are useful, but
arrays, traces, spectra, images, nested structures, file-backed data, domain-specific records, notes,
metadata, and opaque user-defined values must remain first-class. DataFrames are supported but not
privileged as the foundation.

The system should preserve arbitrary user code as a core capability. Structure is added progressively
where it helps with inspection, visualization, caching, reproducibility, automation, and packaging.

Low-level protocols, browser frontends, multiplayer sessions, alternative backends, and aggressive
binary trimming are long-term options, to be evaluated once the workspace model, command layer,
generic plotter, and extension-package model are mature. The near-term priority is a strong
Julia-centered foundation that leaves those directions open.

## 4. Core Principles

- **Julia-native first.** The most complete experience exists in Julia, where the GUI, Makie
  integration, package system, and extension model are most natural.
- **Arbitrary user code is supported.** Scientific work involves custom data structures, irregular
  files, special corrections, device conventions, notes, and one-off processing. The app provides
  inspection, caching, visualization, and structure *around* arbitrary code rather than forbidding
  it.
- **Progressive structure.** Users begin with callbacks and arbitrary data. When data is recognized
  as arrays, tables, spectra, images, signals, parameter sweeps, nested records, Makie figures, or
  file collections, the app offers richer generic tools. Unknown data stays inspectable.
- **Visualization is a capability layered over data, not a data model.** A visualizer may accept an
  array, nested structure, table, image, spectrum, collection of traces, graph, file-backed data,
  simulation result, or domain-specific structure.
- **One command model behind every interface (GUI/API parity).** GUI actions, REPL commands, CLI
  commands, Python calls, and agent/MCP tools operate on the same underlying workspace operations.
  That shared set of operations is the API (see §5, §9). The GUI is an interface to the workspace
  model, not a separate implementation.
- **Non-blocking REPL.** Starting the browser must not consume the Julia session. The user keeps the
  same REPL to inspect data, run helpers, change project code, and drive the app.
- **Live workspace.** Source files, project code, selections, and views can update without restarting
  the app. Long work has visible progress and does not block unrelated interaction.
- **Execution is separate from inspection.** Users run arbitrary analysis code, but the app makes it
  easy to stop, inspect, visualize, cache, compare, replay, and export intermediate and final data.
- **Generic tools first.** Built-in visualizers cover common work (tables, X-vs-Y plots, overlays,
  heatmaps, histograms, summaries, fit inspection); project code *extends* those tools rather than
  replacing the browsing and composition system.
- **Persistent work, not exported scripts.** Figure-script export is a transitional path. The durable
  artifact is a saved workflow: a replayable, editable, inspectable set of actions drivable from GUI
  or API.
- **Composable figures.** Figures are editable objects in the app — some backed by live data
  selections, others self-contained without a source root.
- **Coarse-grained interoperability.** Foreign-language APIs exchange whole data structures, files,
  arrays, specifications, handles, or results rather than tiny cross-language calls in loops.

## 5. Conceptual Model

The model rests on a small set of named concepts that apply equally to GUI and API use. Keeping them
distinct is what makes caching, reproducibility, and packaging coherent.

**Project.** The durable thing the user means by "this analysis": data plus the code that interprets
it. In normal use it is represented by one editable project file. That file records the project entry
code, data-source configuration, and generated-state locations such as the cache path. It can be
opened, closed, edited, versioned, and reopened.

**Workspace.** One live opening of a project. It holds the current index, loaded data, cache
connection, background jobs, selections, generated figures, errors, logs, and attached interfaces. It
is transient. One project can back one-or-many workspaces, successive or concurrent. Closing the last
workspace does not delete the project file, data-source config, or cache.

**Data source.** A configured origin the project opens: a directory, database query, instrument run,
stream, or remote dataset. Data-source configuration is core project state. Editing it changes the
project file, but it does not rewrite the Julia script that defines callbacks and plots.

**Package.** A distributable extension — reusable behavior with *no data attached* (loaders,
processors, data types, visualizers, commands, UI windows, examples, docs, templates). A package is
not a project: many projects share one installed package by pinning its version. Domain analysis
packages (§11) are the main example.

**Data / item.** Any value the workspace knows about that a user can inspect, manipulate, or
visualize: arrays, tables, nested structures, domain-specific measurements, file references,
database-backed results, Makie objects, image-like or signal-like data, or opaque user-defined
values. A logical browsable entry discovered from sources is an **item** (e.g. a single measurement),
with stable identity, parameters, stats, and access to its data. The system need not understand every
value fully, but it should identify it, track it in the index, summarize it, and expose actions.
A **Core-supported data shape** is a type the toolkit understands out of the box: caching, summaries,
and default visualizers work with no user code (tables today; vectors, arrays, and image-like values
later). The shape *is* the payload's Julia type — there is no separate detection or declaration step,
and shapes are not exclusive: every visualizer whose accepted types match the payload is offered, so
the same matrix opens as a grid, an image, or a heatmap depending on which visualizer the user opens.
The app is type-agnostic in the failure direction too: a payload nothing matches simply is not cached
and gets only the basic inspector, and a visualizer applied to data it cannot handle fails cleanly on
that view — the user can then register their own support. Shapes whose support needs heavy
dependencies (audio codecs, for example) belong in extension packages even when generic, on the same
host/extension pattern as the GUI itself.

**Selection.** A concrete list of items, or a *rule* that resolves to items as the project changes.
Concrete selections reproduce exactly; live rules keep views attached as items are added, removed, or
changed.

**Command.** A stable operation on a workspace — the foundation for GUI actions, REPL, CLI, Python,
testing, and agent/MCP support. A command has a name, inputs, outputs, errors, help text, and
predictable effects, and is inspectable so an agent can discover what is possible.

**Visualizer, Figure, Workflow, Annotation.** A **visualizer** renders or interacts with a kind of
data (generic built-ins or package-supplied), optionally exposing declarative state. At runtime a
visualizer presents as zero-or-more windows the user can spawn, toggle, minimize, and close; the
current plot window is the first such visualizer (today 1-to-many, since it is the only one).
Registered plot APIs like `register_plot!` are how *one* visualizer accepts project callbacks, not
the contract every visualizer follows. A **figure** is
a composed visual artifact: visualizer outputs, layout, styling, annotations, and optional embedded
data. A **workflow** is the persisted sequence/graph of actions (open, select, transform, visualize,
fit, annotate, export) — intent, not pixels. An **annotation** is user-authored metadata on an item
(tags, notes, coordinates, spatial position); a **figure annotation** is a visual mark on a figure
(arrow, region, fit label, text). A tag describes an item; an arrow describes a figure. These are
detailed in §8.

A note on the registration API: `define_project` / `register_*` is one convenient way to describe a
data pipeline with ordinary values and functions. Its reader is the only required callback; identity,
labels, collection placement, metadata, processing, and analysis have useful defaults. The type API
is the fundamental layer for packages that model sources and items directly with multiple dispatch.
Both paths drive the same workspace and commands.

### Scriptable shape (illustrative)

The API should feel like normal Julia composition, not an app-control protocol:

```julia
config    = ProjectConfig(sources=[DirectorySource("data")])
project   = define_project("ruo2.toml"; config)
workspace = open_workspace(project)           # one live opening
browser   = open_browser(workspace; wait=false)   # non-blocking; REPL stays free

items = query_items(workspace, "polarity = 'up'")   # concrete list, or a live rule
select_items!(workspace, items)

fig = Figure(workspace)
plot!(fig, items, X(:voltage_V), Y(:polarization_uCcm2); visualizer=LinePlot)
fit = linear_fit!(fig, items, X(:voltage_V), Y(:current_A))
annotate!(fig, Arrow(fit; label="linear region"))

save_workflow("pund_review.dbflow", fig)
```

The GUI issues the same kinds of actions when the user clicks, drags, or edits plot options. The saved
workflow stores those actions in package-owned terms, not private widget state.

A script that calls `define_project(project_file; config)` creates or reopens the project file, then
registers callbacks and plots. Rerunning the script points at the same project file. If the file
already exists and the passed config disagrees with it, `define_project` should fail unless the caller
explicitly asks to overwrite.

## 6. Near-Term Architecture: the DataBrowser package family

The package is split along **function**, because function and dependencies coincide — the plotter
needs GLMakie, the cache needs DuckDB, the loaders need file libraries — so splitting by what code
*does* also isolates heavy dependencies (fast recompiles) and puts package edges exactly where
extensions plug in. The organizing idea: **`DataBrowserAPI` is what every interface and extension
talks to; the engine packages implement it; the frontends call it.**

```
DataBrowserAPI     the shared surface every interface calls and every extension implements:
                   command declarations (open_workspace, scan, load, select, materialize, plot,
                   export, query, register_*, …), the extension definitions (what a source / loader /
                   visualizer / command is), and the Project value type. Payload-agnostic; tiny deps.
DataBrowserProfiling  instrumentation / traces (leaf: traces don't reference item identity)
DataBrowserAnnotations → API : tags / notes / spatial-layout model, attached to item identities
DataBrowserSources → API : file discovery, loading, metadata extraction, data summaries
DataBrowserCache   → API : persistence (DuckDB where it fits); keyed by PROJECT, not workspace
DataBrowserCore    → API, Sources, Cache, Annotations : the workspace, the index, background work,
                   command execution, project save/load, provenance, REPL, and the runtime that
                   drives a project's callbacks through interpret → process → analyze; base support
                   for common data shapes such as tables, vectors, arrays, and image-like data; shapes
                   needing heavy dependencies (audio, say) live in extension packages
DataBrowserGUI     → API, Core (+ CImGui/GLFW) : the host shell, panels, browser, shared widgets,
                   the window registry, and lightweight visualizers for Core-supported data shapes
DataBrowserPlots   → API, Core, GUI (+ GLMakie) : a first-party/default GUI extension for Makie
                   embedding, high-level plot APIs, and enhanced default visualizers
DataBrowserCLI     → API, Core : a command-line frontend over the same commands (may come later)
DataBrowser        → GUI, default extensions : umbrella / install target; re-exports the API
                   and wires defaults
```

The dependency graph is acyclic with `DataBrowserAPI` at the bottom and the GUI host/extensions at
the top. Heavy dependencies live where their feature lives — DuckDB in `DataBrowserCache`, CImGui/GLFW
in `DataBrowserGUI`, GLMakie in `DataBrowserPlots` — so headless engine work never recompiles the GUI
or GPU stack, and the GUI shell can exist without importing Makie.

`DataBrowserAPI` stays a clean leaf only under one discipline: something belongs in it only if it is
shared by two or more packages **and** is declaration or data — abstract extension types, lightweight
value types, generic-function declarations — with no behavior and no heavy dependencies. Method
implementations and live machinery (the `Workspace`, the index, the scheduler) live in the engine
packages that add methods to those declarations. If the API starts absorbing behavior or heavy types,
that is a signal a boundary is wrong, not that the rule should bend.

See [plans/roadmap.md](plans/roadmap.md) for the concrete boundaries, current-code mapping, and
phasing.

## 7. Registration and Callback Model

The registration API defines named or default pipelines over ordinary Julia data. `read` performs
the per-source-item interpretation once and returns either one data value, one `DataItem`, or
several data items. `DataItem` construction supplies cheap sibling identity, labels, collection placement,
and metadata at the same time data is split. Deferred `process` and `analyze` callbacks receive data
and effective metadata. A registration name provides stable replacement and associates extension or
collection operations. It is not stored as a property of the user's data.

The type API models `AbstractDataSource`, `AbstractDataSourceItem`, and `AbstractDataItem` values
directly. Concrete item types survive materialization, processing, and visualization. Their behavior
comes from multiple dispatch, and item methods have defaults for source-derived identity, labels,
collections, empty metadata, identity processing, and source-backed storage.

Project code stays close to natural scripts: identify data, load it, compute results, and define
domain views. The workspace owns cache files, UI state, background jobs, and package lifecycle.

Callbacks need not be serializable. The system records what it can: callback name, package origin,
signature or declared capability, input data, parameters, outputs, logs, errors, cache keys, and
package version.

Registration remains simple: one reader is a complete pipeline. Additional callbacks express only
the stages the data requires. Packages can instead expose their own constructors and project APIs over
the type layer, including domain-specific configuration and concrete item types.

Python callbacks are supported where realistic, especially for loading and processing. Python
visualizers may return supported data structures, plot specifications, or Makie-compatible objects via
JuliaCall. Full cross-language registration of new Makie recipes is not a near-term requirement.
Most of this section is up to future investigation and direction.

## 8. Generic Plotter, Figures, and Workflows

The generic Makie plotter (`DataBrowserPlots`) is a major near-term priority and a default GUI
extension: it makes the app useful before users write custom visualizers, while proving the same
window registry external packages will use. The base GUI still owns lightweight inspectors for
Core-supported data shapes, such as tables. The Makie plotter should be data-oriented in the
scientific sense: what kind of data is this, what can be inspected, what axes or dimensions exist,
what metadata exists, what natural views are available, what comparisons are meaningful?

The first Makie visualizers should cover: X-vs-Y scatter and line plots; overlays collected by item,
device, tag, or parameter; enhanced table plots; two-dimensional arrays and image-like data; heatmaps
for gridded or pivoted tables; simple summaries and histograms; scalar metadata; nested structures;
fit views for common models (starting with linear fits); and Makie figures — good enough that common
exploration works without custom draw functions. Visualizers consume data through package data-access
functions; they do not read source files, cache files, or GUI state directly.

Progressive enhancement: unknown data gets a basic inspector; recognized data gets generic
visualizers; domain-specific data gets specialized visualizers; packaged modules add the
highest-quality experience. Project-specific visualizers look like extensions of the same system —
they declare what data they need and draw into a figure or panel, and the package owns selection,
composition, persistence, and background work.

**Editable state.** A visualizer exposes editable state — selected x/y data, visible traces,
normalization, offsets, filters, colormaps, axis scales, labels, annotations, linked selections,
export options — which should be inspectable and eventually scriptable.

**Figures, live and self-contained.** A *live* figure references a workspace and one or more
selections and updates when matching items, data, or project code change. A *self-contained* figure
stores enough layout, plotted values, style, and annotations to reopen without source files, and stays
editable in the app. **Figure annotations** are visual marks on a figure — linear fits and fit labels,
arrows, regions of interest, text callouts attached to points or axes — and are distinct from item
annotations (tags/notes on a device or measurement).

**Workflows.** A workflow captures intent rather than pixels: open project → select by rule → load or
process → create visualizer → map columns to axes → set collection axis and style → add fit → annotate
→ export. The first requirement is a stable action model both GUI and API use; a text DSL can come
later as a representation of that model. Workflows support both concrete selections and live rules.
A workflow records the whole session surface: if three visualizer windows are open, each with its own
selection and state, replaying the workflow reopens all three.

Open design questions: whether a self-contained figure stores a data snapshot, only rendered marks, or
both; how live figures degrade when the source project is unavailable; and whether figure annotations
live in workflow files, figure files, or a shared store. See
[plans/plotting-api-design.md](plans/plotting-api-design.md) and
[plans/spatial-browser.md](plans/spatial-browser.md).

## 9. REPL, CLI, and Agent Readiness

The REPL and CLI are first-class interfaces to the same app, not separate utilities. The key
abstraction is the command model.

Every important user action has a command form declared in `DataBrowserAPI`. The GUI calls commands
internally; the REPL (in `DataBrowserCore`) exposes them as Julia functions; the CLI exposes them as
shell commands; an MCP server exposes them as agent tools. Because they are all the same commands,
the interfaces stay consistent by construction.

Commands should exist for opening projects, scanning files, loading data, listing data, inspecting
summaries, running processors, creating visualizers, modifying visualizer state, exporting figures,
exporting data, querying cache state, and showing logs. Commands are inspectable so an agent can
discover what is possible. Agent readiness does not require the agent to operate the GUI: the GUI is
for humans; the command/workspace model is for scripts and agents.

## 10. Python Packaging

The Python module (`databrowser` on PyPI) uses JuliaCall rather than a C protocol. This preserves the
full Julia runtime and gives Python users a high-level experience without prematurely freezing the
internal API.

It should be pip-installable, provide a Python entry point, resolve the Julia package family, start or
attach to a workspace, and expose Pythonic wrappers over the same commands. Python users should be
able to provide custom data sources and processing callbacks. They do not need full support for
registering new Makie recipes in the near term. The Python API should not pretend to be a native
Python rewrite — it is a Python interface to a Julia system, explicit in implementation, smooth for
the user.

## 11. Analysis Packages

Analysis packages are a central product concept. An analysis package is a Julia package that extends
the toolkit for a specific domain, instrument, file format, or workflow. It may include loaders, data
types, processing pipelines, visualizers, UI windows, templates, examples, and documentation. It
depends on `DataBrowserAPI` (to register behavior), on `DataBrowserGUI` if it adds windows, and on
`DataBrowserPlots` if it adds Makie-based visualizers.

The base app makes these packages valuable by providing common infrastructure for file browsing,
caching, data inspection, plotting, UI layout, export, command execution, logging, workspace state,
and agent access. A good analysis package feels like a specialized app built on the general toolkit;
users can install, inspect, modify, and compose it with other packages.

An analysis package is **not** a project: it is reusable, data-free behavior. A user does XPS work by
creating a *project* (their data + state) whose environment *pins* the XPS package version. This is
the likely commercial direction — high-quality, open, modifiable domain packages that save time while
keeping trust high.

## 12. Caching and Persistence

**The cache belongs to the project, not the workspace.** Its identity is a function of the project's
definitions + the data + the parameters, never a workspace instance and never merely "what the last
scan saw." That is what lets it survive closing and reopening, and lets GUI, REPL, Python, and agent
interfaces share the same work. Because a project's code can change while a workspace is live, the
project definition carries a version/fingerprint, and cache identity = project-definition fingerprint
+ data identity + parameters.

DuckDB is used where it fits — metadata, indexes, tabular summaries, cache catalogs, provenance
records, queryable results — without forcing all scientific data into tables. It should be
the primary storage and caching backend for the system. It should handle metadata, indexes,
tabular summaries, cache catalogs, provenance records, and most cached results in a unified way.
The cache should remain conceptually simple and centralized. While it should be possible to store
certain data outside of DuckDB when necessary (for example, large binary data or specialized formats),
this should be the exception rather than the norm. Additional storage backends should be introduced
deliberately and sparingly, and primarily as internal extensions rather than a broad collection
of overlapping formats.

**State ownership.** A clear owner for each kind of state is what makes workflows and a non-blocking
REPL work. `ui_state` is a GUI rendering detail; it must not become the durable owner of workspace
data, loaded data, workflow actions, or figure contents. The workspace object is the package-owned
thing the GUI observes and mutates through the same operations available to Julia callers:

| State | Owner | Persistence |
|---|---|---|
| Workspace state | Workspace object | Recomputed from source, cache, and source metadata. |
| Data cache | Cache layer | DuckDB (+ side files) outside the source root; keyed by project. |
| GUI state | Browser window | Mostly transient; per-dataset state in `databrowser.toml`. App-wide prefs TBD with project-file persistence. |
| Item annotations | Source root | Hand-editable tags, notes, coordinates, spatial positions. |
| Workflow state | Workflow model | Saved workflow file, replayable from GUI or API. |
| Figure state | Figure / workflow model | Produced from workflows or exported separately. |

Normal project persistence is a small text project file plus generated cache. Export/archive
workflows may later package a Julia environment, data copies, and cache, but that is separate from
ordinary open/save. See [cache.md](cache.md) for the current cache design and
[plans/project-persistence.md](plans/project-persistence.md) for the saved-project model.

## 13. Packaging the Julia Application

Near-term packaging prioritizes reliability over minimal binary size. The main target is still a Julia
package family for Julia users; a secondary target is a bundled desktop app for non-technical users,
likely via PackageCompiler or related app-bundling tooling.

The full GUI app should not be designed around `juliac --trim` yet: GUI + GLMakie + runtime callbacks
+ dynamic packages + user extension code is not a good first trimming target. But the codebase should
be structured so smaller headless components (command-line workers, fixed import/export tools,
specific analysis kernels) may become trim-friendly later — which the functional split already
supports, since those components would depend only on `DataBrowserCore` (or lower), not the GUI.

## 14. Long-Term Architecture Options

- **Internal protocol / ABI.** Later, but internal and command-oriented: it exposes workspaces,
  commands, data handles, buffers, errors, and metadata — not arbitrary Julia internals.
- **Client-server.** Realistic once the command model is stable. The server owns the Julia workspace;
  clients could be the desktop GUI, a browser GUI, Python, MATLAB, scripts, agents, or multiplayer
  frontends.
- **Browser frontend.** Plausible once the backend is command-driven; it should not require rewriting
  the scientific core.
- **Multiplayer collaboration.** A long-term extension of client-server: workspace synchronization,
  permissions, conflict handling, shared views, careful state design.
- **Rust.** Not a near-term plan. Possibly useful later for low-level components, packaging helpers,
  file watchers, binary protocols, or performance kernels. A full backend rewrite only makes sense if
  the Julia runtime became the main blocker.
- **GPU acceleration.** Operator-based: built-in operations gain CPU, threaded, distributed, and GPU
  implementations where appropriate. Could be interesting and high-risk/high-reward work.
- **Editable figure bundles.** A strong long-term direction: a saved figure should contain more than
  pixels — figure state, data snapshots or enough data access information, visualizer identity,
  package versions, annotations, layout, and enough provenance to reopen and modify it.

## 15. Roadmap

The concrete, current split is [plans/roadmap.md](plans/roadmap.md). The longer product arc:

1. **Architectural cleanup.** Split the package into the functional family of §6; pin the workspace /
   project / package / command vocabulary; establish `DataBrowserAPI` as the command-and-definitions
   surface.
2. **Generic plotter.** Strong default visualizers for common data shapes so the app is useful without
   custom draw functions.
3. **Command unification.** Move important GUI and REPL actions onto the shared command layer; add
   command discovery, help text, structured inputs/outputs, and errors.
4. **Project persistence and cache discipline.** Make projects reopenable; track provenance, cache
   entries, package versions, and command history well enough for trust and reproducibility.
5. **Python package.** Expose the workspace model and major commands via JuliaCall; support
   Python-side loaders and processors; keep the interface high-level and coarse-grained.
6. **Analysis packages.** Define the extension-package interface and build one or two serious domain
   packages (XPS, ellipsometry) as proof points.
7. **Agent support.** Add an MCP layer over the same commands; prioritize inspection, command
   discovery, data summaries, plot creation, exports, and reproducible execution.
8. **Packaging.** Produce a reliable bundled GUI app; evaluate PackageCompiler, app bundlers, and
   selective `juliac` use for smaller headless components.
9. **Long-term protocol evaluation.** Once the command model is stable, evaluate an internal protocol,
   client-server mode, browser frontend, multiplayer GUI, or external-language clients.

## 16. Key Design Risks

- **Too generic to feel useful.** Countered by making the generic plotter genuinely good and building
  real analysis packages early.
- **Arbitrary callbacks defeating inspection/caching/replay.** Countered with progressive metadata,
  command logging, data summaries, and optional structure around callbacks.
- **Over-investing in language parity too early.** Julia stays the full-power extension language
  near-term; Python gets a high-quality wrapper; other languages wait for the command model and
  possible protocol to settle.
- **Chasing small binaries before the product is stable.** Reliable distribution comes first;
  aggressive trimming is evaluated later for smaller headless components.
- **Exposing too many software-engineering concepts to domain experts.** Countered by keeping the
  default experience instrument-like — the workspace, data, plots, and commands are the surface the
  user sees; packages, environments, fingerprints, and the command protocol stay behind it until a
  user chooses to go deeper.
