# GUI Internals

> The CImGui frontend. This doc describes GUI state and interaction contracts, not a map of files.

## Entry point and frame loop

`start_browser(root_path; project)` creates one `BrowserState` and runs the render loop. External
project code passes its project object directly; the browser then uses that project for every source
root opened in the window. Without that argument, the browser uses its saved bundled-project
preference. Docking layout is configured once at startup: left side for navigation and information,
right side for plot-oriented work.

## State boundary

Render functions receive one `BrowserState`. Its `workspace` field is the browser's single reference
to the open workspace. `PlotState` owns plot windows and choices, `FigureScriptState` isolates the
interface that Workflow will replace, `TableInspectorState` owns the generic file-inspection window,
and `PerformanceState` owns samples used only for diagnostics.

The workspace owns the source, item index, selected collection and item identities, loaded cache,
loaded item memory, and the state of source-scan and cache work. Those values must not be copied
into browser state.

Source-scan and cache work each have independent state, progress, errors, cancellation, and task
ownership. This prevents one operation from overwriting the status shown for the other.

## Panels

| Panel | Role |
|---|---|
| Hierarchy tree | Multi-select tree, primary navigation. |
| Plot Area | Main plot window with plot-kind chooser, Live toggle, Detach, Export, and Help. |
| Information | Device modal and figure scripts. |
| Table Inspector | Generic delimited-text preview and quick X/Y plot. |

### Selection flow (tree → plot)

1. A tree or future spatial-browser interaction writes collection paths and item keys to
   `workspace.selection`.
2. Each frame resolves those stable ids against `workspace.index` and applies tag visibility.
3. Plot windows with `Live` enabled use that visible selection and redraw when it changes.

A new panel changes selection only through the workspace. It does not maintain its own selected
objects.

## Plot rendering

Plots render directly when their inputs change. `PlotViewState` stores stable item keys,
plot type, figure, errors, and Live setting. It resolves those ids from the workspace index when
rendering, so browser state never keeps a second copy of item records. The main and detached
plots use the same type and the same rendering path.

Each plot window owns its plot kind and Live setting. The main plot starts with Live enabled, so it follows the browser selection. Detached plot windows start with Live disabled, so they keep the items they were created with unless the user enables Live in that window.

While the app runs, plot choices are stored as `item kind => plot type`. Project-local
persistence writes the type names to `measurementbrowser.toml`. When the main Live plot sees one
item kind, it restores that kind's last plot choice.

## Context menus (right-click)

Two existing popups, both attached to selected browser items:

- Tree nodes: "Mark Bad" / "Unmark Bad".
- Item rows: "Open Plot in New Window", mark/unmark bad.

Both apply to multi-selection via `_set_collections_bad!` / `_set_items_bad!`.

## Table Inspector

`inspect_table(path)` opens an arbitrary delimited text file and returns a bounded preview. It keeps
the file's own columns, detected delimiter, optional header row, first data row, approximate row
count, and warning text. It does not rename columns, infer item meaning, call project data
readers, or write cache entries.

The browser exposes this through `Inspect -> Table Inspector`. With `Live` enabled, the window
follows the current item selection and opens that source path. Opening a file manually turns
`Live` off. Its quick plot is only a visual check of two preview columns.

## Annotations in the GUI

Tree nodes and item rows colour their text from the dominant effective tag for that item.
Untagged items render in the default colour. The `bad` tag's catalog entry is added on demand
whenever the user marks something bad. Visibility filtering uses the same effective-tag lookup over
the item key and its ancestors.

This is currently the browser's only annotation interface. Coordinates, spatial layout, and
notes are not yet shown or edited by the GUI.

## MakieIntegration constraint

Makie integration keeps **one Makie screen per `title_id`** in a global context. Multiple
independent canvases require either a unique `title_id` per canvas or a non-Makie renderer for the
secondary surface.
