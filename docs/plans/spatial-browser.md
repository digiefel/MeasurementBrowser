# Spatial Browser + Annotation System — High-Level Plan

## Context

`MeasurementBrowser` currently navigates devices through a tree panel (`Gui.jl:2542 _render_hierarchy_tree_panel`). It works, but doesn't reflect the **spatial** reality of measurements on a chip — where devices physically sit, which sites are clustered, which chips are in flight together. The user wants a hybrid file-browser + GDS-style canvas where:

- Positioned items (devices, sites, sub-sites) render at hierarchical XY coordinates from `devices_info.txt`.
- Unpositioned items (chips) get an initial grid layout but are then freely draggable like icons in a macOS file browser; a "reset layout" button reapplies the grid on demand.
- A "level" indicator + `+/-` controls let the user descend/ascend the hierarchy, with deeper levels rendering as dots, shallower levels as bounding boxes.
- Tags (generalizing the current "bad" flag) drive per-node color via a priority-ordered catalog and inherit to descendants.
- Per-node notes inherit ancestor notes as read-only sections, merged into one window.
- All of this dovetails with the existing tree panel (which stays) and the existing measurement plot panel via shared selection state.

Alongside this feature, `Gui.jl` (3,973 LOC, single module) is split into a small set of files for maintainability — driven by the new feature surface, not as a separate refactor.

## Decisions locked in

- **Renderer**: ImDrawList (CImGui primitives). No Makie second canvas.
- **Storage**: focused files at the source root, lean and human-editable — `devices_info.txt` extended with `x_um/y_um/w_um/h_um`; new `layout.txt`, `tags.txt` (catalog + assignments in one file), `notes.txt`.
- **Tree panel**: stays. Spatial browser is an additional dockable panel sharing `ui_state[:selected_device_paths]`.
- **Level semantics**: depth from root. Level 0 = top containers; level N = "focusable" depth, deeper hidden, shallower rendered as bboxes.

## High-level architecture

### Data layer (new module: `src/Annotations/`)

A new internal subpackage (mirrors the structure of `DataLoader`, `DataAnalysis`, `DataPlotter`) housing all on-disk metadata that isn't filenames or measurement contents.

```
src/Annotations/
├── Project.toml
├── src/
│   ├── Annotations.jl        — module entry, re-exports
│   ├── Coords.jl             — read/write x/y/w/h from devices_info.txt; bbox computation
│   ├── Layout.jl             — layout.txt: per-path world-XY for unpositioned containers
│   ├── Tags.jl               — tags.txt: catalog + per-path assignments in one file
│   ├── Notes.jl              — notes.txt with [path]+fenced sections; ancestor merging
│   └── BadRegistryShim.jl    — BadRegistry expressed as a built-in tag for migration
```

**Key APIs (sketch):**

- `Coords.node_position(node) :: Maybe{(x, y)}` — explicit if in `devices_info.txt`, `nothing` otherwise.
- `Coords.bounding_box(hierarchy, path) :: Rect` — computed from descendant positions; cached per scan.
- `Layout.load(root) :: Dict{String, (x, y)}` / `Layout.save(root, dict)`.
- `Layout.reset!(state, paths)` — re-seeds positions for the given paths using a default grid, called by the UI's "reset layout" button.
- `Tags.load(root) :: TagState` (parses both the catalog and assignments from a single `tags.txt`).
- `Tags.catalog(state) :: Vector{TagDef(name, color, priority)}`.
- `Tags.assignments(state) :: Dict{String, Set{String}}` (path → tag names; explicit only).
- `Tags.effective(path, assignments) :: Set{String}` (own tags ∪ inherited from ancestors).
- `Tags.dominant_color(path, catalog, assignments) :: RGBA` (highest-priority tag wins).
- `Notes.read_section(root, path) :: String` — own section only.
- `Notes.merged_view(root, path) :: Vector{Section{path, text, editable}}` — ancestors read-only + self editable.
- `Notes.write_section!(root, path, text)`.

**Storage formats (human-editable):**

All files live at the source root, kept lean so they read well in any text editor.

- `devices_info.txt` — extended with optional `x_um, y_um, w_um, h_um` columns (path-prefix matching as today). Positioned entries set `x_um, y_um`; container entries can additionally set `w_um, h_um` to override computed bbox.
- `layout.txt` — `<path>\t<x_um>\t<y_um>` lines for the user-arranged unpositioned containers.
- `tags.txt` — single file with two top sections, `[catalog]` and `[assignments]`. Catalog rows: `<name>\t<color_hex>\t<priority>`. Assignment rows: `<path>\t<tag_name>`. Ships with a default `bad` catalog entry; the legacy `bad_measurements` file is auto-migrated on first save.
- `notes.txt` — fenced sections, robust against `[brackets]` appearing inside note bodies:
  ```
  [ChipB]
  ```
  Oxygen flow: 5%
  ```
  [ChipB/SiteVI]
  ```
  dust particle visible in the middle
  ```
  ```
  Each section is `[<path>]` followed by a triple-backtick-fenced body. Parser keys on the fences, not on the brackets, so users can write `[anything]` freely inside.

### UI layer (`Gui.jl` split)

`Gui.jl` (3,973 LOC) decomposes into a `Gui/` directory. Target: ~7 files, no further. Each file exposes a small surface; state stays in the existing `ui_state::Dict` (typed access helpers added but no big-bang struct refactor).

```
src/Gui.jl                      — module entry; re-exports create_window_and_run_loop
src/Gui/
├── State.jl                    — ui_state init helpers + typed accessors
├── Layout.jl                   — _setup_docking_layout!, top menu, frame loop
├── TreePanel.jl                — current tree (move from Gui.jl as-is)
├── PlotPanel.jl                — plot area + plot job queue + combined plots
├── SpatialBrowser.jl           — NEW. Pan/zoom canvas, hit-testing, level controls
├── AnnotationsUI.jl            — context menus, tag catalog editor, notes window
└── InfoModal.jl                — existing device modal + figure scripts modal
```

The split is driven by the natural seams the new feature exposes (tree, plot, spatial, annotations, info), not by lines-of-code targets.

### Spatial browser internals (`SpatialBrowser.jl`)

Render every frame inside a CImGui child window:

1. **Camera**: pan offset (world μm) + zoom (px per μm). Mouse wheel zooms at cursor; right-drag (or middle-drag) pans.
2. **Level state**: `ui_state[:spatial_level] :: Int`. `+/-` keys and on-canvas buttons increment/decrement. Display "level: N" in a corner overlay.
3. **Visibility & rendering** — for each node `n` at path depth `d`:
   - `d < level` and `has_children`: render its bounding box (transparent fill, thin outline). Click-transparent unless `d == level`.
   - `d == level`: render as the **focusable** primitive — bbox if it has children, dot otherwise. Hit-target.
   - `d > level`: skip (hidden).
4. **Labels**: every focusable node renders its name (last path segment) at the top-left of its bbox/dot in a tiny font. Drawn after the shape so it's always on top; LOD-culled when the node footprint is below a few pixels.
5. **Color**: `Tags.dominant_color(path, ...)`; default neutral if no tags. Computed once per frame per visible node, cached when assignments don't change.
6. **Selection**: shared with tree via `ui_state[:selected_device_paths]`. Click toggles, ctrl-click adds, shift-click range-selects (in spatial it's "all in current marquee path", deferring the exact semantics to implementation), drag-on-empty starts a marquee.
7. **Drag**: only nodes whose `Coords.node_position` is `nothing` (i.e., chips) and only when `d == level`. Dragging updates `Layout` in memory; release writes to disk.
8. **Right-click**: opens context menu with **Notes…**, **Tags ▸**, **Info…**, **Mark/unmark bad** (compatibility shortcut).

**Performance**: viewport-cull via the camera's world rect against each node's bbox. Quadtree only if profiling shows we need it; for thousands of nodes ImDrawList is fine.

### Selection sync

`ui_state[:selected_device_paths]` is the single source of truth. Both tree (`Gui.jl:2788 _update_multi_selection!`) and spatial browser write to it; both read from it for highlight rendering. The existing `_apply_visible_selection!` (`Gui.jl:654`) continues to drive the measurement panel.

### Notes window UI

A standalone, modeless ImGui window — not anchored to any node. Opened from the right-click menu of any node; multiple windows can be open at once (one per opened node), independently movable, dockable, and closable.

- Window title is the node's full path so multiple windows are distinguishable.
- Body lists each ancestor section in path order, then the node itself.
- Each section: bolded title (the segment name, e.g. **ChipB**), a horizontal rule below the body, body text.
- Ancestor bodies are read-only and rendered slightly desaturated.
- Own body is editable (ImGui multiline input). Save-on-blur or explicit Save button — implementation choice deferred.
- Window state lives in `ui_state[:notes_windows] :: Vector{NotesWindowState}` so reopening the same node focuses the existing window instead of creating a duplicate.

## Phased implementation roadmap

Each phase is a separate spawn: small, reviewable, behavior-checked before moving on.

1. **P1 — Gui.jl split, no behavior change.** Move code to `Gui/` files, keep public API identical, run smoke test (`julia --project start.jl test/fixtures/...`). Critical: no logic edits.
2. **P2 — `Annotations` subpackage, data model only.** Implement `Coords`, `Layout`, `Tags`, `Notes`, `BadRegistryShim`. Tests with fixture files. Existing `BadRegistry` keeps working via shim.
3. **P3 — Tagging supersedes BadRegistry.** Wire `Tags` (catalog + assignments) into the existing tree panel and measurement filtering: the bad-toggle context menu, the bad styling in the tree, and `_device_is_visible` / `_measurement_is_visible` all read from `Tags` instead of `BadRegistry`. On startup, `Tags.load` migrates any legacy `bad_measurements` file into `tags.txt` (catalog gets a `bad` entry if missing; each line becomes an assignment). `BadRegistry` is removed once parity is verified — the new system is a strict superset, not a parallel one. No spatial UI yet; this is purely a data-layer swap with the existing UI.
4. **P4 — Spatial browser shell.** New panel, pan/zoom, render only the unpositioned-roots layer (chips with default grid initialization + "reset layout" button), drag-to-rearrange, persistence. No levels yet. Tag colors already work because P3 wired them.
5. **P5 — Levels + bbox rendering.** Compute hierarchical positions, level state, +/- controls, render bboxes vs dots per the rules. Hit-testing across levels. Marquee selection. Labels at top-left.
6. **P6 — Notes UI.** Right-click → notes window. Ancestor merging. Editable own section. Multi-window support.
7. **P7 — Polish.** Filters (toggle-show by tag), keyboard shortcuts, save-state in project preferences, info modal cleanup.

Each phase ends with a manual verification run + git commit.

## Critical files (touched or created)

**Modified:**
- [src/Gui.jl](src/Gui.jl) — split into `Gui/` files (P1).
- [src/MeasurementBrowser.jl](src/MeasurementBrowser.jl) — `include` paths for the new layout.
- [src/DeviceParser.jl](src/DeviceParser.jl:442) — `_load_scan_metadata` extended for new optional columns (`x_um, y_um, w_um, h_um`); no behavior change for files that don't use them.
- [src/projects/RuO2/Cache.jl:12](src/projects/RuO2/Cache.jl:12) — already reserves `:x, :y`; ensure cache writes flow through.
- [Project.toml](Project.toml) — add `Annotations` as a path-deps subpackage.

**New:**
- `src/Gui/{State,Layout,TreePanel,PlotPanel,SpatialBrowser,AnnotationsUI,InfoModal}.jl`
- `src/Annotations/` (full subpackage as outlined).
- Test fixtures under `test/fixtures/spatial/` with sample `devices_info.txt`, `tags.txt`, `notes.txt`, `layout.txt`.

## Verification

End-to-end smoke after each phase:

```bash
julia --project -e 'using Pkg; Pkg.instantiate()'
julia --project=src/Annotations -e 'using Pkg; Pkg.instantiate()'   # P2 onwards
julia --project -e 'using Pkg; Pkg.test()'
julia --project start.jl test/fixtures/spatial   # manual UI check
```

Per-phase manual checks:
- **P1**: launch UI, confirm tree/plot/menus/modal all behave identically.
- **P2**: unit tests against fixture files for each storage format.
- **P3**: launch a project that has an existing `bad_measurements` file → confirm it's migrated to `tags.txt`, `bad_measurements` is gone (or kept as a one-time backup), all previously-bad items still render bad in the tree, mark/unmark still works, filtering still works. Also: assign a non-bad tag manually in `tags.txt` and confirm it's visible in the tree styling.
- **P4**: open spatial panel, see chips initialized as a grid, drag → reopen → positions persisted; "reset layout" button re-grids them; tag colors already flow through.
- **P5**: +/- changes level; bboxes/dots render per rules with labels at top-left; marquee selects; selection mirrors tree.
- **P6**: walk the example workflow from the user prompt verbatim; ancestor sections appear read-only, own section edits persist to `notes.txt`; opening notes for two different nodes yields two independent windows.

## Deferred decisions / TBDs

- **Save-on-blur vs explicit Save button** for notes editing.
- **Marquee selection semantics across levels** — does a marquee at level N select only level-N items, or also descendants? Default: level-N only; implementation will revisit.
- **Tag catalog editor placement** — separate modal or inline in `AnnotationsUI`. Will decide in P5.
- **Quadtree** for hit-testing — only if P4 profiling shows we need it.
- **"Info window (TBD)"** from the user spec — kept TBD; existing device modal can be the placeholder until we know what's wanted.
