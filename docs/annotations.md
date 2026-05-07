# Annotations Subpackage

> Per-path metadata that isn't filenames or measurement contents — coordinates, layout, tags, notes. Lives at [src/Annotations/](../src/Annotations/) as a path-deps subpackage.

## Purpose

`Annotations` owns the on-disk side of the annotation system. Each module pairs a small, hand-editable text file at the source root with a focused load/save API. Path keys are slash-joined strings (`"RuO2test/A9/VI/D1"`), shared with the rest of the codebase via `device_path_key`.

The subpackage is consumed via `using Annotations` and namespaced access (`Annotations.Tags.load(root)`, `Annotations.Notes.read_section(root, path)`, etc.). Submodule names are re-exported from the top-level `Annotations` module so the namespace works from outside.

## Module map

| Module | File | Concern |
|---|---|---|
| `Coords` | [src/Annotations/src/Coords.jl](../src/Annotations/src/Coords.jl) | Pulls `x_um, y_um, w_um, h_um` out of a parsed `devices_info` table; computes axis-aligned bounding boxes. |
| `Layout` | [src/Annotations/src/Layout.jl](../src/Annotations/src/Layout.jl) | Reads/writes `layout.txt` (user-dragged container positions); generates default grid layouts. |
| `Tags` | [src/Annotations/src/Tags.jl](../src/Annotations/src/Tags.jl) | Reads/writes `tags.txt` (catalog + per-path assignments); resolves effective tags and dominant color. |
| `Notes` | [src/Annotations/src/Notes.jl](../src/Annotations/src/Notes.jl) | Reads/writes `notes.txt` fenced sections; merges ancestor sections for display. |

## Public API

### `Coords`

- `read_positions(devices_info) -> Dict{String, (x, y)}` — rows that supply both `x_um` and `y_um`.
- `read_overrides(devices_info) -> Dict{String, (w, h)}` — rows that supply both `w_um` and `h_um`.
- `bounding_box(positions, descendant_paths; override=nothing) -> Union{Rect, Nothing}` — minima/maxima over the positions of `descendant_paths`. Missing paths are skipped. Returns `nothing` when no descendants resolve. `override=(w, h)` keeps the lower-left corner from the descendants but uses the given size for the width/height.
- `Rect(x, y, w, h)` — value type; `(x, y)` is the lower-left corner, in micrometres.

### `Layout`

- `load(root) -> PositionMap` — `Dict{String, (x, y)}`. Missing file returns an empty map. Malformed rows raise `LayoutParseError`.
- `save(root, positions)` — sorted by path. An empty map removes the file.
- `reset!(positions, paths; cols=nothing, spacing_um=200.0, origin=(0.0, 0.0))` — overwrites entries for `paths` with a row-major grid. `cols` defaults to `ceil(sqrt(n))`.

### `Tags`

- `TagDef(name, color::NTuple{3,UInt8}, priority::Int)` — single catalog entry.
- `TagState(catalog::Vector{TagDef}, assignments::Dict{String, Set{String}}, measurement_assignments::Dict{String, Set{String}})` — full state. `assignments` is device-path-keyed; `measurement_assignments` is measurement-ID-keyed. `TagState()` is the empty state.
- `load(root) -> TagState`. Reads `tags.txt` if present, then reads `bad_measurements` if present and merges its entries as `bad` assignments into both maps. See [storage.md](storage.md#bad_measurements) and [storage.md](storage.md#tags.txt). `load` never writes to disk.
- `save(root, state)` — writes only `tags.txt`. Empty state removes the file.
- `effective(state, path, ancestor_paths) -> Set{String}` — union of `path`'s device-path assignments with assignments on every entry of `ancestor_paths`. Does not include measurement-ID tags.
- `assigned_to_measurement(state, measurement_id) -> Set{String}` — explicit tags on `measurement_id`; empty set if none. Measurement-ID tags are not ancestor-walked; callers compose `effective` and `assigned_to_measurement` to get the full applicable set for a measurement.
- `dominant_color(state, effective_tags) -> Union{Nothing, NTuple{3,UInt8}}` — highest-priority hit's color among catalog entries whose name is in `effective_tags`. Returns `nothing` for an empty input or no catalog matches.

### `Notes`

- `read_section(root, path) -> String` — body for `path`, or `""` when absent.
- `merged_view(root, path, ancestor_paths) -> Vector{NamedTuple{(:path, :body, :editable)}}` — ancestor sections in input order (each `editable=false`), then `path` itself (`editable=true`, included even when empty).
- `write_section!(root, path, body)` — replaces the section in place if present, otherwise appends. Other sections preserved in original order.

## File formats

See [storage.md](storage.md) for the on-disk shape of `layout.txt`, `tags.txt`, `notes.txt`, plus the `x_um, y_um, w_um, h_um` columns on `devices_info.txt`.

## Inheritance semantics

Annotations attach to specific paths. Two modules expose ancestor inheritance at lookup time; the caller supplies the ancestor list (closest-first or any order — the modules don't care).

- **`Tags.effective`** — set union of own tags and every ancestor's tags. No precedence among ancestors; tags are membership, not values. Use `dominant_color` to pick a single color out of the resulting set via the catalog's `priority` field.
- **`Notes.merged_view`** — each ancestor section is included read-only in the order supplied; the focal node's section is appended last and marked editable. Ancestors that have no section are skipped silently.

Neither module rewrites stored state when it inherits — assignments and note bodies stay anchored to the path they were authored against.

## Where to look

| Concern | Location |
|---|---|
| Module entry point | [src/Annotations/src/Annotations.jl](../src/Annotations/src/Annotations.jl) |
| Subpackage `Project.toml` | [src/Annotations/Project.toml](../src/Annotations/Project.toml) |
| Tests | [test/test_annotations.jl](../test/test_annotations.jl) |
| Test fixtures | [test/fixtures/annotations/](../test/fixtures/annotations) |
| File formats | [storage.md](storage.md) |
| Path key conventions | [data-model.md](data-model.md#identity--path-keys) |
