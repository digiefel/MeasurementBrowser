# Annotations

Annotations are user-editable information attached to collections and items: spatial
coordinates, saved positions, tags, and notes. They are separate from item data, generated
cache data, browser state, and annotations drawn on figures.

## Purpose

The implementation is the path-dependency package
[`DataBrowserAnnotations`](../lib/DataBrowserAnnotations/). Each namespace reads or writes one focused source-root
format. Path keys are
slash-joined strings such as `"RuO2test/A9/VI/D1"`, shared with the rest of the package through
`collection_path_key`.

## Current API

Code imports `DataBrowserAnnotations` and uses four namespaces:

| API | Purpose |
|---|---|
| `Annotations.Coords` | Reads `x_um`, `y_um`, `w_um`, and `h_um` from parsed `metadata.txt` metadata and computes bounding boxes. |
| `Annotations.Layout` | Reads and writes user-arranged positions in `layout.txt`. |
| `Annotations.Tags` | Reads and writes tag definitions and collection/item assignments in `tags.txt`. |
| `Annotations.Notes` | Reads and writes inherited per-path notes in `notes.txt`. |

## Current GUI

Only tags are connected to the browser today. Collection and item context menus show
`Mark Bad` and `Unmark Bad`. The menu bar provides `Show Bad`, and effective tag colors style tree
and item text. A malformed `tags.txt` is shown as a tag error.

Coordinates, saved spatial layout, and notes currently have no GUI. They exist for the planned
spatial browser and should not be described as existing user-facing features.

## Detailed API

### `Coords`

- `read_positions(metadata) -> Dict{String, (x, y)}` — rows that supply both `x_um` and `y_um`.
- `read_overrides(metadata) -> Dict{String, (w, h)}` — rows that supply both `w_um` and `h_um`.
- `bounding_box(positions, descendant_paths; override=nothing) -> Union{Rect, Nothing}` — minima/maxima over the positions of `descendant_paths`. Missing paths are skipped. Returns `nothing` when no descendants resolve. `override=(w, h)` keeps the lower-left corner from the descendants but uses the given size for the width/height.
- `Rect(x, y, w, h)` — value type; `(x, y)` is the lower-left corner, in micrometres.

### `Layout`

- `load(root) -> PositionMap` — `Dict{String, (x, y)}`. Missing file returns an empty map. Malformed rows raise `LayoutParseError`.
- `save(root, positions)` — sorted by path. An empty map removes the file.
- `reset!(positions, paths; cols=nothing, spacing_um=200.0, origin=(0.0, 0.0))` — overwrites entries for `paths` with a row-major grid. `cols` defaults to `ceil(sqrt(n))`.

### `Tags`

- `TagDef(name, color::NTuple{3,UInt8}, priority::Int)` — single catalog entry.
- `TagState(catalog::Vector{TagDef}, assignments::Dict{String, Set{String}})` — full state. `assignments` holds all explicitly attached tags, keyed by any string: collection paths (slash-joined segments, e.g. `"RuO2test/A9/VI/D1"`) and item ids share the same map and never collide. `TagState()` is the empty state.
- `load(root) -> TagState`. Reads `tags.txt` when present. Missing file returns an empty state.
- `save(root, state)` — writes `tags.txt`. Empty state removes the file.
- `effective(state, key, ancestor_keys) -> Set{String}` — union of `key`'s own assignments with assignments on every entry of `ancestor_keys`. To get the full applicable tag set for an item, call `effective(state, id, [collection_path; collection_ancestors...])`: the item id and the collection-path ancestors are looked up uniformly in the same map.
- `dominant_color(state, effective_tags) -> Union{Nothing, NTuple{3,UInt8}}` — highest-priority hit's color among catalog entries whose name is in `effective_tags`. Returns `nothing` for an empty input or no catalog matches.

### `Notes`

- `read_section(root, path) -> String` — body for `path`, or `""` when absent.
- `merged_view(root, path, ancestor_paths) -> Vector{NamedTuple{(:path, :body, :editable)}}` — ancestor sections in input order (each `editable=false`), then `path` itself (`editable=true`, included even when empty).
- `write_section!(root, path, body)` — replaces the section in place if present, otherwise appends. Other sections preserved in original order.

## File formats

See [storage.md](storage.md) for the on-disk shape of `layout.txt`, `tags.txt`, `notes.txt`, plus the `x_um, y_um, w_um, h_um` columns on `metadata.txt`.

## Inheritance semantics

Tags and notes attach to specific paths. Both expose ancestor inheritance at lookup time; the caller
supplies the ancestor list.

- **`Tags.effective`** — set union of own tags and every ancestor's tags. Keys are arbitrary strings; collection-path keys and item keys are looked up in the same `assignments` map. No precedence among keys; tags are membership, not values. Use `dominant_color` to pick a single color out of the resulting set via the catalog's `priority` field.
- **`Notes.merged_view`** — each ancestor section is included read-only in the order supplied; the focal node's section is appended last and marked editable. Ancestors that have no section are skipped silently.

Neither module rewrites stored state when it inherits — assignments and note bodies stay anchored to the path they were authored against.

See [storage.md](storage.md) for file formats and
[data-model.md](data-model.md#identity) for identity conventions.
