# On-Disk Metadata

> Files at the **source root** (the measurement folder the user opens, not the repo). All hand-editable plain text. Each file owns one orthogonal concern; resist the urge to merge them.

## devices_info.txt

CSV. First column is a path (segment-joined with `/`); remaining columns are arbitrary key/value parameters merged into matching nodes.

```
device_path,area_um2,t_HZO_nm,notes,active
D1,50.0,7.0,all exact D1 units,true
A9,,legacy_site,
RuO2test/A9/VI/D1,10.0,7.0,exact match wins,true
RuO2test_A10/VI/25um/D1,25.0,7.0,nested scope,
```

**Path-prefix matching** ([DeviceParser.jl:477](../src/DeviceParser.jl)): a row applies to every descendant whose location starts with that key. Longer (more specific) keys override shorter ones. So `D1` sets defaults for every `D1` leaf; `RuO2test/A9/VI/D1` overrides for that exact leaf.

Field types are inferred (bool / int / float / date / string).

## bad_measurements

Flat list, no extension, no header:

```
# comments allowed
device RuO2test/A9/VI/D1
measurement <measurement_id>
```

Loaded by `BadRegistry` ([src/BadRegistry.jl](../src/BadRegistry.jl)). Loaded at project init in [src/Gui/BadAndStyling.jl](../src/Gui/BadAndStyling.jl) via `_load_bad_registry_for_root!`; persisted via `save_bad_registry`.

`Annotations.Tags.load` reads `bad_measurements` whenever it is present, regardless of whether `tags.txt` also exists. Entries are merged as `bad` assignments: `device <path>` lines go into `assignments` (device-path-keyed map); `measurement <id>` lines go into `measurement_assignments` (measurement-ID-keyed map). If the catalog has no `bad` entry, one is added with color `(0xff, 0x30, 0x30)` and priority `100`. `load` does not write to disk. After `load → save`, the merged state is encoded in `tags.txt`, so subsequent `load → save` cycles produce byte-identical output as long as `bad_measurements` has not grown.

## layout.txt

Tab-separated, one record per line, no header:

```
<path>	<x_um>	<y_um>
```

Lines starting with `#` and blank lines are ignored. Saved sorted by path. Loaded and saved by [src/Annotations/src/Layout.jl](../src/Annotations/src/Layout.jl) (`Annotations.Layout.load` / `Annotations.Layout.save`). Saving an empty map removes the file. Malformed rows raise `LayoutParseError`.

## tags.txt

Two bracket-headed sections, tab-separated rows:

```
[catalog]
bad	ff3030	100
todo	30c0ff	50

[assignments]
device	RuO2test/A9/VI/D1	bad
measurement	abc123hash	bad
device	RuO2test/A10/VI	todo
```

Catalog rows: `<name>\t<color_hex_rrggbb>\t<priority>`. Assignment rows: `<kind>\t<key>\t<tag_name>`, where `<kind>` is either `device` (device-path key) or `measurement` (measurement-ID key). Unknown kind tokens raise `TagsParseError`. Fields are tab-separated on write; whitespace-tolerant on read. Lines starting with `#` and blank lines are ignored. Loaded and saved by [src/Annotations/src/Tags.jl](../src/Annotations/src/Tags.jl) (`Annotations.Tags.load` / `Annotations.Tags.save`). Saving empty state removes the file. Malformed rows raise `TagsParseError`.

## notes.txt

Fenced sections, robust against `[brackets]` inside note bodies:

````
[ChipB]
```
Oxygen flow: 5%
```
[ChipB/SiteVI]
```
dust particle visible
[observed in second pass]
```
````

Each section is `[<path>]` on its own line, followed by an opening triple-backtick fence, then body lines, then a closing fence. The parser keys on the fences, not the brackets, so `[anything]` inside a body is preserved verbatim. The trailing newline before the closing fence is trimmed; other whitespace is preserved. Loaded and saved by [src/Annotations/src/Notes.jl](../src/Annotations/src/Notes.jl) (`Annotations.Notes.read_section`, `Annotations.Notes.merged_view`, `Annotations.Notes.write_section!`). Malformed sections raise `NotesParseError`.

## File ownership matrix

| File | Concern | Read by | Written by |
|---|---|---|---|
| `devices_info.txt` | Per-path parameters | `_load_scan_metadata` ([DeviceParser.jl:497](../src/DeviceParser.jl)) | Hand-edited |
| `bad_measurements` | Bad flag | `load_bad_registry`; `Annotations.Tags.load` (merged as `bad` assignments) | `save_bad_registry` |
| `layout.txt` | User-arranged XY positions | `Annotations.Layout.load` | `Annotations.Layout.save` |
| `tags.txt` | Tag catalog + assignments | `Annotations.Tags.load` | `Annotations.Tags.save` |
| `notes.txt` | Per-path note bodies | `Annotations.Notes.read_section`, `Annotations.Notes.merged_view` | `Annotations.Notes.write_section!` |

## Conventions for metadata file

- Live at source root. No subdirectories unless there's volume reason.
- Prefer line-oriented text. CSV if there are columns; key/value with simple delimiters otherwise.
- Use `device_path_key(location)` (slash-joined) for any path-keyed entry.
- Provide a load + save function pair in a focused module/file.
- Update this doc and any relevant cross-references in [ARCHITECTURE.md](ARCHITECTURE.md).
- Tolerate missing/malformed files quietly on read; fail loudly on write.

## Cache layer (separate from source-root metadata)

The HDF5 cache lives elsewhere and is regenerated from CSVs. The RuO2 cache schema reserves `:x, :y` slots on devices ([src/projects/RuO2/Cache.jl:12](../src/projects/RuO2/Cache.jl)).
