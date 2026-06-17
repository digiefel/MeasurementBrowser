# Annotation Storage

> Plain-text files describing the formats below. Each file owns one orthogonal concern; resist the urge to merge them.

These files annotate collections and items. Figure annotations are stored separately.

> **Location (current, stopgap):** these files live **next to the cache**, keyed by `source_id`, so a
> non-filesystem `AbstractDataSource` has somewhere to persist them. The regression is that they are no
> longer hand-editable at a data folder; the intended fix makes annotation storage a source capability
> (file-backed sources expose hand-editable source-root files). See the TODO in [cache.md](cache.md).
> The formats below are unchanged regardless of where the files live.

## device_info.txt

CSV. First column is a path (segment-joined with `/`); remaining columns are arbitrary key/value parameters merged into matching nodes.

```
collection_path,area_um2,t_HZO_nm,notes,active
D1,50.0,7.0,all exact D1 units,true
A9,,legacy_site,
RuO2test/A9/VI/D1,10.0,7.0,exact match wins,true
RuO2test_A10/VI/25um/D1,25.0,7.0,nested scope,
```

**Path-fragment matching**: a row applies when its slash-separated path exactly matches a
contiguous fragment of the parsed device location. Longer matches are merged later, so more specific
rows override shorter ones. `D1` applies to every device path containing a `D1` segment, while
`RuO2test/A9/VI/D1` applies to that exact sequence.

Field types are inferred (bool / int / float / date / string).

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
RuO2test/A9/VI/D1                                          bad
/Users/davide/data/RuO2test/A9/VI/D1/3V FE PUND.csv        bad
RuO2test/A10/VI                                            todo
```

Catalog rows: `<name>\t<color_hex_rrggbb>\t<priority>`. Assignment rows: `<key>\t<tag_name>`. Keys are either device-path strings (slash-joined segments, e.g. `RuO2test/A9/VI/D1`) or measurement-ID strings (absolute filesystem paths, optionally with `#cycle=N` / `#split=X` suffixes). The two namespaces never overlap, so no kind prefix is written. Fields are tab-separated on write; whitespace-tolerant on read. Lines starting with `#` and blank lines are ignored. Loaded and saved by [src/Annotations/src/Tags.jl](../src/Annotations/src/Tags.jl) (`Annotations.Tags.load` / `Annotations.Tags.save`). Saving empty state removes the file. Malformed rows raise `TagsParseError`.

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
| `device_info.txt` | Per-collection parameters | `ItemIndex.load_scan_metadata` | Hand-edited |
| `layout.txt` | User-arranged XY positions | `Annotations.Layout.load` | `Annotations.Layout.save` |
| `tags.txt` | Tag catalog + assignments | `Annotations.Tags.load` | `Annotations.Tags.save` |
| `notes.txt` | Per-path note bodies | `Annotations.Notes.read_section`, `Annotations.Notes.merged_view` | `Annotations.Notes.write_section!` |

## Conventions for metadata file

- Live at source root. No subdirectories unless there's volume reason.
- Prefer line-oriented text. CSV if there are columns; key/value with simple delimiters otherwise.
- Use `collection_path_key(collection)` (slash-joined) for any path-keyed entry.
- Provide a load + save function pair in a focused module/file.
- Update this doc and any relevant cross-references in [ARCHITECTURE.md](ARCHITECTURE.md).
- Tolerate missing/malformed files quietly on read; fail loudly on write.

## Cache layer (separate from source-root metadata)

The HDF5 cache is generated data, not hand-edited source-root metadata. It lives outside the
measurement folder and is documented in [cache.md](cache.md).
