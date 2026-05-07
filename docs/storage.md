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

## File ownership matrix

| File | Concern | Read by | Written by |
|---|---|---|---|
| `devices_info.txt` | Per-path parameters | `_load_scan_metadata` ([DeviceParser.jl:497](../src/DeviceParser.jl)) | Hand-edited |
| `bad_measurements` | Bad flag | `load_bad_registry` | `save_bad_registry` |

## Conventions for metadata file

- Live at source root. No subdirectories unless there's volume reason.
- Prefer line-oriented text. CSV if there are columns; key/value with simple delimiters otherwise.
- Use `device_path_key(location)` (slash-joined) for any path-keyed entry.
- Provide a load + save function pair in a focused module/file.
- Update this doc and any relevant cross-references in [ARCHITECTURE.md](ARCHITECTURE.md).
- Tolerate missing/malformed files quietly on read; fail loudly on write.

## Cache layer (separate from source-root metadata)

The HDF5 cache lives elsewhere and is regenerated from CSVs. The RuO2 cache schema reserves `:x, :y` slots on devices ([src/projects/RuO2/Cache.jl:12](../src/projects/RuO2/Cache.jl)).
