# Data Model

> Types and identity for measurements and the hierarchy tree. Defined in [src/DeviceParser.jl](../src/DeviceParser.jl).

## Core types

```julia
struct DeviceInfo
    location::Vector{String}           # ["RuO2test", "A9", "VI", "D1"]
    parameters::Dict{Symbol,Any}       # merged from device_info.txt (area_um2, t_HZO_nm, …)
end

struct MeasurementInfo
    id::String                         # filesystem-stable ID
    filename::String
    filepath::String
    clean_title::String
    measurement_kind::Symbol           # :pund, :iv, :breakdown, :tlm4p, :wakeup, :cvsweep
    timestamp::Union{DateTime,Nothing}
    device_info::DeviceInfo
    parameters::Dict{Symbol,Any}       # parsed from filename (voltage_V, frequency_Hz, …)
    wakeup_pulse_count::Union{Int,Nothing}
end

struct HierarchyNode
    name::String                       # one path segment
    kind::Symbol                       # :root | :level (intermediate) | :leaf (last segment)
    children::Vector{HierarchyNode}
    measurements::Vector{MeasurementInfo}   # populated only on :leaf nodes
end

struct MeasurementHierarchy
    root::HierarchyNode
    all_measurements::Vector{MeasurementInfo}
    root_path::String
    index::Dict{Tuple{Vararg{String}}, HierarchyNode}   # path tuple → node
    has_device_metadata::Bool
    project::AbstractProject
    skipped_count::Int
end
```

## Variable hierarchy depth

Depth is **not** uniform. One branch may be `RuO2test/A9/VI/D1` (4 segments) while another is `RuO2test_A10/VI/25um/D1` (also 4) and another could be 3 or 5. Don't assume "level 3 = device".

- All measurements attach at `:leaf` nodes (the last segment).
- Intermediate nodes are `:level`.
- A leaf is **the last segment of its path**, not "depth N from root."

## Identity / path keys

Two equivalent representations:

- **Tuple** (used for `index` lookup): `("RuO2test", "A9", "VI", "D1")`.
- **String** (used in on-disk files): `"RuO2test/A9/VI/D1"` — built by `device_path_key(location)` ([DeviceParser.jl:65](../src/DeviceParser.jl)).
- **Round-trip**: `device_path_tuple("RuO2test/A9/VI/D1")` parses + validates ([DeviceParser.jl:68](../src/DeviceParser.jl)).

The slash-joined string is the canonical identifier in any on-disk metadata file. New metadata systems should reuse it.

## Path-prefix matching

`device_info.txt` rows are keyed by path; metadata applies to **all descendants whose location is prefixed by that key**, with longer (more specific) keys overriding shorter ones. See [scanning.md](scanning.md) for the merge logic.

Implication: hierarchical metadata is stored once at the highest applicable level. Inheritance is "lookup-time" — children don't physically carry parent metadata in their structs.

## Virtual measurement expansion

Some single CSVs become multiple `MeasurementInfo` entries during scan:

- **Breakdown files** → `expand_multi_device` splits multi-device sweeps into per-device entries.
- **PUND fatigue files** → `expand_pund_fatigue` creates one virtual `:pund` per cycle, storing `parameters[:fatigue_cycle]`.

Virtual measurements share a filepath but get distinct `id`s.

## Detection ordering matters

`detect_measurement_kind` ([DeviceParser.jl](../src/DeviceParser.jl)) checks filename patterns in order. Specific patterns must come **before** general ones — `pund_fatigue` before `pund`, etc. Adding a new kind: insert it at the correct precedence.

## Where to look

| Concern | File |
|---|---|
| Type definitions | [src/DeviceParser.jl](../src/DeviceParser.jl) lines 58–106 |
| Hierarchy construction | [src/DeviceParser.jl](../src/DeviceParser.jl) lines 283–344 |
| Path key helpers | [src/DeviceParser.jl](../src/DeviceParser.jl) lines 65–74 |
| `scan_source` | [src/DeviceParser.jl](../src/DeviceParser.jl) lines 630–676 |
| Metadata merge (path-prefix) | [src/DeviceParser.jl](../src/DeviceParser.jl) lines 477–493 |
| Fixtures showing tree shape | [test/test_scan_directory_progress.jl](../test/test_scan_directory_progress.jl) |
