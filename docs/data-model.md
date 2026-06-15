# Data Model

`MeasurementIndex` defines physical source files, logical measurements, and the device hierarchy.
Projects interpret indexed source files; package code stores and presents the resulting measurements.

## Core types

```julia
struct DeviceInfo
    location::Vector{String}           # ["RuO2test", "A9", "VI", "D1"]
    parameters::Dict{Symbol,Any}       # merged from device_info.txt (area_um2, t_HZO_nm, …)
end

struct MeasurementInfo
    unique_id::String                  # stable identity for one logical measurement
    filename::String
    filepath::String
    clean_title::String
    measurement_kind::Symbol
    timestamp::Union{DateTime,Nothing}
    device_info::DeviceInfo
    parameters::Dict{Symbol,Any}       # acquisition settings known while interpreting the file
    stats::Dict{Symbol,Any}            # values computed after the required context is available
end

struct SourceFile
    unique_id::String
    filepath::String
    filename::String
    timestamp::Union{DateTime,Nothing}
    header_summary::Dict{String,String}
    fingerprint::FileFingerprint
    measurements::Vector{MeasurementInfo}
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
    project::Project
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
- **String** (used in on-disk files): `"RuO2test/A9/VI/D1"` — built by `device_path_key(location)`.
- **Round-trip**: `device_path_tuple("RuO2test/A9/VI/D1")` parses and validates the string form.

The slash-joined string is the canonical identifier in any on-disk metadata file. New metadata systems should reuse it.

## Path-prefix matching

`device_info.txt` rows are keyed by path; metadata applies to **all descendants whose location is prefixed by that key**, with longer (more specific) keys overriding shorter ones.

The on-disk format is documented in [storage.md](storage.md). The scan applies matching parameters
before measurements enter the hierarchy.

Implication: hierarchical metadata is stored once at the highest applicable level. Inheritance is "lookup-time" — children don't physically carry parent metadata in their structs.

## Virtual measurement expansion

A single source file may yield several `MeasurementInfo` entries when the project's `measurements`
callback returns a vector with more than one entry. They share a filepath but get distinct
`unique_id` values, which the engine mints from filepath + kind + the `parameters` that distinguish
siblings.

Expansion is purely a project concern. For example, a project may split a multi-device sweep into one
entry per device, or expand a fatigue file into one entry per cycle (storing the cycle number in
`parameters`).

Measurement kinds and filename rules belong to project implementations. The package stores the
resulting `Symbol` without assigning project-specific meaning to it.

## Direct I-V data

DC current-voltage measurements use `:i` for current in amperes and `:v` for voltage in volts.
Voltage-driven sweeps, four-terminal measurements, breakdown measurements, and their visualizers
share these names so the same data can be reused across compatible views.
