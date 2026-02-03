# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

```bash
# Install dependencies
julia --project -e 'using Pkg; Pkg.instantiate()'

# Launch the GUI browser
julia --project start.jl /path/to/measurements

# Run tests
julia --project -e 'using Pkg; Pkg.test()'

# Precompile without launching UI
julia --project -e 'using Pkg; Pkg.precompile()'
```

When modifying subpackages (DataLoader, DataAnalysis, DataPlotter), also instantiate their local projects: `julia --project=src/DataLoader -e 'using Pkg; Pkg.instantiate()'`.

## Architecture

MeasurementBrowser is a Julia GUI application for browsing/analyzing ferroelectric and semiconductor device measurements (RuO2 test chips). It combines GLMakie plotting with CImGui for the UI.

### Module Dependency Graph

```
MeasurementBrowser.jl (root module)
  ├── DeviceParser.jl     — filename parsing, device hierarchy, measurement detection
  ├── Gui.jl              — CImGui UI (tree panel, plot panel, combined plots)
  └── MakieIntegration.jl — GLMakie↔CImGui bridge (included by Gui.jl)

Gui.jl uses:
  ├── DataPlotter/        — figure_for_file(), figure_for_files(), plot functions
  │     ├── uses DataLoader/   — CSV readers (read_fe_pund, read_iv_sweep, etc.)
  │     └── uses DataAnalysis/ — analyze_pund, analyze_tlm_combined, etc.
```

### Data Flow

1. **Scan:** `scan_directory()` walks a folder tree, creating `MeasurementInfo` structs from CSV filenames via regex patterns (`RuO2test_CHIP_TYPE_GEOMETRY_...`). Some files expand into multiple virtual measurements (breakdown → per-device, PUND fatigue → per-cycle).

2. **Hierarchy:** Measurements organize into a `MeasurementHierarchy` tree (chip → layer → device type → geometry). Optional `device_info.txt` at root merges metadata (area, thickness) into devices.

3. **Plot:** Single file → `figure_for_file(path, kind)` routes to the appropriate reader + plotter. Multiple files → `figure_for_files(paths, combined_kind)` for combined analysis (TLM, fatigue).

### Key Types (DeviceParser.jl)

- `MeasurementInfo` — one measurement: filename, filepath, kind (`:pund`, `:iv`, `:tlm4p`, `:breakdown`, `:wakeup`, `:pund_fatigue`), timestamp, device_info, parameters
- `DeviceInfo` — device location path + metadata parameters dict
- `HierarchyNode` / `MeasurementHierarchy` — browsable tree structure

### Measurement Kinds and Their Modules

| Kind | Detected by filename | Reader (DataLoader) | Analyzer (DataAnalysis) | Plotter (DataPlotter) |
|------|---------------------|---------------------|------------------------|-----------------------|
| `:pund` | "fe pund" / "fepund" | `read_fe_pund` | `analyze_pund` | `plot_fe_pund` |
| `:pund_fatigue` | "pund_fatigue" / "pund fatigue" | `read_pund_fatigue_cycle` | `analyze_pund` (trough fallback) | expands to virtual `:pund` entries |
| `:iv` | "i_v sweep" / "iv sweep" | `read_iv_sweep` | — | `plot_iv_sweep_single` |
| `:tlm4p` | "tlm_4p" / "tlm" | `read_tlm_4p` | `analyze_tlm_combined` | `plot_tlm_4p`, `plot_tlm_combined` |
| `:breakdown` | "break" / "breakdown" | `read_iv_sweep` | `analyze_breakdown` | `plot_iv_sweep_single` |
| `:wakeup` | "wakeup" | `read_wakeup` | — | `plot_wakeup` |

### Combined Plot Types

Selectable in the GUI when multiple measurements are selected:
- `:tlm_analysis` — width-normalized resistance vs length
- `:tlm_temperature` — sheet resistance vs temperature
- `:pund_fatigue` — P-E curve evolution + remnant polarization vs fatigue cycles

### Virtual Measurement Expansion

Some single CSV files expand into multiple `MeasurementInfo` entries during scanning:
- **Breakdown files** → `expand_multi_device()` splits multi-device breakdown measurements
- **PUND Fatigue files** → `expand_pund_fatigue()` creates one virtual `:pund` entry per cycle, storing cycle number in `parameters[:fatigue_cycle]`

## Coding Conventions

- Julia 1.12; 4-space indentation; ~100-char lines
- `snake_case` functions/variables, `UpperCamelCase` types, `ALL_CAPS` constants
- Don't catch errors; let failures surface
- Docstrings for public APIs; brief comments only where non-obvious
- Short imperative commit messages matching existing history style
- File parsing depends on `RuO2test_...` regex patterns — update regexes and tests together

## Testing

Tests in `test/` with fixture CSVs. For plot/GUI changes, assert metadata/labels and ensure figure creation doesn't error rather than pixel-perfect checks.

## Notes

- UI requires OpenGL-capable environment (GLMakie + CImGui)
- Scanning large measurement trees is slow; point at a small subfolder or `test/` fixtures during development
- `detect_measurement_kind()` ordering matters: more specific patterns (e.g. "pund_fatigue") must come before general ones ("pund")
