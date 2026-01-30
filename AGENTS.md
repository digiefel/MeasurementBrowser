# Repository Guidelines

## Project Structure & Modules
- Root `Project.toml` defines the MeasurementBrowser package; entry script is `start.jl`.
- `src/` houses core code: `MeasurementBrowser.jl` wires modules, `DeviceParser.jl` parses filenames, `Gui.jl` + `MakieIntegration.jl` drive the GLMakie/CImGui UI, and helper packages live in `src/DataLoader`, `src/DataPlotter`, `src/DataAnalysis` (each with its own `Project.toml`).
- Measurement data are CSV files; `start.jl` defaults to a local path—pass a directory argument to override.

## Setup, Build, Run
- Install dependencies: `julia --project -e 'using Pkg; instantiate()'` (repeat with `--project=src/DataLoader` etc. when hacking subpackages).
- Launch the browser: `julia --project start.jl /path/to/measurements`.
- Precompile without opening the UI: `julia --project -e 'using Pkg; precompile()'`.

## Coding Style & Naming
- Target Julia 1.12; 4-space indentation; aim for ~100-char lines.
- Use `snake_case` for functions/variables, `UpperCamelCase` for structs/types, `ALL_CAPS` for constants; keep side effects explicit.
- Prefer docstrings for public APIs and brief comments only where UI/GL hooks are non-obvious.
- Find simple, stable, idiomatic solutions; avoid shortcuts.
- Consider whether a cleaner solution emerges by refactoring first.
- Don’t catch errors; let failures surface.

## Testing Guidelines
- Tests live in `test/` with fixtures (sample CSVs) plus `runtests.jl` and feature-focused files.
- Run all tests: `julia --project -e 'using Pkg; test()'`.
- Add targeted `@testset`s per feature; keep fixtures small/deterministic and align filenames with parsing regexes.
- For plot/GUI changes, assert metadata/labels and ensure figure creation does not error rather than pixel-perfect checks.

## Commit Workflow
- Commits use short, imperative summaries like existing history (`fix negative filters`, `adjust dependencies`); keep scope tight and cohesive.
- Describe behavior changes in the commit body when needed and note test commands run; reference relevant files/modules.

## Notes & Gotchas
- UI depends on GLMakie + CImGui; ensure an OpenGL-capable environment when running locally.
- File parsing hinges on `RuO2test_...` patterns; update regexes and tests together if naming schemes shift.
- Scanning very large measurement trees is slow; while developing, point `start.jl` at a small sample folder or a copy of `test/` fixtures to shorten load times.
