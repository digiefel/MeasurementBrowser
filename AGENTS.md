# Repository Guidelines

## Product Goal
MeasurementBrowser should be a better working environment than opening a Julia REPL and scripting
the same project by hand. Without the app, a user would find files, load them, parse them, reshape
or extract useful tables, compute values, and write figure code. The app exists to make that loop
faster, easier to inspect, and less repetitive while preserving the same mental model.

Project code should therefore stay close to the script a project maintainer would naturally write:
identify the measurements in a source file, load the data for a measurement, and render or summarize
that data. Project code should not know about cache files, background jobs, UI state, or package
lifecycle details. The package owns those concerns so that project implementations remain small and
focused on data interpretation and presentation.

The long-term direction is a live measurement workspace with instant feedback. Source data should be
watched continuously; new or changed files should appear without restarting the app; views should be
able to stay attached to selections, devices, or matching rules as the project changes. Built-in
visual inspection tools should cover common experimental tasks, and custom project code should extend
those tools without replacing the package's shared browsing, caching, selection, composition, and
reload behavior.

The project boundary should stay simple enough to support hot reloading, interpretation, or
out-of-process implementations later. Avoid designs that require project authors to understand the
Julia package internals, compile-time wiring, cache storage, worker orchestration, or GUI machinery.
The package should compete with the direct scripting loop by being faster and more interactive, not
by asking the project author to adopt a larger framework.

## Read The Docs First
Treat `docs/ARCHITECTURE.md` as the entry point before touching code. It links to the focused docs
for the data model, cache, GUI, storage, annotations, and figure scripts. Current-state docs describe
how the app works now; `docs/plans/` describes intended changes. Do not copy structure summaries
into this file when they belong in the docs.

When a change affects documented structure, update the matching doc in the same commit. If the code
and docs disagree, fix the disagreement instead of adding caveats or transition notes.

## Working Rules
Use Julia 1.12 style: 4-space indentation, `snake_case` functions and variables, `UpperCamelCase`
types, and concise public docstrings with signatures and return types. Prefer small, explicit APIs.
Project-facing functions should not expose cache controls, background-job controls, UI state, or
other package machinery.

Find simple, stable, idiomatic solutions. Refactor first when that makes the change smaller or more
obvious. Do not keep compatibility paths, fallback behavior, or legacy symbols unless the user
explicitly asks for them. Do not hide failures that should be fixed; surface errors clearly.

Use `rg` for search. Use `apply_patch` for manual edits. Do not revert user changes or unrelated
dirty files.

## Commands
Install the root package with `julia --project -e 'using Pkg; instantiate()'`. If you change a
path-dependency package, also instantiate that package's environment:
`src/DataAnalysis` or `src/Annotations`. `src/DataLoader` still exists for old generic CSV readers
but is deprecated; do not add new project readers there.

Launch the browser with an explicit source-data root:
`julia --project start.jl /path/to/measurements`. Do not rely on the default path in `start.jl`;
it is a local convenience path.

Run the full test suite with `julia --project -e 'using Pkg; Pkg.test()'`.

Tests live in `test/`. Keep fixtures small and deterministic. For GUI or plotting changes, test the
data, labels, and figure creation behavior that can be checked reliably without pixel-perfect image
assertions.

Commits should be small, cohesive, and imperative. Mention behavior changes in the commit body when
needed, and include the tests run.
