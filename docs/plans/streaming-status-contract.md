# Streaming Status Contract

## Purpose

`WorkspaceStatus` is the one watcher-facing summary of background work. The GUI, scripts, and future
workflow tools should not inspect scan jobs, processing queues, DuckDB state, or analysis tasks
directly.

```julia
struct WorkspaceStatus
    level::Symbol
    label::String
    detail::String
    busy::Bool
    progress::Union{Nothing,Float32}
    errors::Vector{WorkspaceError}
end
```

The workspace owns detailed mutable job state. `poll_workspace!` consumes completed events, publishes
replacement index/statistics values, and rebuilds this small snapshot only while work is active or
state changed. Idle frames reuse the previous value.

## User Model

One workspace update overlaps three understandable activities:

- **Scanning** discovers source items, compares fingerprints, and interprets changed items.
- **Caching and processing** commits interpreted data, produces missing processed data, and records
  deliberate cache opt-outs.
- **Summarizing** computes item and collection statistics.

The status label names the activity currently most useful to the user. The detail line carries counts
or the current source item. A determinate fraction is shown only when both completed and total work
are known. Counts never move backward.

## Publication Rules

- Cached hierarchy state may appear before source discovery finishes.
- Source records appear only after their source-item cache group commits.
- Cacheable processed values reach selected views only after their processed write commits.
- Deliberately uncached processed values reach current waiters after processing finishes.
- Processing and statistics failures are published before their job can appear idle.
- Results from a canceled or older workspace generation are ignored.
- A rescan builds replacement values without mutating a hierarchy still read by analysis or the GUI.

The defensive copy made before a cached hierarchy is handed to both refresh and GUI state remains in
place because it prevents a demonstrated threaded crash.

## Errors

Errors attach to the smallest useful identity: source item, logical item, or collection path. The
status snapshot contains the current list in stable order. A failure does not turn unrelated valid
records or cache stages into misses.

Generated cache corruption is repaired automatically. Source-identity conflicts and failures that
prevent the workspace from continuing remain visible because they require a real decision or code
fix.

## GUI Contract

The GUI uses only:

- `status.level` for color;
- `status.label` for the short activity name;
- `status.detail` for one concise explanation;
- `status.busy` for activity indication;
- `status.progress` for the optional progress bar;
- `status.errors` for the error list.

It does not infer progress from DuckDB row counts or reach into `scan`, `processing`, `analysis`, or
cache transaction internals. This keeps batching, worker counts, and storage changes behind the
workspace boundary.
