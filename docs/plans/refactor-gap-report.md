# Refactor Gap Report

## Purpose

This report lists the remaining differences between the current branch and the architecture in
[ARCHITECTURE.md](../ARCHITECTURE.md) and [workspace-vision.md](workspace-vision.md).

## Completed Foundation

The package is now divided into internal `Project`, `MeasurementIndex`, `Cache`, `Workspace`,
`Visualization`, and `Browser` modules. `DeviceParser.jl` and the root `ProjectCache.jl` no longer
exist.

One `Workspace` owns the open project, source root, measurement index, selection, cache state,
direct and processed data held in memory, scan work, cache work, errors, and background tasks.
Opening a workspace starts one source scan and one cache load. Both can progressively populate the
same measurement index.

Direct measurement data follows one path:

```text
workspace memory
  -> valid cache entry
  -> project source reader
```

Processed measurement data follows the same path and is cached separately. Debug plots bypass both
memory and cache. Cache repair uses the completed source scan rather than launching another scan.

The browser has typed state for its controls, filters, plot views, persistence, and rendering.
Normal plots call `setup_plot` and `plot_data!` directly; the old plot-job layer has been removed.

## Remaining Work

### External Projects

RuO2 and TASE are still compiled inside MeasurementBrowser. They should become ordinary external
Julia project code that imports the public MeasurementBrowser API. They do not need to be packages.
Tests should define small test projects instead of depending on bundled projects for package
contracts.

### Public API

The root export list is smaller, but the project API still needs a final review after projects move
out of the package. Internal scan and cache operations must remain internal. Public names should
describe only the work a project author or Julia user performs.

### Visualization

The current API supports project-defined `PlotKind` types through `setup_plot` and `plot_data!`.
Generic table views, column plots, composition, figure annotations, and live matching rules do not
exist yet. These should use the same operations from the GUI and Julia.

### Workflows

Figure scripts remain only as deprecated code. There is no saved workflow model yet. Workflows must
eventually store selections, processing, visualizers, figure layout, annotations, and live matching
rules without exposing browser widget state.

### Live Source Updates

Workspace startup scans once. Continuous source watching is not implemented. A watcher should feed
new, changed, and deleted files through the existing measurement-index and cache-repair paths rather
than create a second architecture.

### Analysis Progress

File interpretation is progressive and failure-tolerant. Cross-file measurement statistics still
run after file interpretation because they require the complete measurement set. This phase needs
clear progress and should preserve every successful measurement when individual analyses fail.

### Annotations

Tags are connected to the tree and measurement browser. Coordinates, saved spatial layout, and
notes have storage APIs but no browser interface. The spatial browser plan owns that remaining work.

### Code Surface

RuO2 plotting and analysis remain the largest project-specific code surface. `FigureScripts.jl`
also remains large despite being deprecated. Each replacement should delete the code it supersedes;
no compatibility paths should remain.

### Documentation And Function Quality

Non-trivial functions still exist without useful docstrings, especially in older project code.
Private helpers should be removed when direct code is clearer. The target is fewer functions and
less code, not the same logic compressed into fewer lines.

### Performance Verification

Fixture tests verify behavior but not large-project performance. Before merging, measure:

```text
cache startup
source scan
cross-file analysis
one-file cache repair
cached direct-data reads
cached processed-data reads
shutdown during active work
```

## Work Order

1. Finish the workspace transition and remove remaining browser-owned project state.
2. Move RuO2 and TASE outside the core package.
3. Define generic table and column visualizers.
4. Replace figure scripts with workflows and delete the deprecated implementation.
5. Add continuous source watching through the existing workspace operations.
6. Add the spatial browser and remaining annotation interfaces.
7. Measure the real RuO2 and TASE projects and remove demonstrated bottlenecks.
8. Continue the docstring, helper, function-count, and code-line review.
