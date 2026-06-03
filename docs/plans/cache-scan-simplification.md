# Cache and Scan Simplification

## Context

The intended behavior is simple. When the app opens a project, it should load whatever cached
measurements are available and show them quickly. At the same time, it should scan the filesystem in
the background. The cache gives the UI a fast first tree; the filesystem scan checks whether that
tree still matches the source files.

The current code does not model this as one project load. `_open_project_path!` starts
`_launch_project_reload_job!`, which tries to read the HDF5 cache through `load_project_cache`.
`scan_source` is a separate path, usually reached through manual scan or cache-build actions. Both
paths can produce a full `MeasurementHierarchy`, and the UI applies whichever result is active
instead of treating cache and filesystem data as two observations of the same project.

The rule for correctness is file-based. A cached file entry is usable only when the source
fingerprint stored in that entry matches the current filesystem fingerprint for the same normalized
path. If the source file exists and the fingerprint differs, the cache entry must be rebuilt. If the
source file exists and no cache entry exists, it must be added. If the cache entry exists and the
source file is gone, it must be removed. Cache write time can help explain unusual cases, but it
does not prove the cache is valid.

Projects should not know about the cache. A project should describe how to parse source files,
produce measurements, and load/analyze/draw plots. The cache layer should store and reuse project
results without requiring project-specific cache files such as `projects/RuO2/Cache.jl`.
No project hook should write HDF5 groups, declare cache schemas, compare cache freshness, or decide
cache repair policy.

The opened `AbstractProject` should carry package-managed runtime state: root path, cache identity,
indexes, loaded source data, and loaded measurement data. Project scripts should receive the project
object and semantic package accessors, not a separate cache object.

The project API names for this are `index_source_file` for measurement discovery,
`load_source_data` for source-table loading, and `data_of_measurements` for package-owned data
retrieval after measurements exist.

This document is the top-level refactor outline. Each item should be expanded only after the parent
level is agreed. If an item becomes too large for this page, its details move into a separate file.

## Plan

1. Verify current behavior and remove stale docs
2. Measure cache and scan costs
3. Define one record per source file and one record per cache entry
4. Define reconciliation rules
5. Remove project-specific cache code
6. Redesign cache layout and read path
7. Build unified project loading
8. Implement surgical cache repair
9. Simplify UI status and controls

## Details Files

No detail files yet.
