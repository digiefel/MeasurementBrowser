# Project Cache

Source scanning is inherently slow: the package must find every source file, interpret its
measurements, and finish project analysis before the final measurement tree is known. The cache
stores the result of that work so the previous tree can appear immediately while a new scan checks
the filesystem.

The source files remain authoritative. The cache is generated package data, stored outside the
project, and never exposed to project code.

## Startup

Opening a project starts cache loading and source scanning together.

Cache loading reads one compact index containing the previous `SourceScan`. It can populate the
workspace index immediately with the previous hierarchy, parameters, stats, skipped-file count,
source fingerprints, and analysis failures.

The filesystem scan independently streams current measurements into the same workspace index and
eventually produces the authoritative `SourceScan`. The package compares that scan with the cached
index and updates the cache when files were added, changed, deleted, or produced a different
analysis result.

## Identity

Each source root has one deterministic cache id:

```text
<root-folder>-<path-digest>
```

Cache files live under:

```text
DEPOT_PATH[1]/measurementbrowser/cache/<project>/<cache-id>.h5
```

The cached identity records the cache id, project name, normalized source root, and cache path. A
cache is accepted only when all identity fields and the schema version match.

## Contents

The HDF5 file contains:

| Entry | Purpose |
|---|---|
| `schema_version` attribute | Rejects old cache layouts before deserialization |
| `/index` | One uncompressed serialized `ProjectCacheIndex` for fast startup |
| `/direct/<file>/<measurement>` | Lazily cached data returned by `read_item_data` |
| `/processed/<file>/<measurement>` | Lazily cached data returned by `process_item_data` |

`ProjectCacheIndex` contains the completed source scan, a path lookup for its `SourceFile` entries,
and analysis errors grouped by physical source file.

Direct and processed data are separate because they have different meanings and lifetimes. Both are
grouped by physical source file, so changing or deleting one file invalidates its payloads with two
group deletions regardless of how many logical measurements it contains.

## Freshness

A cached source file matches the filesystem only when its normalized path, size, and modification
time match. The final scan comparison reports:

| Count | Meaning |
|---|---|
| fresh | cached file matches the current scan and has no analysis error |
| stale | fingerprint or analysis result changed |
| new | present only in the current scan |
| deleted | present only in the cache |
| errors | current source files whose analysis failed |

Analysis failure does not remove measurements from the tree. It is stored as a valid cache
difference so one failed file cannot invalidate the rest of the project.

## Updates

A normal update uses the completed source scan already owned by the workspace. It never starts
another scan. It replaces the compact index, preserves payloads for unchanged source files, and
deletes payload groups for changed, deleted, or newly failing files. A rebuild replaces the complete
HDF5 file.

Index writing does not load measurement data. Direct and processed data are cached only when a
caller requests them.

## Measurement Data

`read_item_data(workspace, measurements)` returns one direct dataframe per measurement, in
the requested order. It checks workspace memory, then reads all available cached payloads with one
HDF5 open. Missing data is loaded through the workspace project's `load_source_data` method, kept in
memory, and written back together with one HDF5 open.

`process_item_data(workspace, measurements)` follows the same path for processed data.
Missing entries are produced by the project's
`process_item_data(workspace, measurement)` method and cached separately.

Before reading or writing a payload, the package checks the current source fingerprint against the
in-memory cache index. Stale payloads are never returned.

The loaded cache index belongs to its workspace. There is no process-wide active cache.

All HDF5 access uses one process-wide lock. Background cache updates therefore cannot overlap plot
reads or lazy payload writes.
