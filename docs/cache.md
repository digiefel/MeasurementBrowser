# HDF5 Cache

The HDF5 cache exists because the source tree can be slow to scan. Opening a project should not have
to wait for every CSV to be discovered, parsed, interpreted, and analyzed before the browser can show
useful measurements.

The cache is the fast path. It stores enough generated project data to repopulate the browser tree
without reparsing every CSV on startup. The CSV files remain the source of truth. A cached file entry
is usable only when its stored source fingerprint still matches the current file on disk.

The intended user experience is: open a project, show cached measurements quickly when they exist,
scan the filesystem in the background, and repair stale cache entries without user action.

The cache is not source-root metadata and it is not meant to be hand-edited.

Current implementation lives in [src/ProjectCache.jl](../src/ProjectCache.jl). Projects do not own
cache storage files.

## Location and Identity

Cache files live under the first Julia depot:

```text
DEPOT_PATH[1]/measurementbrowser/cache/<project>/<cache_id>.h5
```

`project_cache_id(root_path)` derives a deterministic id from the normalized project root: a slug of
the folder name plus a short SHA1 digest of the full path.

`project_cache_identity(cache_id, project, root_path)` binds together the cache id, project name,
normalized root path, and final cache file path. Cache operations check that the identity project
and root match the source or cache being used.

## Source Fingerprints

Each cached CSV entry stores a `FileFingerprint`:

```julia
struct FileFingerprint
    path::String
    size_bytes::Int64
    mtime_ns::Int64
end
```

The cache treats a source file as fresh only when the normalized path, size, and mtime match the
stored fingerprint. If any of those differ, the cache entry is stale for that file.

## Cache File Structure

Top-level groups:

| Path | Purpose |
|---|---|
| `/meta` | cache schema version, project name, cache id, root path, device metadata flag, update time |
| `/files` | one group per cached source CSV |
| `/indexes` | compact startup index rebuilt after cache writes |

Each `/files/<file_key>` group stores:

| Entry | Purpose |
|---|---|
| `fingerprint` fields | source path, size, mtime |
| `source_file_unique_id` | source file id used during indexing |
| `measurement_keys` | ordered keys for measurement groups |
| `measurements/<measurement_key>` | optional payload group for that logical measurement |
| `measurements/<measurement_key>/data` | optional cached dataframe for that logical measurement |
| `status` | `ok`, `skipped`, or `error` |

Startup cache loading does not walk the per-file measurement groups.

`analysis_error` is a file status for a source file that was found and indexed, but where one or
more logical measurements failed analysis. The measurements stay visible; the failure is reported
as cache/source status.

`/indexes/startup_blob` stores the browser startup data as one versioned serialized record:
measurements, source file statuses, file errors, and skipped-file count. Startup uses the blob to
rebuild the browser tree when a project opens.

The startup blob is the only cache copy of measurement metadata, including parameters and stats.
Per-file measurement groups exist so lazy dataframe payloads have a stable place to live.

## Loading the Cache

`load_project_cache(identity)` reads cached browser metadata without walking the source CSV tree and
without walking every cached file group.

It currently reads:

1. `/meta` for validation.
2. `/meta/has_device_metadata`.
3. `/indexes/startup_blob`.

It returns a `ProjectCacheSnapshot`, containing:

| Field | Meaning |
|---|---|
| `identity` | cache identity used for the read |
| `hierarchy` | `MeasurementHierarchy` rebuilt from cached measurements |
| `status` | cache status counts |
| `errors` | per-file cache transform errors |

Startup cache load does not read dataframe plot payloads. Plot payloads are read later when plot jobs
request them.

Caches without the current schema version and compact index are out of date. They should be rebuilt.
Semantic fields do not decide cache validity.

## Writing and Updating the Cache

`write_project_cache!(identity, source; mode)` writes cache entries from a `SourceScan`.

Update mode is for compatible cache files. A compatible cache has the current schema version and the
same project/root identity.

Supported modes:

| Mode | Behavior |
|---|---|
| `:build` | create a cache from the source scan |
| `:rebuild` | replace the existing cache contents |
| `:update` | rewrite only new or stale source files and remove deleted files |

During write, each source file is compared with cached fingerprints. Unchanged files are kept.
Changed or new files are rewritten. Cached files missing from the source scan are deleted.

After file updates, `/indexes` and `/meta` are rebuilt. The index is the cache startup path;
per-file groups are for fingerprint checks, cache repair, and lazy measurement-data payloads.

Cache build and update do not eagerly write every dataframe payload. They write the browser tree,
source fingerprints, measurement metadata, and status. Measurement data is cached when
`data_of_measurements` has to load it.

## Cache Status

`ProjectCacheStatus` stores count summaries:

| Field | Meaning |
|---|---|
| `total_files` | CSV files in the current source scan |
| `cached_files` | files currently represented in cache |
| `fresh_files` | source files whose cache fingerprint matches and did not error |
| `stale_files` | cached files whose source fingerprint changed |
| `new_files` | source files missing from cache |
| `deleted_files` | cached files missing from source |
| `error_files` | cached files with transform errors |

`error_files` includes cache transform errors and analysis errors found during the source scan.

When a cache is loaded without a source scan, the app can count cached entries but cannot know which
source files are fresh, stale, new, or deleted. That requires filesystem scanning.

## Measurement Data

Normal plotting code asks for data through `data_of_measurements(project, measurements)`. That is
where cached measurement data is used. Project code does not read cache storage directly.

The package returns cached dataframe data only when the cached file fingerprint still matches the
current source file. If the cache has no data, the file changed, or the cache entry errored,
`data_of_measurements` opens the source file through `load_source_data`.

When `data_of_measurements` loads source data for a measurement whose cache entry is current, it
stores the returned dataframe in that measurement's cache entry for later calls.

## UI Behavior

Opening a project starts cache loading and source scanning together. Cache loading may apply a
cached `MeasurementHierarchy` before the source scan finishes. Source scanning streams interpreted
measurements as files complete, then applies the final scan result after project stats are computed.

When the final source scan is available, the UI compares source fingerprints with cached
fingerprints. Missing, stale, deleted, or analysis-error cache entries queue a cache update
automatically.
