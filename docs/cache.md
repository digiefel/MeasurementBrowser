# HDF5 Cache

The HDF5 cache exists because the source tree can be slow to scan. Opening a project should not have
to wait for every CSV to be discovered, parsed, interpreted, and analyzed before the browser can show
useful measurements.

The cache is the fast path. It stores enough generated project data to repopulate the browser tree
and serve normal plot jobs without reparsing every CSV on startup. The CSV files remain the source of
truth. A cached file entry is usable only when its stored source fingerprint still matches the
current file on disk.

The intended user experience is: open a project, show cached measurements quickly when they exist,
then check the filesystem in the background and repair stale cache entries. The current
implementation supports the cache read and cache update pieces, but cache loading and source
scanning are still separate UI/code paths.

The cache is not source-root metadata and it is not meant to be hand-edited.

Current implementation lives mostly in [src/ProjectCache.jl](../src/ProjectCache.jl). RuO2-specific
payload storage lives in [src/projects/RuO2/Cache.jl](../src/projects/RuO2/Cache.jl).

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
| `/meta` | project name, cache id, root path, device metadata flag, update time |
| `/semantics` | project-defined semantic field names |
| `/files` | one group per cached source CSV |
| `/indexes` | compact startup index rebuilt after cache writes |

Each `/files/<file_key>` group stores:

| Entry | Purpose |
|---|---|
| `fingerprint` fields | source path, size, mtime |
| `source_file_unique_id` | source file id used during indexing |
| `measurement_keys` | ordered keys for measurement groups |
| `measurements/<measurement_key>` | serialized `MeasurementInfo` metadata |
| `status` | `ok`, `skipped`, or `error` |
| project payloads | project-specific cached plot data |

Measurement groups still store metadata next to project payloads, but startup cache loading does not
walk those groups.

`/indexes` stores the browser metadata as compact arrays: source file statuses, file errors,
measurement ids, filenames, filepaths, clean titles, kinds, timestamps, device paths, device
parameters, measurement parameters, and measurement stats. This is the data used to rebuild the
browser tree when a project opens.

## Loading the Cache

`load_project_cache(identity)` reads cached browser metadata without walking the source CSV tree and
without walking every cached file group.

It currently reads:

1. `/meta` for validation.
2. `/semantics`.
3. `/meta/has_device_metadata`.
4. `/indexes` compact arrays.

It returns a `ProjectCacheSnapshot`, containing:

| Field | Meaning |
|---|---|
| `identity` | cache identity used for the read |
| `hierarchy` | `MeasurementHierarchy` rebuilt from cached measurements |
| `status` | cache status counts |
| `semantic_fields` | project semantic fields stored in cache |
| `errors` | per-file cache transform errors |

Startup cache load does not read dataframe plot payloads. Plot payloads are read later when plot jobs
request them.

Caches without the current compact index are considered out of date and must be updated or rebuilt.

## Writing and Updating the Cache

`write_project_cache!(identity, source; mode)` writes cache entries from a `SourceScan`.

Supported modes:

| Mode | Behavior |
|---|---|
| `:build` | create a cache from the source scan |
| `:rebuild` | replace the existing cache contents |
| `:update` | rewrite only new or stale source files and remove deleted files |

During write, each source file is compared with cached fingerprints. Unchanged files are kept.
Changed or new files are rewritten. Cached files missing from the source scan are deleted.

After file updates, `/indexes`, `/meta`, and `/semantics` are rebuilt. The index is the cache
startup path; per-file groups are for payload lookup and cache repair.

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

When a cache is loaded without a source scan, the app can count cached entries but cannot know which
source files are fresh, stale, new, or deleted. That requires filesystem scanning.

## Plot Payloads

Normal single-measurement plot jobs require an active cache identity and read cached payloads through
`_measurement_group_for_cached_plot`.

RuO2 regular measurements store analyzed plot payloads under each measurement group:

| Entry | Purpose |
|---|---|
| `plot/source = "analyzed"` | marks a normal analyzed payload |
| `plot/df` | cached dataframe |
| `plot/scalars` | scalar analysis fields |
| `plot/arrays` | vector analysis fields |

RuO2 PUND fatigue files are stored differently. The full fatigue dataframe is cached once at the
file level under `signals/pund_fatigue`, and each virtual measurement stores a lightweight marker.
When plotting one cycle, the cached full dataframe is read and the requested cycle is selected.

Debug plot mode bypasses cached plot payloads and reads source CSV data directly.

## UI Behavior

Opening a project currently starts cache reload first through `_open_project_path!` and
`_launch_project_reload_job!`. On success, the UI applies the cached `MeasurementHierarchy`. Source
scanning is a separate path triggered by scan/rescan or cache build/update actions.

Because cache loading and source scanning are separate paths today, the UI can show cached
measurements before checking whether the source tree still matches the cache.
