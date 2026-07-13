# Registration API

`register_item!` connects ordinary Julia functions to the DataBrowser pipeline. Only `read` is
required. Add other callbacks when a source contains several items, needs processing, or needs
project-specific description.

Read [How DataBrowser works](pipeline.md) first. This chapter specifies the callback interface for
that pipeline.

## Complete callback shape

The complete registration has this data flow:

```julia
register_item!(project, registration_name;
    detect = (file::SourceFile) -> accepted::Bool,
    read = (file::SourceFile) -> loaded_data::LoadedData,
    entries = (loaded_data::LoadedData, metadata::Dict) -> items::Vector,
    label = (data::ItemData, metadata::Dict) -> label::String,
    collection = (data::ItemData, metadata::Dict) -> path::Vector{String},
    id = (data::ItemData, metadata::Dict) -> key,
    process = (data::ItemData, metadata::Dict) -> processed_data::ProcessedData,
    analyze = (processed_data::ProcessedData, metadata::Dict) -> additional_metadata::Dict,
)
```

`LoadedData`, `ItemData`, and `ProcessedData` stand for concrete types chosen by the project. They
are not DataBrowser types, and each stage may use a different type.

The registration name is optional. It identifies a registration inside the project; it is not a
property attached to every item.

## Returning data with metadata

`read` and each value returned by `entries` may return either:

- an ordinary Julia value, interpreted as data with no new metadata; or
- `(data=value, metadata=Dict(...))`.

For example, `entries` can conceptually return two channels as:

```julia
[
    (data=channel_a, metadata=Dict(:channel => "A")),
    (data=channel_b, metadata=Dict(:channel => "B")),
]
```

The outer named tuple has a fixed role. Its `data` value continues through the pipeline and its
`metadata` dictionary is merged with metadata already known for the item. DataBrowser keeps its
identity, scheduling, and cache records private.

## `read`

```julia
read = (file::SourceFile) -> loaded_data::LoadedData
```

`read` loads one accepted source. It runs once for that source revision. Parsing a complete file,
reading a shared header, decompression, and calculations required by every item from the source
belong here.

When `entries` is omitted, the result of `read` is one item. A vector is not implicitly
expanded: `Vector{Float64}` remains one item's data.

Return `(data=..., metadata=Dict(...))` when reading also discovers metadata shared by every item in the
source. That metadata is supplied to `entries` and inherited by its returned items.

## `entries`

```julia
entries = (loaded_data::LoadedData, metadata::Dict) -> items::Vector
```

`entries` interprets a loaded source as zero or more items. Each returned element is ordinary
item data or `(data=..., metadata=Dict(...))`.

This stage may be expensive. It runs once after `read` and should compute item metadata during the
same interpretation when that avoids another pass over the loaded source.

Omit `entries` when one source produces one item. Return an empty vector when an accepted source
contains no usable items.

## `label`

```julia
label = (data::ItemData, metadata::Dict) -> label::String
```

`label` returns the text shown for an item. It runs while the item is published and should be cheap.
When omitted, DataBrowser derives a label from the source.

## `collection`

```julia
collection = (data::ItemData, metadata::Dict) -> path::Vector{String}
```

`collection` returns the successive groups containing an item, from broadest to most specific. It
runs while the item is published and should be cheap.

For example, `["Sample A", "Device 3"]` places an item under `Sample A → Device 3`. Directory-backed
projects default to the source file's directory relative to the workspace root.

See [Metadata and collections](metadata-and-collections.md) for the full grouping model.

## `id`

```julia
id = (data::ItemData, metadata::Dict) -> key
```

DataBrowser generates an integer sibling key when `id` is omitted. Provide a channel name, cycle
number, database key, or another stable domain value only when identity must survive inserted or
reordered siblings.

## `process`

```julia
process = (data::ItemData, metadata::Dict) -> processed_data::ProcessedData
```

`process` transforms one item's data into the value consumed by inspectors, visualizers, and
`analyze`. It is optional; without it, the original data continues unchanged.

Processing is deferred and can be cached. Cleaning, normalization, conversion to a useful
representation, and other repeatable per-item transformations belong here.

## `analyze`

```julia
analyze = (processed_data::ProcessedData, metadata::Dict) -> additional_metadata::Dict
```

`analyze` derives searchable metadata from processed data. Its dictionary is merged into the item's
metadata and does not replace the processed value.

Analysis is deferred and can be cached. Expensive statistics, quality checks, extracted parameters,
and summaries used for filtering or querying belong here.

## Detection

```julia
detect = (file::SourceFile) -> accepted::Bool
```

`detect` decides whether a registration handles a source. It should normally inspect information
already discovered with the source, such as its filename, and avoid opening the file.

Registrations are tried in registration order. The first match handles the source, so specific
detectors belong before general ones.

## Registration identity and repeated execution

A project with one interpretation normally uses the unnamed form:

```julia
register_item!(project;
    read = (file::SourceFile) -> loaded_data::LoadedData,
)
```

A project has one unnamed registration slot. Calling the unnamed form again replaces its previous
contents. Rerunning the project script therefore updates the registration rather than duplicating
it.

Name registrations when a project has several interpretations:

```julia
register_item!(project, :spectrum;
    detect = (file::SourceFile) -> accepted::Bool,
    read = (file::SourceFile) -> spectrum::SpectrumData,
)

register_item!(project, :image;
    detect = (file::SourceFile) -> accepted::Bool,
    read = (file::SourceFile) -> image::ImageData,
)
```

Calling `register_item!` again with the same name replaces that registration in place. This gives
Revise and repeated script execution stable update semantics. The name also connects collection
operations to the corresponding registered pipeline. Items do not carry a public `kind` symbol.

## `SourceFile`

Directory-backed workspaces discover files before registrations run. Each callback receives a
`SourceFile` with information already collected during discovery:

| Property | Meaning |
|---|---|
| `filepath` | complete path used to open the file |
| `filename` | final path component used for recognition and filename metadata |
| `timestamp` | discovered file timestamp, when available |
| `fingerprint` | value used to detect source changes |

The source object lets DataBrowser perform filesystem discovery once and keeps that work out of
project callbacks. Its metadata dictionary contains `:filename` and, when a timestamp was
discovered, `:timestamp`. This dictionary is passed to `entries` and the other item callbacks.

## Collection operations

Item callbacks operate on one item at a time. Collection operations perform work that depends on a
related group of items:

```julia
register_collection_analysis!(project, registration_name;
    process = (data::Vector, metadata::Vector{<:Dict}) -> processed_data::Vector,
    analyze = (processed_data::Vector, metadata::Vector{<:Dict}) -> collection_metadata::Dict,
)
```

Collection `process` receives members after their item-level processing and returns one data value
per member in the same order. Collection `analyze` returns metadata describing the collection
itself. Both wait for their required member results, run in the background, and can be cached.
