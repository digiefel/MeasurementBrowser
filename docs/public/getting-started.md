# Your first project

The smallest project treats every discovered file as one item. DataBrowser supplies source
metadata, a source-derived label, a directory-derived collection, and an identity.

The complete pipeline is therefore:

```text
SourceFile -> read -> one data item
```

For a directory of CSV files, `read` receives one `SourceFile` and returns one `DataFrame`:

```julia
using CSV
using DataBrowser
using DataFrames: DataFrame

project = define_project("CSV tables")

register_item!(project;
    read = (file::SourceFile) -> CSV.read(file.filepath, DataFrame),
)

workspace = open_workspace(project, "/path/to/csv-files")
open_browser(workspace)
```

There is no item type to define and no metadata container to construct. `SourceFile` metadata
includes `:filename` and, when available, `:timestamp`. A returned
`Vector{Float64}`, `DataFrame`, image, or domain object is one item's data.

## What DataBrowser does with this definition

For every discovered file, DataBrowser:

1. calls `read` once;
2. treats the returned value as one item;
3. publishes the item using its source filename and directory;
4. processes requests in the background;
5. persists reusable workspace state and results.

The project code only explains how to read the data.

## Adding metadata during reading

When reading naturally reveals source-wide metadata, return it with the data:

```julia
register_item!(project;
    read = (file::SourceFile) -> (
        data=CSV.read(file.filepath, DataFrame),
        metadata=Dict(:acquired_at => file.timestamp),
    ),
)
```

The outer named tuple tells DataBrowser which value is data and which `Dict` is metadata. The
`DataFrame` remains the item's data; it is not converted to a public wrapper type.

## When one file contains several items

Add `entries` when users should select parts of a loaded source independently:

```julia
register_item!(project;
    read = (file::SourceFile) -> table::DataFrame,
    entries = (table::DataFrame, metadata::Dict) -> items::Vector,
)
```

`read` still runs once. `entries` receives its result and returns one value per item. The
[Registration API](registration.md) explains the accepted return shapes and every optional stage.
