# API reference

This page summarizes the public project interface. The guide chapters explain the meaning and timing
of each operation.

## Project construction

```julia
define_project(name::AbstractString; description::AbstractString="")::Project
```

Creates an empty project definition.

```julia
register_item!(project, [registration_name];
    read,
    detect=nothing,
    entries=nothing,
    label=nothing,
    collection=nothing,
    id=nothing,
    process=nothing,
    analyze=nothing,
)::Project
```

Registers or replaces one item pipeline. See [Registration API](registration.md).

```julia
register_collection_analysis!(project, registration_name;
    process=nothing,
    analyze=nothing,
)::Project
```

Registers or replaces collection-level operations for one registered pipeline.

## Workspace lifecycle

```julia
open_workspace(project, source)::Workspace
open_workspace(project, root_path::AbstractString)::Workspace
open_browser(workspace)
close_workspace!(workspace)::Nothing
wait_workspace_idle!(workspace; timeout=60)::Workspace
workspace_status(workspace)::WorkspaceStatus
```

See [Workspaces](workspaces.md).

## Workspace items

```julia
query_items(workspace)::Vector{String}
query_items(workspace, predicate::AbstractString)::Vector{String}
select_items!(workspace, ids)::Nothing
materialize_items(workspace, ids)::Vector{<:AbstractDataItem}
materialize_items(workspace)::Vector{<:AbstractDataItem}
read_item_data(workspace, ids)::Vector
read_item_data(workspace)::Vector
```

The no-argument `query_items` returns every currently published stable item ID. The no-argument
materialization and data-access forms use the current selection. Predicate queries filter committed
effective metadata in the project cache.

## Registration callbacks

```julia
detect = (file::SourceFile) -> accepted::Bool
read = (file::SourceFile) -> loaded_data::LoadedData
entries = (loaded_data::LoadedData, metadata::Dict) -> items::Vector
label = (data::ItemData, metadata::Dict) -> label::String
collection = (data::ItemData, metadata::Dict) -> path::Vector{String}
id = (data::ItemData, metadata::Dict) -> key
process = (data::ItemData, metadata::Dict) -> processed_data::ProcessedData
analyze = (processed_data::ProcessedData, metadata::Dict) -> additional_metadata::Dict
```

## Item interface

```julia
item_data(item::MyItem)::MyData
metadata(item::MyItem)::Dict
item_label(item::MyItem)::String
collection(item::MyItem)::Vector{<:AbstractCollection}
label(collection::AbstractCollection)::String
metadata(collection::AbstractCollection)::Dict
id(item::MyItem)::Any
process(item::MyItem)::MyProcessedItem
analyze(item::MyProcessedItem)::Dict
fingerprint(item::MyItem)::Any
cacheable(item::MyItem)::Bool
```

See [Type API](type-api.md).

## Source interface

```julia
source_id(source::MySource)::String
source_label(source::MySource)::String
source_items(source::MySource)::Vector{MySourceItem}
open_source(source::MySource)::MySource
close_source!(source::MySource)::Nothing
watch_source(source::MySource, on_change; cancel_token)::Nothing
source_open_options(source::MySource)::NamedTuple
```

## Source-item interface

```julia
source_item_id(item::MySourceItem)::String
source_item_label(item::MySourceItem)::String
fingerprint(item::MySourceItem)::Any
source_item_path(item::MySourceItem)::Union{Nothing,String}
source_item_timestamp(item::MySourceItem)::Any
metadata(item::MySourceItem)::Dict

data_items(
    project,
    source::MySource,
    source_item::MySourceItem,
)::Vector{<:AbstractDataItem}
```
