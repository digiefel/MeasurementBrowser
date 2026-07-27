"""
The registration dialect: `define_project`, `register_item!`, `register_collection_analysis!`.

This is a convenience layer written over ordinary `(data, metadata)` values, and nothing more. It
holds no privileged position in the engine: its recipes become methods of the same typed stage
contract in `DataBrowserAPI` that any project implements, so a registered project and an equivalent
typed project traverse the same stages, cache on the same terms, and rehydrate the same way.

Being built purely on the public type API makes this package a living conformance test of that
extension surface — the role `DataBrowserPlots` plays for the GUI extension surface. If something
here needs a private hook, the type API is missing something.
"""
module DataBrowserRecipes

using DataBrowserAPI
using DataBrowserAPI:
    AbstractCollection,
    AbstractDataItem,
    AbstractDataSource,
    AbstractDataSourceItem,
    AbstractProject,
    @timed_dbg
import DataBrowserAPI:
    analyze,
    attach_record,
    entries,
    _has_collection_analysis,
    _has_collection_process,
    collection,
    id,
    item_data,
    item_type,
    kind,
    label,
    metadata,
    process,
    project_description,
    project_name,
    read,
    reconstruct
using DataBrowserAPI.ItemIndex: ItemRecord, MetadataDict, NamedCollection, metadata_dict
import DataBrowserAPI.ItemIndex: named_collection_path

include("carriers.jl")
include("recipes.jl")
include("registration.jl")
include("stages.jl")

export CollectionRecipe,
    ItemRecipe,
    NoMatch,
    Project,
    RegisteredDataItem,
    RegisteredReadResult,
    define_project,
    register_collection_analysis!,
    register_item!

end
