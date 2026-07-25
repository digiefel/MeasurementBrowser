using DataBrowser

function _registered_collection_key(collections, names...)
    id_path = DataBrowserAPI.collection_id_path(
        DataBrowserAPI.ItemIndex.RegisteredCollection.(String.(names)))
    return collections.key_by_id[last(id_path)]
end

_registered_collection_record(collections, names...) =
    collections.records[_registered_collection_key(collections, names...)]

if !isdefined(Main, :TEST_CACHE_DEPOT)
    const TEST_CACHE_DEPOT = mktempdir()
    pushfirst!(DEPOT_PATH, TEST_CACHE_DEPOT)
    atexit(() -> rm(TEST_CACHE_DEPOT; force=true, recursive=true))
end

"""Build the directory source used by high-level project tests."""
function test_source(_project::Project, root_path::AbstractString)
    return DataBrowser.DirectorySource(root_path)
end

"""Block until source and graph work settle, then return the workspace."""
wait_workspace_idle!(workspace; timeout::Real=15) =
    DataBrowserCore.Workspace.wait_workspace_idle!(workspace; timeout)
