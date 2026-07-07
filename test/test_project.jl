using CSV
using DataFrames: DataFrame
using MeasurementBrowser

if !isdefined(Main, :TEST_CACHE_DEPOT)
    const TEST_CACHE_DEPOT = mktempdir()
    pushfirst!(DEPOT_PATH, TEST_CACHE_DEPOT)
    atexit(() -> rm(TEST_CACHE_DEPOT; force=true, recursive=true))
end

"""
A small registry project used across the package tests. One `:table` item per CSV; filenames
beginning with `broken` fail their read so failure handling can be exercised.
"""
function _build_test_project()
    project = define_project("TestProject"; description="Small project used to test package behavior")
    register_item!(
        project,
        :table;
        detect=file -> endswith(lowercase(file.filename), ".csv"),
        read=function (file)
            startswith(file.filename, "broken") && error("Broken test source file")
            return CSV.read(file.filepath, DataFrame)
        end,
        entries=function (file, data)
            name = splitext(file.filename)[1]
            return [DataItem(
                kind=:table,
                collection=["test", name],
                label="Table $name",
                data=data,
            )]
        end,
        process=function (item)
            data = item.data
            processed = DataFrame(data)
            processed.processed = data.x .+ data.y
            return DataItem(item, processed)
        end,
    )
    return project
end

const TEST_PROJECT = _build_test_project()

"""Build the directory source used by high-level project tests."""
function test_source(_project::Project, root_path::AbstractString)
    return MeasurementBrowser.DirectorySource(root_path)
end

"""Write one small source table for scan and cache tests."""
function write_test_source(path::AbstractString, offset::Real=0)::String
    write(path, "x,y\n$(offset + 1),$(offset + 2)\n$(offset + 3),$(offset + 4)\n")
    return String(path)
end

"""Block until source and graph work settle, then return the workspace."""
wait_workspace_idle!(workspace; timeout::Real=15) =
    MeasurementBrowser.Workspace.wait_workspace_idle!(workspace; timeout)
