using CSV
using DataFrames: DataFrame
using MeasurementBrowser

"""
A small registry project used across the package tests. One `:table` measurement per CSV; filenames
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
        entries=function (file, _data)
            name = splitext(file.filename)[1]
            return [DataItem(
                kind=:table,
                collection=["test", name],
                label="Table $name",
            )]
        end,
        process=function (_item, data)
            processed = DataFrame(data)
            processed.processed = data.x .+ data.y
            return processed
        end,
    )
    return project
end

const TEST_PROJECT = _build_test_project()

"""Build the private callback-adapter source used by high-level project tests."""
function test_source(project::Project, root_path::AbstractString)
    return MeasurementBrowser.RegisteredProjectSource(project, root_path)
end

"""Scan a registered test project through the source-first engine API."""
function scan_test_source(project::Project, root_path::AbstractString; kwargs...)
    return MeasurementBrowser.scan_source(project, test_source(project, root_path); kwargs...)
end

"""Write one small source table for scan and cache tests."""
function write_test_source(path::AbstractString, offset::Real=0)::String
    write(path, "x,y\n$(offset + 1),$(offset + 2)\n$(offset + 3),$(offset + 4)\n")
    return String(path)
end
