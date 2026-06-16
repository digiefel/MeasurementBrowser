using CSV
using DataFrames: DataFrame
using MeasurementBrowser

"""
A small registry project used across the package tests. One `:table` measurement per CSV; filenames
beginning with `broken` fail their read so failure handling can be exercised.
"""
function _build_test_project()
    project = define_project("TestProject"; description="Small project used to test package behavior")
    register_measurement!(
        project,
        :table;
        detect=file -> endswith(lowercase(file.filename), ".csv"),
        read=function (file)
            startswith(file.filename, "broken") && error("Broken test source file")
            return CSV.read(file.filepath, DataFrame)
        end,
        measurements=function (file, _data)
            name = splitext(file.filename)[1]
            return [ItemRecord(
                filepath=file.filepath,
                kind=:table,
                collection=["test", name],
                timestamp=file.timestamp,
                clean_title="Table $name",
            )]
        end,
        process=function (_measurement, data)
            processed = DataFrame(data)
            processed.processed = data.x .+ data.y
            return processed
        end,
    )
    return project
end

const TEST_PROJECT = _build_test_project()

"""Write one small source table for scan and cache tests."""
function write_test_source(path::AbstractString, offset::Real=0)::String
    write(path, "x,y\n$(offset + 1),$(offset + 2)\n$(offset + 3),$(offset + 4)\n")
    return String(path)
end
