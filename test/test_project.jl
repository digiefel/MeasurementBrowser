using CSV
using DataFrames: DataFrame
using MeasurementBrowser

import MeasurementBrowser:
    detect_kind,
    display_label,
    interpret_file,
    kind_label,
    load_source_data,
    parse_device_info,
    process_measurement_data,
    project_description,
    project_name

struct TestProject <: AbstractProject end

const TEST_PROJECT = TestProject()

project_name(::TestProject)::String = "TestProject"
project_description(::TestProject)::String = "Small project used to test package behavior"

"""Return one device path from the test filename."""
function parse_device_info(::TestProject, file::SourceFile)::DeviceInfo
    return DeviceInfo(["test", splitext(file.filename)[1]])
end

detect_kind(::TestProject, filename::String)::Symbol =
    endswith(lowercase(filename), ".csv") ? :table : :unknown

kind_label(::TestProject, kind::Symbol)::String = kind === :table ? "Table" : "Unknown"

display_label(::TestProject, measurement::MeasurementInfo)::String = measurement.clean_title

"""Interpret one test CSV, with `broken` filenames reserved for failure tests."""
function interpret_file(project::TestProject, file::SourceFile)::Vector{MeasurementInfo}
    startswith(file.filename, "broken") && error("Broken test source file")
    detect_kind(project, file.filename) === :unknown && return MeasurementInfo[]
    device = parse_device_info(project, file)
    return [MeasurementInfo(
        filepath=file.filepath,
        measurement_kind=:table,
        device_info=device,
        timestamp=file.timestamp,
        clean_title="Table $(last(device.location))",
    )]
end

"""Read one test CSV as its direct measurement table."""
function load_source_data(
    ::TestProject,
    source_file::SourceFile;
    measurement::Union{Nothing,MeasurementInfo}=nothing,
)::DataFrame
    return CSV.read(source_file.filepath, DataFrame)
end

"""Add one deterministic column used to verify processed-data caching."""
function process_measurement_data(
    workspace::MeasurementBrowser.Workspace.Workspace{TestProject},
    measurement::MeasurementInfo,
)::DataFrame
    data = only(read_measurement_data(workspace, [measurement]))
    processed = DataFrame(data)
    processed.processed = data.x .+ data.y
    return processed
end

"""Write one small source table for scan and cache tests."""
function write_test_source(path::AbstractString, offset::Real=0)::String
    write(path, "x,y\n$(offset + 1),$(offset + 2)\n$(offset + 3),$(offset + 4)\n")
    return String(path)
end
