using MeasurementBrowser
using CSV
using DataFrames: DataFrame, nrow
using Test

const MB = MeasurementBrowser

@testset "scan skip on unchanged files" begin
    mktempdir() do dir
        write_test_source(joinpath(dir, "first.csv"))
        write_test_source(joinpath(dir, "second.csv"), 10)

        source = MB.scan_source(dir; project=TEST_PROJECT)
        cached_files = Dict(file.filepath => file for file in source.files)

        # Nothing changed: the scan must short-circuit to the very same cached object.
        skipped = MB.scan_source(
            dir;
            project=TEST_PROJECT,
            cached_files=cached_files,
            cached_source=source,
        )
        @test skipped === source

        # A changed fingerprint forces a real scan that returns a fresh object.
        write_test_source(joinpath(dir, "second.csv"), 99)
        rescanned = MB.scan_source(
            dir;
            project=TEST_PROJECT,
            cached_files=cached_files,
            cached_source=source,
        )
        @test rescanned !== source
        @test length(rescanned.files) == 2
    end
end

@testset "per-kind scan profile" begin
    mktempdir() do dir
        for i in 1:4
            write(joinpath(dir, "m$i.csv"), "x,y\n1,2\n3,4\n")
        end

        project = MB.define_project("ProfileProject")
        MB.register_measurement!(
            project,
            :table;
            detect=file -> endswith(file.filename, ".csv"),
            read=file -> DataFrame(CSV.File(file.filepath)),
            measurements=(file, data) -> [MB.MeasurementInfo(
                filepath=file.filepath,
                measurement_kind=:table,
                device_info=MB.DeviceInfo(["dev", splitext(file.filename)[1]]),
                timestamp=file.timestamp,
                clean_title=file.filename,
            )],
            stats=(mi, data) -> Dict{Symbol,Any}(:rows => nrow(data)),
        )

        source = MB.scan_source(dir; project=project)
        rows = MB.scan_profile_summary(project)
        @test length(rows) == 1
        row = only(rows)
        @test row.kind == :table
        @test row.files == 4
        @test row.measurements == 4
        @test row.read_seconds >= 0
        @test row.stats_seconds >= 0

        # A cache hit does no per-kind work, so the profile resets to empty.
        cached_files = Dict(file.filepath => file for file in source.files)
        MB.scan_source(
            dir;
            project=project,
            cached_files=cached_files,
            cached_source=source,
        )
        @test isempty(MB.scan_profile_summary(project))
    end
end
