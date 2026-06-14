using MeasurementBrowser
using CSV
using DataFrames: DataFrame, nrow
using Serialization: serialize, deserialize
using Test

const MB = MeasurementBrowser

"""Build a small registry project whose CSV reader expands one measurement per file."""
function _profile_project()
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
    return project
end

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

        project = _profile_project()

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

@testset "RegistryProject serialization drops transient state" begin
    project = _profile_project()
    # Dirty the transient fields the cache must never persist.
    project.read_cache["/tmp/x.csv"] = DataFrame(a=[1])
    project.scan_profile[:table] = MB.KindProfile(2, 5, 1.0, 0.5)

    io = IOBuffer()
    serialize(io, project)
    seekstart(io)
    restored = deserialize(io)

    @test restored isa MB.RegistryProject
    @test restored.name == project.name
    @test length(restored.recipes) == 1
    @test restored.recipes[1].kind == :table
    # Transient state is rebuilt empty rather than carried across the cache boundary.
    @test isempty(restored.read_cache)
    @test isempty(restored.scan_profile)

    # The same project is reachable twice from a cached scan (source.project and hierarchy.project);
    # serialization must dedup so both references resolve to one object after load.
    shared = IOBuffer()
    serialize(shared, (project, project))
    seekstart(shared)
    a, b = deserialize(shared)
    @test a === b
end
