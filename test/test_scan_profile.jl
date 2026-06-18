using MeasurementBrowser
using CSV
using DataFrames: DataFrame, nrow
using Serialization: serialize, deserialize
using Test

const MB = MeasurementBrowser

"""Build a small registry project whose CSV reader expands one item per file."""
function _profile_project()
    project = MB.define_project("ProfileProject")
    MB.register_item!(
        project,
        :table;
        detect=file -> endswith(file.filename, ".csv"),
        read=file -> DataFrame(CSV.File(file.filepath)),
        entries=(file, data) -> [DataItem(
            kind=:table,
            collection=["dev", splitext(file.filename)[1]],
            label=file.filename,
        )],
        stats=(item, data) -> Dict{Symbol,Any}(:rows => nrow(data)),
    )
    return project
end

@testset "scan skip on unchanged files" begin
    mktempdir() do dir
        write_test_source(joinpath(dir, "first.csv"))
        write_test_source(joinpath(dir, "second.csv"), 10)

        source = scan_test_source(TEST_PROJECT, dir)

        # Nothing changed: the scan must short-circuit to the very same cached object.
        skipped = MB.scan_source(
            TEST_PROJECT,
            test_source(TEST_PROJECT, dir);
            cached_source=source,
        )
        @test skipped === source

        # A changed fingerprint forces a real scan that returns a fresh object.
        write_test_source(joinpath(dir, "second.csv"), 99)
        rescanned = MB.scan_source(
            TEST_PROJECT,
            test_source(TEST_PROJECT, dir);
            cached_source=source,
        )
        @test rescanned !== source
        @test length(rescanned.source_item_fingerprints) == 2
    end
end

@testset "per-kind scan profile" begin
    mktempdir() do dir
        for i in 1:4
            write(joinpath(dir, "m$i.csv"), "x,y\n1,2\n3,4\n")
        end

        project = _profile_project()

        source = scan_test_source(project, dir)
        rows = MB.scan_profile_summary(project)
        @test length(rows) == 1
        row = only(rows)
        @test row.kind == :table
        @test row.source_items == 4
        @test row.items == 4
        @test row.read_seconds >= 0
        @test row.stats_seconds >= 0

        # A cache hit does no per-kind work, so the profile resets to empty.
        MB.scan_source(
            project,
            test_source(project, dir);
            cached_source=source,
        )
        @test isempty(MB.scan_profile_summary(project))
    end
end

@testset "Project serialization drops transient state" begin
    project = _profile_project()
    # Dirty the transient fields the cache must never persist.
    push!(project.stat_failures, ("/tmp/x.csv", "id", "boom"))
    project.scan_profile[:table] = MB.KindProfile(2, 5, 1.0, 0.5)

    io = IOBuffer()
    serialize(io, project)
    seekstart(io)
    restored = deserialize(io)

    @test restored isa MB.Project
    @test restored.name == project.name
    @test length(restored.recipes) == 1
    @test restored.recipes[1].kind == :table
    # Transient state is rebuilt empty rather than carried across the cache boundary.
    @test isempty(restored.stat_failures)
    @test isempty(restored.scan_profile)

    # Project serialization should preserve shared references.
    shared = IOBuffer()
    serialize(shared, (project, project))
    seekstart(shared)
    a, b = deserialize(shared)
    @test a === b
end

@testset "per-item stat failure surfaces through the source scan" begin
    mktempdir() do dir
        write(joinpath(dir, "ok.csv"), "x,y\n1,2\n")
        write(joinpath(dir, "bad.csv"), "x,y\n1,2\n")

        project = MB.define_project("FailureProject")
        MB.register_item!(
            project,
            :table;
            detect=file -> endswith(file.filename, ".csv"),
            read=file -> DataFrame(CSV.File(file.filepath)),
            entries=(file, data) -> [DataItem(
                kind=:table,
                collection=["dev", splitext(file.filename)[1]],
                label=file.filename,
            )],
            # Stats throw only for bad.csv; ok.csv still gets its stats computed.
            stats=function (item, data)
                startswith(item.label, "bad") && error("stats blew up")
                return Dict{Symbol,Any}(:rows => nrow(data))
            end,
        )

        source = scan_test_source(project, dir)
        @test length(source.analysis_failures) == 1
        failure = only(source.analysis_failures)
        @test basename(failure.source_item_id) == "bad.csv"
        @test occursin("step=stats", failure.message)

        ok = only(item for item in source.hierarchy.all_items if endswith(item.source_item_path, "ok.csv"))
        @test ok.stats[:rows] == 1
    end
end
