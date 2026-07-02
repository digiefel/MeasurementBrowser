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
            data=data,
        )],
        stats=item -> Dict{Symbol,Any}(:rows => nrow(item.data)),
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

@testset "scan cache compares effective parameters" begin
    mktempdir() do dir
        write_test_source(joinpath(dir, "first.csv"))
        write(joinpath(dir, "metadata.txt"), "collection_path,wafer\ntest,A\n")

        source = scan_test_source(TEST_PROJECT, dir)

        write(joinpath(dir, "metadata.txt"), "collection_path,wafer\ntest,A\nunused,B\n")
        same = MB.scan_source(
            TEST_PROJECT,
            test_source(TEST_PROJECT, dir);
            cached_source=source,
        )
        @test same === source

        write(joinpath(dir, "metadata.txt"), "collection_path,wafer\ntest,B\n")
        changed = MB.scan_source(
            TEST_PROJECT,
            test_source(TEST_PROJECT, dir);
            cached_source=source,
        )
        @test changed !== source
        @test only(changed.hierarchy.all_items).parameters[:wafer] == "B"
    end
end

@testset "per-kind scan profile" begin
    mktempdir() do dir
        for i in 1:4
            write(joinpath(dir, "m$i.csv"), "x,y\n1,2\n3,4\n")
        end
        mkpath(joinpath(dir, "sub"))
        write(joinpath(dir, "sub", "m5.csv"), "x,y\n1,2\n3,4\n")

        project = _profile_project()

        source = scan_test_source(project, dir)
        rows = MB.scan_profile_summary(project)
        @test length(rows) == 1
        row = only(rows)
        @test row.kind == :table
        @test row.source_items == 5
        @test row.items == 5
        @test row.detect_seconds >= 0
        @test row.read_seconds >= 0
        @test row.entries_seconds >= 0
        @test row.process_seconds >= 0
        @test row.stats_seconds >= 0
        @test row.total_seconds >= row.read_seconds

        source_rows = MB.scan_source_profile(project)
        @test length(source_rows) == 5
        @test all(source_row -> source_row.kind == :table, source_rows)
        @test all(source_row -> source_row.items == 1, source_rows)
        @test all(source_row -> !isempty(source_row.thread_ids), source_rows)
        @test Set(row.source_item_label for row in source_rows) ==
            Set([("m$i.csv" for i in 1:4)..., joinpath("sub", "m5.csv")])
        @test Set(row.source_item_path for row in source_rows) ==
            Set([(joinpath(dir, "m$i.csv") for i in 1:4)..., joinpath(dir, "sub", "m5.csv")])
        @test all(row -> !isabspath(row.source_item_label), source_rows)
        @test issorted(source_rows; by=row -> row.total_seconds, rev=true)

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
    project.scan_profile["source.csv"] = MB.Projects.SourceItemProfile("source.csv")

    io = IOBuffer()
    serialize(io, project)
    seekstart(io)
    restored = deserialize(io)

    @test restored isa MB.Project
    @test restored.name == project.name
    @test length(restored.recipes) == 1
    @test restored.recipes[1].kind == :table
    # Transient state is rebuilt empty rather than carried across the cache boundary.
    @test isempty(restored.scan_profile)

    # Project serialization should preserve shared references.
    shared = IOBuffer()
    serialize(shared, (project, project))
    seekstart(shared)
    a, b = deserialize(shared)
    @test a === b
end

@testset "workspace stats see effective parameters" begin
    mktempdir() do dir
        write(joinpath(dir, "ok.csv"), "x,y\n1,2\n")
        write(joinpath(dir, "metadata.txt"), "collection_path,area_um2\ndev,42\n")

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
                data=data,
            )],
            stats=function (item)
                return Dict{Symbol,Any}(
                    :rows => nrow(item.data),
                    :area_um2 => item.parameters[:area_um2],
                )
            end,
        )

        workspace = MB.open_workspace(project, test_source(project, dir))
        try
            wait_workspace_idle!(workspace)

            @test workspace.scan.state == :done
            @test isempty(workspace.index.analysis_errors)
            ok = only(workspace.index.hierarchy.all_items)
            @test workspace.index.item_stats[ok.id][:rows] == 1
            @test workspace.index.item_stats[ok.id][:area_um2] == 42
        finally
            MB.close_workspace!(workspace)
        end
    end
end

@testset "expanded source is read once for scan, stats, and cache" begin
    mktempdir() do dir
        path = joinpath(dir, "expanded.csv")
        CSV.write(path, DataFrame(
            part=repeat(1:100; inner=10),
            x=collect(1.0:1000.0),
            y=collect(1001.0:2000.0),
        ))

        read_count = Threads.Atomic{Int}(0)
        project = MB.define_project("ExpandedSource")
        MB.register_item!(
            project,
            :table;
            detect=file -> endswith(file.filename, ".csv"),
            read=function (file)
                Threads.atomic_add!(read_count, 1)
                return CSV.read(file.filepath, DataFrame)
            end,
            entries=function (_file, data)
                items = MB.AbstractDataItem[]
                sizehint!(items, 100)
                for part in 1:100
                    mask = data.part .== part
                    push!(items, DataItem(
                        kind=:table,
                        collection=["expanded"],
                        parameters=Dict{Symbol,Any}(:part => part),
                        data=DataFrame(x=data.x[mask], y=data.y[mask]),
                    ))
                end
                return items
            end,
            process=item -> DataItem(
                item,
                DataFrame(
                    x=item.data.x,
                    y=item.data.y,
                    sum=item.data.x .+ item.data.y,
                ),
            ),
            stats=item -> Dict{Symbol,Any}(
                :rows => nrow(item.data),
                :sum_max => maximum(item.data.sum),
            ),
        )

        workspace = MB.open_workspace(project, test_source(project, dir))
        try
            wait_workspace_idle!(workspace)

            @test workspace.scan.state == :done
            @test read_count[] == 1
            records = workspace.index.hierarchy.all_items
            @test length(records) == 100
            @test all(record -> workspace.index.item_stats[record.id][:rows] == 10, records)

            loaded = MB.Workspace.materialize_items(workspace, records)
            @test length(loaded) == 100
            @test all(item -> nrow(item.data) == 10, loaded)
            @test read_count[] == 1
        finally
            MB.close_workspace!(workspace)
        end
    end
end

@testset "incremental reopen re-reads only changed source items" begin
    mktempdir() do dir
        for i in 1:5
            CSV.write(joinpath(dir, "f$i.csv"), DataFrame(x=Float64.(1:3), y=Float64.(4:6)))
        end

        read_count = Threads.Atomic{Int}(0)
        project = MB.define_project("Reopen")
        MB.register_item!(
            project,
            :table;
            detect=file -> endswith(file.filename, ".csv"),
            read=function (file)
                Threads.atomic_add!(read_count, 1)
                return CSV.read(file.filepath, DataFrame)
            end,
            entries=(file, data) -> [DataItem(
                kind=:table, collection=[splitext(file.filename)[1]], data=data)],
            stats=function (item)
                item.collection == ["f2"] && error("persistent f2 analysis failure")
                return Dict{Symbol,Any}(:rows => nrow(item.data))
            end,
        )

        function run_scan!()
            ws = MB.open_workspace(project, test_source(project, dir))
            wait_workspace_idle!(ws)
            @test ws.scan.state in (:done, :unchanged)
            return ws
        end

        # Cold build reads every file.
        ws1 = run_scan!()
        identity = ws1.cache.identity
        @test read_count[] == 5
        @test length(ws1.index.hierarchy.all_items) == 5
        @test ws1.cache.status.new_source_items == 5
        failed_item_id = only(
            record.id
            for record in ws1.index.hierarchy.all_items
            if record.collection == ["f2"]
        )
        @test haskey(ws1.index.analysis_errors, failed_item_id)

        # A missing fingerprint cannot prove that a source item is unchanged.
        cached_source = ws1.index.source
        unverifiable = Dict{String,Any}(
            id => nothing for id in keys(cached_source.source_item_fingerprints))
        unverifiable_source = MB.SourceScan(
            cached_source.source_id,
            cached_source.source_label,
            unverifiable,
            cached_source.hierarchy,
            cached_source.analysis_failures,
        )
        @test !MB.ItemIndex.source_unchanged(
            test_source(project, dir), unverifiable, unverifiable_source)
        _, reusable = MB.ItemIndex._incremental_reuse_plan(unverifiable_source, unverifiable)
        @test isempty(reusable)
        MB.close_workspace!(ws1)

        try
            # Reopen with nothing changed reads no files.
            Threads.atomic_xchg!(read_count, 0)
            ws2 = run_scan!()
            @test read_count[] == 0
            @test length(ws2.index.hierarchy.all_items) == 5
            @test ws2.cache.status.new_source_items == 0
            @test ws2.cache.status.stale_source_items == 0
            @test ws2.cache.status.deleted_source_items == 0
            @test haskey(ws2.index.analysis_errors, failed_item_id)
            MB.close_workspace!(ws2)

            # Changing one file re-reads only that file.
            CSV.write(joinpath(dir, "f3.csv"), DataFrame(x=Float64.(1:7), y=Float64.(8:14)))
            Threads.atomic_xchg!(read_count, 0)
            ws3 = run_scan!()
            @test read_count[] == 1
            @test ws3.cache.status.stale_source_items == 1
            @test ws3.cache.status.new_source_items == 0
            changed = only(
                r for r in ws3.index.hierarchy.all_items if r.collection == ["f3"])
            @test ws3.index.item_stats[changed.id][:rows] == 7
            @test haskey(ws3.index.analysis_errors, failed_item_id)
            MB.close_workspace!(ws3)

            # Adding a file reads only the new file.
            CSV.write(joinpath(dir, "f6.csv"), DataFrame(x=Float64.(1:2), y=Float64.(3:4)))
            Threads.atomic_xchg!(read_count, 0)
            ws4 = run_scan!()
            @test read_count[] == 1
            @test ws4.cache.status.new_source_items == 1
            @test length(ws4.index.hierarchy.all_items) == 6
            MB.close_workspace!(ws4)

            # Deleting a file reads nothing and drops its records.
            rm(joinpath(dir, "f6.csv"))
            Threads.atomic_xchg!(read_count, 0)
            ws5 = run_scan!()
            @test read_count[] == 0
            @test ws5.cache.status.deleted_source_items == 1
            @test length(ws5.index.hierarchy.all_items) == 5
            MB.close_workspace!(ws5)
        finally
            rm(dirname(identity.cache_path); force=true, recursive=true)
        end
    end
end
