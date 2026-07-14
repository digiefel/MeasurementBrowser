using DataBrowser
using CSV
using DataFrames: DataFrame, nrow
using Test

const MB = DataBrowser

"""Build a small registry project whose CSV reader expands one item per file."""
function _profile_project()
    project = MB.define_project("ProfileProject")
    MB.register_item!(
        project,
        :table;
        detect=file -> endswith(file.filename, ".csv"),
        read=file -> DataFrame(CSV.File(file.filepath)),
        collection=(_data, metadata) ->
            ["dev", splitext(metadata[:filename])[1]],
        label=(_data, metadata) -> metadata[:filename],
        analyze=(data, _metadata) -> Dict{Symbol,Any}(:rows => nrow(data)),
    )
    return project
end

@testset "reopen applies collection parameters changed between sessions" begin
    mktempdir() do dir
        write(joinpath(dir, "ok.csv"), "x,y\n1,2\n")
        write(joinpath(dir, "metadata.txt"), "collection_path,area_um2\ndev,42\n")

        project = MB.define_project("ReopenParameters")
        MB.register_item!(
            project,
            :table;
            detect=file -> endswith(file.filename, ".csv"),
            read=file -> DataFrame(CSV.File(file.filepath)),
            collection=(_data, metadata) ->
                ["dev", splitext(metadata[:filename])[1]],
            analyze=(_data, metadata) ->
                Dict{Symbol,Any}(:area_um2 => metadata[:area_um2]),
        )

        workspace = MB.open_workspace(
            project, test_source(project, dir); background_processing=true)
        identity = workspace.cache.identity
        try
            wait_workspace_idle!(workspace)
            record = only(DataBrowserAPI.ItemIndex.all_items(workspace.index.hierarchy))
            @test workspace.index.item_metadata[record.id][:area_um2] == 42
        finally
            MB.close_workspace!(workspace)
        end

        # The metadata file changed while the workspace was closed: reopen diffs
        # source_collection_metadata and recomputes affected items.
        write(joinpath(dir, "metadata.txt"), "collection_path,area_um2\ndev,43\n")
        reopened = MB.open_workspace(
            project, test_source(project, dir); background_processing=true)
        try
            wait_workspace_idle!(reopened)
            record = only(DataBrowserAPI.ItemIndex.all_items(reopened.index.hierarchy))
            @test DataBrowserAPI.ItemIndex.effective_metadata(
                reopened.index.hierarchy, record)[:area_um2] == 43
            @test reopened.index.item_metadata[record.id][:area_um2] == 43
        finally
            MB.close_workspace!(reopened)
            rm(dirname(identity.cache_path); force=true, recursive=true)
        end
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

        workspace = MB.open_workspace(
            project, test_source(project, dir); background_processing=true)
        identity = workspace.cache.identity
        wait_workspace_idle!(workspace)
        MB.close_workspace!(workspace)
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
        @test row.analyze_seconds >= 0
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

        # A cache-hit reopen does no per-kind work, so the profile resets to empty.
        try
            reopened = MB.open_workspace(
                project, test_source(project, dir); background_processing=true)
            wait_workspace_idle!(reopened)
            MB.close_workspace!(reopened)
            @test isempty(MB.scan_profile_summary(project))
        finally
            rm(dirname(identity.cache_path); force=true, recursive=true)
        end
    end
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
            collection=(_data, metadata) ->
                ["dev", splitext(metadata[:filename])[1]],
            label=(_data, metadata) -> metadata[:filename],
            analyze=function (data, metadata)
                return Dict{Symbol,Any}(
                    :rows => nrow(data),
                    :area_um2 => metadata[:area_um2],
                )
            end,
        )

        workspace = MB.open_workspace(
            project, test_source(project, dir); background_processing=true)
        try
            wait_workspace_idle!(workspace)

            @test workspace.scan.state == :done
            @test isempty(workspace.index.analysis_errors)
            ok = only(DataBrowserAPI.ItemIndex.all_items(workspace.index.hierarchy))
            @test workspace.index.item_metadata[ok.id][:rows] == 1
            @test workspace.index.item_metadata[ok.id][:area_um2] == 42
        finally
            MB.close_workspace!(workspace)
        end

        reopened = MB.open_workspace(
            project, test_source(project, dir); background_processing=true)
        try
            wait_workspace_idle!(reopened)
            ok = only(DataBrowserAPI.ItemIndex.all_items(reopened.index.hierarchy))
            effective = DataBrowserAPI.ItemIndex.effective_metadata(reopened.index.hierarchy, ok)
            @test effective[:area_um2] == 42
            @test reopened.index.item_metadata[ok.id][:area_um2] == 42
        finally
            MB.close_workspace!(reopened)
        end
    end
end

@testset "DirectorySource metadata changes invalidate only affected items" begin
    mktempdir() do dir
        write(joinpath(dir, "a.csv"), "x\n2\n")
        write(joinpath(dir, "b.csv"), "x\n3\n")
        metadata_path = joinpath(dir, "device_info.txt")
        write(metadata_path, "collection_path,scale\ndev/a,2\ndev/b,4\n")

        reads = Dict("a.csv" => Threads.Atomic{Int}(0), "b.csv" => Threads.Atomic{Int}(0))
        entries = Dict("a.csv" => Threads.Atomic{Int}(0), "b.csv" => Threads.Atomic{Int}(0))
        processes = Dict("a" => Threads.Atomic{Int}(0), "b" => Threads.Atomic{Int}(0))
        project = MB.define_project("MetadataWatcher_$(basename(dir))")
        MB.register_item!(
            project,
            :table;
            detect=file -> endswith(file.filename, ".csv"),
            read=function (file)
                Threads.atomic_add!(reads[file.filename], 1)
                return DataFrame(CSV.File(file.filepath))
            end,
            entries=function (data, metadata)
                Threads.atomic_add!(entries[metadata[:filename]], 1)
                name = splitext(metadata[:filename])[1]
                item_metadata = name == "b" ? Dict{Symbol,Any}(:scale => 20) :
                    Dict{Symbol,Any}()
                return [(data=data, metadata=item_metadata)]
            end,
            collection=(_data, metadata) ->
                ["dev", splitext(metadata[:filename])[1]],
            process=function (data, metadata)
                name = splitext(metadata[:filename])[1]
                Threads.atomic_add!(processes[name], 1)
                scale = get(metadata, :scale, 1)
                return DataFrame(y=data.x .* scale)
            end,
            analyze=(data, _metadata) -> Dict{Symbol,Any}(:maximum => maximum(data.y)),
        )

        workspace = MB.open_workspace(
            project,
            DirectorySource(dir; metadata_file="device_info.txt");
            cache=false,
        )
        try
            wait_workspace_idle!(workspace)
            records = sort(DataBrowserAPI.ItemIndex.all_items(workspace.index.hierarchy); by=record -> record.item_label)
            DataBrowserCore.Workspace.materialize_items(workspace, records)
            wait_workspace_idle!(workspace)
            @test [counter[] for counter in values(reads)] == [1, 1]
            @test [counter[] for counter in values(entries)] == [1, 1]
            @test processes["a"][] == 1
            @test processes["b"][] == 1

            write(metadata_path, "collection_path,scale\ndev/a,3\ndev/b,5\n")
            @test Base.timedwait(
                () -> workspace.index.hierarchy.index[("dev", "a")].metadata[:scale] == 3,
                5,
            ) === :ok
            records = sort(DataBrowserAPI.ItemIndex.all_items(workspace.index.hierarchy); by=record -> record.item_label)
            loaded = DataBrowserCore.Workspace.materialize_items(workspace, records)
            wait_workspace_idle!(workspace)

            @test only(
                item_data(item).y for (record, item) in zip(records, loaded)
                if last(record.collection) == "a"
            ) == [6]
            @test only(
                item_data(item).y for (record, item) in zip(records, loaded)
                if last(record.collection) == "b"
            ) == [60]
            @test all(counter[] == 1 for counter in values(reads))
            @test all(counter[] == 1 for counter in values(entries))
            @test processes["a"][] == 2
            @test processes["b"][] == 1

            write(metadata_path, "collection_path,scale\ndev/b,5\n")
            @test Base.timedwait(
                () -> !haskey(
                    workspace.index.hierarchy.index[("dev", "a")].metadata,
                    :scale,
                ),
                5,
            ) === :ok
            records = sort(DataBrowserAPI.ItemIndex.all_items(workspace.index.hierarchy); by=record -> record.item_label)
            DataBrowserCore.Workspace.materialize_items(workspace, records)
            wait_workspace_idle!(workspace)
            @test processes["a"][] == 3
            @test processes["b"][] == 1
            @test all(counter[] == 1 for counter in values(reads))

            write(metadata_path, "collection_path,scale\ndev/a,4\ndev/b,5\n")
            @test Base.timedwait(
                () -> get(
                    workspace.index.hierarchy.index[("dev", "a")].metadata,
                    :scale,
                    nothing,
                ) == 4,
                5,
            ) === :ok
            records = sort(DataBrowserAPI.ItemIndex.all_items(workspace.index.hierarchy); by=record -> record.item_label)
            DataBrowserCore.Workspace.materialize_items(workspace, records)
            wait_workspace_idle!(workspace)
            @test processes["a"][] == 4
            @test processes["b"][] == 1
            @test all(counter[] == 1 for counter in values(reads))
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
            entries=function (data, _metadata)
                items = NamedTuple[]
                sizehint!(items, 100)
                for part in 1:100
                    mask = data.part .== part
                    push!(items, (
                        data=DataFrame(x=data.x[mask], y=data.y[mask]),
                        metadata=Dict{Symbol,Any}(:part => part),
                    ))
                end
                return items
            end,
            collection=(_data, _metadata) -> ["expanded"],
            process=(data, _metadata) -> DataFrame(
                x=data.x,
                y=data.y,
                sum=data.x .+ data.y,
            ),
            analyze=(data, _metadata) -> Dict{Symbol,Any}(
                :rows => nrow(data),
                :sum_max => maximum(data.sum),
            ),
        )

        workspace = MB.open_workspace(
            project, test_source(project, dir); background_processing=true)
        try
            wait_workspace_idle!(workspace)

            @test workspace.scan.state == :done
            @test read_count[] == 1
            records = DataBrowserAPI.ItemIndex.all_items(workspace.index.hierarchy)
            @test length(records) == 100
            @test all(record -> workspace.index.item_metadata[record.id][:rows] == 10, records)

            loaded = DataBrowserCore.Workspace.materialize_items(workspace, records)
            @test length(loaded) == 100
            @test all(item -> nrow(item_data(item)) == 10, loaded)
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
            collection=(_data, metadata) -> [splitext(metadata[:filename])[1]],
            analyze=function (data, metadata)
                metadata[:filename] == "f2.csv" &&
                    error("persistent f2 analysis failure")
                return Dict{Symbol,Any}(:rows => nrow(data))
            end,
        )

        function run_scan!()
            ws = MB.open_workspace(
                project, test_source(project, dir); background_processing=true)
            wait_workspace_idle!(ws)
            @test ws.scan.state in (:done, :unchanged)
            return ws
        end

        # Cold build reads every file.
        ws1 = run_scan!()
        identity = ws1.cache.identity
        @test read_count[] == 5
        @test length(DataBrowserAPI.ItemIndex.all_items(ws1.index.hierarchy)) == 5
        @test ws1.cache.status.new_source_items == 5
        failed_item_id = only(
            record.id
            for record in DataBrowserAPI.ItemIndex.all_items(ws1.index.hierarchy)
            if record.collection == ["f2"]
        )
        @test haskey(ws1.index.analysis_errors, failed_item_id)

        MB.close_workspace!(ws1)

        try
            # Reopen with nothing changed reads no files.
            Threads.atomic_xchg!(read_count, 0)
            ws2 = run_scan!()
            @test read_count[] == 0
            @test length(DataBrowserAPI.ItemIndex.all_items(ws2.index.hierarchy)) == 5
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
                r for r in DataBrowserAPI.ItemIndex.all_items(ws3.index.hierarchy) if r.collection == ["f3"])
            @test ws3.index.item_metadata[changed.id][:rows] == 7
            @test haskey(ws3.index.analysis_errors, failed_item_id)
            MB.close_workspace!(ws3)

            # Adding a file reads only the new file.
            CSV.write(joinpath(dir, "f6.csv"), DataFrame(x=Float64.(1:2), y=Float64.(3:4)))
            Threads.atomic_xchg!(read_count, 0)
            ws4 = run_scan!()
            @test read_count[] == 1
            @test ws4.cache.status.new_source_items == 1
            @test length(DataBrowserAPI.ItemIndex.all_items(ws4.index.hierarchy)) == 6
            MB.close_workspace!(ws4)

            # Deleting a file reads nothing and drops its records.
            rm(joinpath(dir, "f6.csv"))
            Threads.atomic_xchg!(read_count, 0)
            ws5 = run_scan!()
            @test read_count[] == 0
            @test ws5.cache.status.deleted_source_items == 1
            @test length(DataBrowserAPI.ItemIndex.all_items(ws5.index.hierarchy)) == 5
            MB.close_workspace!(ws5)
        finally
            rm(dirname(identity.cache_path); force=true, recursive=true)
        end
    end
end
