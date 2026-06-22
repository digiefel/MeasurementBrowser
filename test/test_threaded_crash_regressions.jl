using MeasurementBrowser
using Test

"""Poll a workspace until its scan settles (or a timeout), so tests can assert the steady state."""
function poll_to_done(workspace; timeout=20)
    deadline = time() + timeout
    while time() < deadline
        MeasurementBrowser.Workspace.poll_workspace!(workspace)
        workspace.scan.state in (:done, :unchanged) && return nothing
        sleep(0.01)
    end
    return nothing
end

@testset "threaded workspace smoke test" begin
    @test Threads.nthreads() >= 4
    mktempdir() do root
        write(
            joinpath(root, "metadata.txt"),
            "collection_path,wafer,owner\nroot,W,collection\nroot/g0,W0,group\nroot/g1,W1,group\nroot/g2,W2,group\nroot/g3,W3,group\n",
        )
        for index in 1:48
            write(joinpath(root, "item_$index.txt"), "value=$index\n")
        end

        project = MeasurementBrowser.define_project("Threaded workspace smoke")
        MeasurementBrowser.register_item!(
            project,
            :text;
            detect=file -> endswith(file.filename, ".txt"),
            read=function (file)
                sleep(0.002)
                return read(file.filepath, String)
            end,
            entries=function (file, data)
                index = parse(Int, replace(splitext(file.filename)[1], "item_" => ""))
                group = "g$(index % 4)"
                return [MeasurementBrowser.DataItem(
                    kind=:text,
                    collection=["root", group],
                    label=file.filename,
                    parameters=Dict{Symbol,Any}(:index => index, :owner => "item"),
                    data=data,
                    id=file.filepath * "#text",
                )]
            end,
            stats=item -> Dict{Symbol,Any}(
                :bytes => length(item.data),
                :wafer => item.parameters[:wafer],
                :owner => item.parameters[:owner],
            ),
        )

        workspace = MeasurementBrowser.open_workspace(
            project,
            MeasurementBrowser.DirectorySource(root),
        )
        try
            deadline = time() + 20
            while time() < deadline
                MeasurementBrowser.Workspace.poll_workspace!(workspace)
                if workspace.scan.state in (:done, :unchanged) &&
                   workspace.analysis.state == :done
                    break
                end
                sleep(0.01)
            end

            @test workspace.scan.state in (:done, :unchanged)
            @test workspace.analysis.state == :done
            @test length(workspace.index.items) == 48
            records = collect(values(workspace.index.items))
            @test all(startswith(String(record.parameters[:wafer]), "W") for record in records)
            @test all(record.parameters[:owner] == "item" for record in records)
            @test all(record.stats[:bytes] > 0 for record in records)
            @test all(record.stats[:owner] == "item" for record in records)
        finally
            MeasurementBrowser.close_workspace!(workspace)
        end
    end
end

@testset "parallel expanded-item processing preserves order and failures" begin
    @test Threads.nthreads() >= 4
    mktempdir() do root
        write(joinpath(root, "expanded.txt"), "expanded")
        project = MeasurementBrowser.define_project("Parallel expansion")
        MeasurementBrowser.register_item!(
            project,
            :expanded;
            detect=file -> file.filename == "expanded.txt",
            read=file -> read(file.filepath, String),
            entries=(file, _data) -> MeasurementBrowser.AbstractDataItem[
                MeasurementBrowser.DataItem(
                    kind=:expanded,
                    collection=["expanded"],
                    parameters=Dict{Symbol,Any}(:index => index),
                    data=index,
                    id="$(file.filepath)#$index",
                )
                for index in 1:80
            ],
            process=function (item)
                item.parameters[:index] == 17 && error("process failure 17")
                return MeasurementBrowser.DataItem(item, item.data * 2)
            end,
            stats=function (item)
                item.parameters[:index] == 23 && error("stats failure 23")
                return Dict{Symbol,Any}(:value => item.data)
            end,
        )

        scan = scan_test_source(project, root)
        records = scan.hierarchy.all_items
        @test [record.parameters[:index] for record in records] == collect(1:80)
        @test length(scan.analysis_failures) == 2
        @test any(occursin("process failure 17", failure.message)
            for failure in scan.analysis_failures)
        @test any(occursin("stats failure 23", failure.message)
            for failure in scan.analysis_failures)
        @test !haskey(records[17].stats, :value)
        @test !haskey(records[23].stats, :value)
        @test records[80].stats[:value] == 160
    end
end

@testset "late progressive items never mutate a scanned hierarchy" begin
    # Regression: a rescan streams :items while the previous scan's analysis task still iterates the
    # live hierarchy. Mutating that hierarchy from poll_workspace! raced the analysis task and
    # corrupted the heap. Once index.source is set, progressive items must be ignored; the next
    # :source replaces the index wholesale.
    mktempdir() do root
        write(joinpath(root, "a.txt"), "value=1\n")
        project = MeasurementBrowser.define_project("Late items regression")
        MeasurementBrowser.register_item!(
            project,
            :text;
            detect=file -> endswith(file.filename, ".txt"),
            read=file -> read(file.filepath, String),
            entries=(file, data) -> [MeasurementBrowser.DataItem(
                kind=:text,
                collection=["root"],
                label=file.filename,
                data=data,
                id=file.filepath * "#text",
            )],
        )

        workspace = MeasurementBrowser.open_workspace(
            project, MeasurementBrowser.DirectorySource(root))
        try
            deadline = time() + 20
            while time() < deadline
                MeasurementBrowser.Workspace.poll_workspace!(workspace)
                workspace.scan.state in (:done, :unchanged) && break
                sleep(0.01)
            end
            @test workspace.scan.state in (:done, :unchanged)
            @test workspace.index.source !== nothing
            hierarchy = workspace.index.hierarchy
            item_count = length(workspace.index.items)

            # Inject a late :items event the way a rescan would, with index.source already set.
            channel = Channel{NamedTuple}(Inf)
            workspace.scan.id += 1
            workspace.scan.events = channel
            workspace.scan.state = :scanning
            late = MeasurementBrowser.ItemIndex.ItemRecord(;
                id="late#text",
                source_item_id="late",
                item_label="late",
                kind=:text,
                collection=["root"],
            )
            put!(channel, (kind=:items, job_id=workspace.scan.id, items=[late]))
            MeasurementBrowser.Workspace.poll_workspace!(workspace)

            @test !haskey(workspace.index.items, "late#text")
            @test length(workspace.index.items) == item_count
            @test workspace.index.hierarchy === hierarchy
        finally
            MeasurementBrowser.close_workspace!(workspace)
        end
    end
end

@testset "rescan detaches the previous scan for live, race-free updates" begin
    mktempdir() do root
        write(joinpath(root, "a.txt"), "value=1\n")
        project = MeasurementBrowser.define_project("Rescan detach")
        MeasurementBrowser.register_item!(
            project,
            :text;
            detect=file -> endswith(file.filename, ".txt"),
            read=file -> read(file.filepath, String),
            entries=(file, data) -> [MeasurementBrowser.DataItem(
                kind=:text, collection=["root"], label=file.filename,
                data=data, id=file.filepath * "#text")],
        )
        workspace = MeasurementBrowser.open_workspace(
            project, MeasurementBrowser.DirectorySource(root))
        try
            poll_to_done(workspace)
            @test workspace.index.source !== nothing
            first_count = length(workspace.index.items)

            # A rescan must immediately detach the previous scan (so a still-running analysis keeps
            # its own hierarchy) and clear stale errors, then recover a complete index live.
            workspace.index.analysis_errors["stale"] = "old error"
            MeasurementBrowser.Workspace.scan_source!(workspace)
            @test workspace.index.source === nothing
            @test isempty(workspace.index.analysis_errors)

            poll_to_done(workspace)
            @test workspace.index.source !== nothing
            @test length(workspace.index.items) == first_count
        finally
            MeasurementBrowser.close_workspace!(workspace)
        end
    end
end
