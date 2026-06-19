using MeasurementBrowser
using MeasurementBrowser: ItemRecord
using Test

@testset "threaded scan publishes each source item once" begin
    mktempdir() do dir
        for index in 1:20
            write(joinpath(dir, "item_$index.txt"), "value=$index\n")
        end

        seen = Dict{String,Int}()
        seen_lock = ReentrantLock()
        project = MeasurementBrowser.define_project("Threaded scan test")
        MeasurementBrowser.register_item!(
            project,
            :text;
            detect=file -> endswith(file.filename, ".txt"),
            read=function (file)
                lock(seen_lock) do
                    seen[file.filepath] = get(seen, file.filepath, 0) + 1
                end
                sleep(0.001)
                return read(file.filepath, String)
            end,
            entries=(file, data) -> [MeasurementBrowser.DataItem(
                kind=:text,
                collection=["text"],
                label=file.filename,
                id=file.filepath * "#text",
                data=data,
            )],
        )

        streamed = String[]
        scan = MeasurementBrowser.scan_source(
            project,
            MeasurementBrowser.DirectorySource(dir; metadata_file=nothing);
            on_items=items -> append!(streamed, item.source_item_id for item in items),
        )

        @test length(scan.hierarchy.all_items) == 20
        @test length(streamed) == 20
        @test length(seen) == 20
        @test all(count -> count == 1, values(seen))
        @test Set(streamed) == Set(keys(seen))
    end
end

@testset "scan_source progress and cancel" begin
    mktempdir() do dir
        write_test_source(joinpath(dir, "first.csv"))
        write_test_source(joinpath(dir, "second.csv"), 10)

        source = scan_test_source(TEST_PROJECT, dir)
        @test source isa MeasurementBrowser.SourceScan
        @test source.hierarchy isa MeasurementBrowser.Hierarchy
        @test source.hierarchy.skipped_count == 0

        events = NamedTuple[]
        MeasurementBrowser.scan_source(
            TEST_PROJECT,
            test_source(TEST_PROJECT, dir);
            on_progress=progress -> push!(events, progress),
            count_first=true,
        )
        @test any(event -> event.phase == :counting, events)
        @test any(event -> event.phase == :discovering, events)
        @test any(event -> event.phase == :scanning, events)
        scanning = filter(event -> event.phase == :scanning, events)
        @test maximum(event -> event.processed_source_items, scanning) == 2
        @test all(event -> event.total_source_items == 2, scanning)

        write(joinpath(dir, ".hidden.csv"), "")
        source = scan_test_source(TEST_PROJECT, dir)
        @test length(source.source_item_fingerprints) == 2

        streamed = Vector{ItemRecord}[]
        source = MeasurementBrowser.scan_source(
            TEST_PROJECT,
            test_source(TEST_PROJECT, dir);
            on_items=items -> push!(streamed, copy(items)),
        )
        @test sum(length, streamed) == length(source.hierarchy.all_items)
        @test Set(
            item.id
            for batch in streamed
            for item in batch
        ) == Set(
            item.id
            for item in source.hierarchy.all_items
        )

        workspace = MeasurementBrowser.Workspace.Workspace(
            TEST_PROJECT,
            test_source(TEST_PROJECT, dir),
        )
        old_index = workspace.index.items
        @test MeasurementBrowser.Workspace.append_items!(
            workspace,
            source.hierarchy.all_items[1:1],
        )
        @test old_index !== workspace.index.items
        @test isempty(old_index)
        @test length(workspace.index.items) == 1

        write_test_source(joinpath(dir, "broken.csv"))
        source = scan_test_source(TEST_PROJECT, dir)
        @test length(source.source_item_fingerprints) == 3
        @test length(source.analysis_failures) == 1
        @test only(source.analysis_failures).source_item_id == joinpath(dir, "broken.csv")

        fired = Ref(false)
        @test_throws MeasurementBrowser.JobCancelled MeasurementBrowser.with_cancel(
            () -> begin
                fired[] && return true
                fired[] = true
                return false
            end,
        ) do
            MeasurementBrowser.scan_source(
                TEST_PROJECT,
                test_source(TEST_PROJECT, dir);
                count_first=true,
            )
        end
    end
end

@testset "DirectorySource discovery and metadata" begin
    mktempdir() do dir
        write_test_source(joinpath(dir, "root.csv"))
        mkdir(joinpath(dir, "nested"))
        write_test_source(joinpath(dir, "nested", "nested.csv"), 10)
        write(joinpath(dir, "metadata.txt"), "collection_path,area_um2\ntest/root,12.5\ntest/nested,99\n")
        write(joinpath(dir, "tags.txt"), "ignored")
        write(joinpath(dir, "measurementbrowser.toml"), "ignored")

        source = MeasurementBrowser.DirectorySource(dir)
        scan = MeasurementBrowser.scan_source(TEST_PROJECT, source)
        @test source isa MeasurementBrowser.DirectorySource
        @test length(scan.source_item_fingerprints) == 2
        @test scan.hierarchy.has_collection_metadata
        metadata_by_label = Dict(
            record.item_label => record.collection_metadata
            for record in scan.hierarchy.all_items
        )
        @test metadata_by_label["Table root"][:area_um2] == 12.5
        @test metadata_by_label["Table nested"][:area_um2] == 99

        shallow = MeasurementBrowser.scan_source(
            TEST_PROJECT,
            MeasurementBrowser.DirectorySource(dir; recursive=false),
        )
        @test length(shallow.source_item_fingerprints) == 1
        @test only(shallow.hierarchy.all_items).item_label == "Table root"

        without_metadata = MeasurementBrowser.scan_source(
            TEST_PROJECT,
            MeasurementBrowser.DirectorySource(dir; metadata_file=nothing),
        )
        @test !without_metadata.hierarchy.has_collection_metadata
        @test all(isempty(record.collection_metadata) for record in without_metadata.hierarchy.all_items)
    end
end
