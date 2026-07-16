using DataBrowser
using CancellationTokens: CancellationTokenSource, get_token
using DataBrowserAPI: AbstractDataSource, source_item_noun
using Test

const MBD = DataBrowser

struct TestNounSource <: AbstractDataSource end

"""Open a memory-only workspace on `source` and wait for the build to settle."""
function _open_settled(project, source)
    workspace = MBD.open_workspace(project, source; cache=false)
    wait_workspace_idle!(workspace)
    return workspace
end

@testset "DirectorySource discovery and metadata" begin
    mktempdir() do dir
        write_test_source(joinpath(dir, "root.csv"))
        mkdir(joinpath(dir, "nested"))
        write_test_source(joinpath(dir, "nested", "nested.csv"), 10)
        write(joinpath(dir, ".hidden.csv"), "")
        write(
            joinpath(dir, "metadata.txt"),
            "collection_path,wafer,area_um2,zero,one,flag\n" *
            "test,A,12.5,0,1,yes\n" *
            "test/nested,B,99,0,1,no\n",
        )
        write(joinpath(dir, "tags.txt"), "ignored")
        write(joinpath(dir, "databrowser.toml"), "ignored")

        project = MBD.define_project("Parameter test")
        MBD.register_item!(
            project,
            :table;
            detect=file -> endswith(lowercase(file.filename), ".csv"),
            read=file -> read(file.filepath, String),
            entries=function (data, source_metadata)
                name = splitext(source_metadata[:filename])[1]
                item_metadata = name == "root" ?
                    Dict{Symbol,Any}(:area_um2 => 7.0, :local_only => true) :
                    Dict{Symbol,Any}()
                return [(data=data, metadata=item_metadata)]
            end,
            collection=(_data, metadata) ->
                ["test", splitext(metadata[:filename])[1]],
            label=(_data, metadata) ->
                "Table $(splitext(metadata[:filename])[1])",
        )

        workspace = _open_settled(project, MBD.DirectorySource(dir))
        try
            # Hidden files and sidecars are never source items.
            @test workspace.cache.status.total_source_items == 2
            collections = workspace.index.collections
            root_collection = _registered_collection_record(collections, "test", "root")
            nested_collection = _registered_collection_record(collections, "test", "nested")
            root_metadata = DataBrowserAPI.ItemIndex.collection_metadata(
                collections, root_collection.key)
            nested_metadata = DataBrowserAPI.ItemIndex.collection_metadata(
                collections, nested_collection.key)
            @test root_metadata[:wafer] == "A"
            @test root_metadata[:area_um2] == 12.5
            @test root_metadata[:zero] === 0
            @test root_metadata[:one] === 1
            @test root_metadata[:flag] === true
            @test nested_metadata[:wafer] == "B"
            @test nested_metadata[:area_um2] == 99
            @test nested_metadata[:zero] === 0
            @test nested_metadata[:one] === 1
            @test nested_metadata[:flag] === false

            metadata_by_label = Dict(
                record.item_label =>
                    DataBrowserAPI.ItemIndex.effective_metadata(collections, record)
                for record in values(workspace.index.items)
            )
            @test metadata_by_label["Table root"][:wafer] == "A"
            @test metadata_by_label["Table root"][:area_um2] == 7.0
            @test metadata_by_label["Table root"][:local_only]
            @test metadata_by_label["Table nested"][:wafer] == "B"
            @test metadata_by_label["Table nested"][:area_um2] == 99
        finally
            MBD.close_workspace!(workspace)
        end

        shallow = _open_settled(TEST_PROJECT, MBD.DirectorySource(dir; recursive=false))
        try
            @test shallow.cache.status.total_source_items == 1
            @test only(collect(values(shallow.index.items))).item_label == "Table root"
        finally
            MBD.close_workspace!(shallow)
        end

        without_metadata = _open_settled(
            TEST_PROJECT, MBD.DirectorySource(dir; metadata_file=nothing))
        try
            collections = without_metadata.index.collections
            @test all(isempty(collection_record.own_metadata)
                for collection_record in values(collections.records))
            @test all(
                record.metadata == Dict{Symbol,Any}(:filename => basename(record.source_item_id))
                for record in values(without_metadata.index.items)
            )
        finally
            MBD.close_workspace!(without_metadata)
        end

        custom_path = joinpath(dir, "device_info.txt")
        write(custom_path, "collection_path,wafer\ntest,C\n")
        custom = MBD.DirectorySource(dir; metadata_file="device_info.txt")
        MBD.open_source(custom)
        try
            @test source_item_noun(custom) == "source files"
            @test source_item_noun(TestNounSource()) == "source items"
            @test custom.has_metadata
            cancel_token = get_token(CancellationTokenSource())
            @test all(
                file -> file.filepath != custom_path,
                MBD.source_items(custom; cancel_token),
            )
        finally
            MBD.close_source!(custom)
        end
    end
end

@testset "malformed live metadata is a recoverable source error" begin
    mktempdir() do dir
        write_test_source(joinpath(dir, "a.csv"))
        metadata_path = joinpath(dir, "metadata.txt")
        write(metadata_path, "collection_path,wafer\ntest,A\n")

        workspace = MBD.open_workspace(
            TEST_PROJECT, MBD.DirectorySource(dir); cache=false)
        try
            wait_workspace_idle!(workspace)
            @test isempty(workspace.source_error)

            write(metadata_path, "no_parameter_columns\n")
            @test Base.timedwait(() -> !isempty(workspace.source_error), 5) === :ok
            DataBrowserCore.Workspace.refresh_status!(workspace)
            status = DataBrowserCore.Workspace.workspace_status(workspace)
            @test status.level === :error
            @test status.label == "Source Error"

            # A corrected file clears the failure without restarting anything.
            write(metadata_path, "collection_path,wafer\ntest,B\n")
            @test Base.timedwait(() -> isempty(workspace.source_error), 5) === :ok
            @test Base.timedwait(
                () -> get(
                    _registered_collection_record(
                        workspace.index.collections, "test").own_metadata,
                    :wafer,
                    nothing,
                ) == "B",
                5,
            ) === :ok
        finally
            MBD.close_workspace!(workspace)
        end
    end
end

@testset "directory watcher publishes file changes" begin
    mktempdir() do dir
        write_test_source(joinpath(dir, "a.csv"))
        workspace = MBD.open_workspace(
            TEST_PROJECT, MBD.DirectorySource(dir); cache=false)
        try
            wait_workspace_idle!(workspace)
            @test length(workspace.index.items) == 1

            new_path = joinpath(dir, "b.csv")
            write_test_source(new_path, 10)
            @test Base.timedwait(() -> length(workspace.index.items) == 2, 5) === :ok

            rm(new_path)
            @test Base.timedwait(() -> length(workspace.index.items) == 1, 5) === :ok
        finally
            MBD.close_workspace!(workspace)
        end
    end
end
