using MeasurementBrowser
using Test

const MBD = MeasurementBrowser

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
            "collection_path,wafer,area_um2\ntest,A,12.5\ntest/nested,B,99\n",
        )
        write(joinpath(dir, "tags.txt"), "ignored")
        write(joinpath(dir, "measurementbrowser.toml"), "ignored")

        project = MBD.define_project("Parameter test")
        MBD.register_item!(
            project,
            :table;
            detect=file -> endswith(lowercase(file.filename), ".csv"),
            read=file -> read(file.filepath, String),
            entries=function (file, data)
                name = splitext(file.filename)[1]
                metadata = name == "root" ?
                    Dict{Symbol,Any}(:area_um2 => 7.0, :local_only => true) :
                    Dict{Symbol,Any}()
                return [MBD.DataItem(
                    kind=:table,
                    collection=["test", name],
                    label="Table $name",
                    metadata=metadata,
                    data=data,
                )]
            end,
        )

        workspace = _open_settled(project, MBD.DirectorySource(dir))
        try
            # Hidden files and sidecars are never source items.
            @test workspace.cache.status.total_source_items == 2
            hierarchy = workspace.index.hierarchy
            @test hierarchy.has_collection_metadata
            root_node = hierarchy.index[("test", "root")]
            nested_node = hierarchy.index[("test", "nested")]
            @test root_node.metadata[:wafer] == "A"
            @test root_node.metadata[:area_um2] == 12.5
            @test nested_node.metadata[:wafer] == "B"
            @test nested_node.metadata[:area_um2] == 99

            metadata_by_label = Dict(
                record.item_label =>
                    MBD.ItemIndex.effective_metadata(hierarchy, record)
                for record in hierarchy.all_items
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
            @test only(shallow.index.hierarchy.all_items).item_label == "Table root"
        finally
            MBD.close_workspace!(shallow)
        end

        without_metadata = _open_settled(
            TEST_PROJECT, MBD.DirectorySource(dir; metadata_file=nothing))
        try
            hierarchy = without_metadata.index.hierarchy
            @test !hierarchy.has_collection_metadata
            @test all(isempty(node.metadata) for node in values(hierarchy.index))
            @test all(isempty(record.metadata) for record in hierarchy.all_items)
        finally
            MBD.close_workspace!(without_metadata)
        end

        custom_path = joinpath(dir, "device_info.txt")
        write(custom_path, "collection_path,wafer\ntest,C\n")
        custom = MBD.DirectorySource(dir; metadata_file="device_info.txt")
        MBD.open_source(custom)
        try
            @test MBD.Projects.has_collection_metadata(custom)
            @test all(
                file -> file.filepath != custom_path,
                MBD.source_items(custom),
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
            MBD.Workspace.refresh_status!(workspace)
            status = MBD.Workspace.workspace_status(workspace)
            @test status.level === :error
            @test status.label == "Source Error"

            # A corrected file clears the failure without restarting anything.
            write(metadata_path, "collection_path,wafer\ntest,B\n")
            @test Base.timedwait(() -> isempty(workspace.source_error), 5) === :ok
            @test Base.timedwait(
                () -> get(
                    workspace.index.hierarchy.index[("test",)].metadata,
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
