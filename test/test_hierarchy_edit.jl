using DataBrowser
using Test

const HE_INDEX = DataBrowserAPI.ItemIndex

function _insert_registered!(index, item_id, names...)
    path = collect(HE_INDEX.RegisteredCollection.(names))
    key = HE_INDEX.resolve_collection_path!(index, HE_INDEX.collection_inputs(path))
    HE_INDEX.insert_item!(index, item_id, key)
    return key
end

@testset "workspace collection publication edits indexes in place" begin
    mktempdir() do dir
        project = DataBrowser.define_project("InPlacePublish_$(basename(dir))")
        workspace = DataBrowserCore.Workspace.Workspace(
            project,
            test_source(project, dir);
            cache=false,
        )
        try
            index = workspace.index
            collections = index.collections
            items = index.items
            path = HE_INDEX.collection_inputs(AbstractCollection[
                HE_INDEX.RegisteredCollection("batch"),
            ])
            for number in 1:2
                record = HE_INDEX.ItemRecord(;
                    source_item_id="source-$number",
                    id="item-$number",
                    item_label="Item $number",
                    kind=:row,
                )
                DataBrowserCore.Workspace.publish_source_item_records!(
                    workspace,
                    record.source_item_id,
                    [record],
                    [path],
                )
                if number == 1
                    collection_key = only(keys(collections.records))
                    HE_INDEX.set_collection_analysis!(
                        collections,
                        collection_key,
                        Dict(:score => 1.0),
                    )
                    DataBrowserCore.Workspace.refresh_collection_metadata_keys!(index)
                    @test index.collection_metadata_keys == [:score]
                end
            end
            @test length(items) == 2
            @test isempty(index.collection_metadata_keys)
            @test HE_INDEX.collection_item_ids(collections, only(keys(collections.records))) ==
                ["item-1", "item-2"]

            existing = items["item-1"]
            DataBrowserCore.Workspace.publish_source_item_records!(
                workspace,
                existing.source_item_id,
                [existing],
                [path],
            )
            @test HE_INDEX.collection_item_ids(collections, only(keys(collections.records))) ==
                ["item-1", "item-2"]

            DataBrowserCore.Workspace.remove_source_item_output!(workspace, "source-1")
            @test collect(keys(items)) == ["item-2"]
        finally
            DataBrowser.close_workspace!(workspace)
        end
    end
end

@testset "collection index edits" begin
    original = HE_INDEX.CollectionIndex("source")
    run_2 = _insert_registered!(original, "a1", "dev-A", "run-2")
    run_10 = _insert_registered!(original, "a2", "dev-A", "run-10")
    run_1 = _insert_registered!(original, "b1", "dev-B", "run-1")
    dev_a = original.records[run_2].parent_key
    @test [original.records[key].label for key in HE_INDEX.sorted_child_keys(original, dev_a)] ==
        ["run-2", "run-10"]

    updated = copy(original)
    HE_INDEX.remove_item!(updated, "a1", run_2)
    HE_INDEX.remove_item!(updated, "a2", run_10)
    run_3 = _insert_registered!(updated, "a3", "dev-A", "run-3")
    updated_dev_a = updated.records[run_3].parent_key

    @test !haskey(updated.records, run_2)
    @test [updated.records[key].label for key in
        HE_INDEX.sorted_child_keys(updated, updated_dev_a)] == ["run-3"]
    @test HE_INDEX.collection_item_ids(updated, run_3) == ["a3"]
    HE_INDEX.remove_item!(updated, "a3", run_3)
    @test !haskey(updated.records, updated_dev_a)
    @test haskey(original.records, dev_a)
end
