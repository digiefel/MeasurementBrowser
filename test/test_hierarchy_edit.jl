using DataBrowser
using Test

const HE_INDEX = DataBrowserAPI.ItemIndex

function _edit_record(;
    id::String,
    source_item_id::String,
    collection::Vector{String},
)::HE_INDEX.ItemRecord
    return HE_INDEX.ItemRecord(;
        id,
        source_item_id,
        item_label=id,
        kind=:iv,
        collection,
    )
end

@testset "copy-on-write hierarchy edits" begin
    mktempdir() do dir
        source = DataBrowser.DirectorySource(dir)
        records = [
            _edit_record(id="a1", source_item_id="file-a", collection=["dev-A", "run-2"]),
            _edit_record(id="a2", source_item_id="file-a", collection=["dev-A", "run-10"]),
            _edit_record(id="b1", source_item_id="file-b", collection=["dev-B", "run-1"]),
        ]
        original = HE_INDEX.Hierarchy(records, source)
        merge!(original.index[("dev-B",)].analysis, Dict(:count => 1))
        # Natural sort: run-2 before run-10.
        @test [child.name for child in original.index[("dev-A",)].children] ==
            ["run-2", "run-10"]

        # Replace file-a's items: one removed collection, one new collection.
        edit = HE_INDEX.edit_hierarchy(original)
        HE_INDEX.remove_records!(
            edit, [record for record in records if record.source_item_id == "file-a"])
        replacement = _edit_record(
            id="a3", source_item_id="file-a", collection=["dev-A", "run-3"])
        HE_INDEX.insert_record!(edit, replacement)
        HE_INDEX.clear_node_analysis!(edit, ("dev-A",))
        updated = HE_INDEX.finish_edit!(edit, source)

        # The original tree is untouched: full contract for concurrent readers.
        @test length(DataBrowserAPI.ItemIndex.all_items(original)) == 3
        @test haskey(original.index, ("dev-A", "run-2"))
        @test [child.name for child in original.index[("dev-A",)].children] ==
            ["run-2", "run-10"]

        # The updated tree reflects the edit: empty nodes pruned, new node in sorted position.
        @test length(DataBrowserAPI.ItemIndex.all_items(updated)) == 2
        @test !haskey(updated.index, ("dev-A", "run-2"))
        @test !haskey(updated.index, ("dev-A", "run-10"))
        @test [child.name for child in updated.index[("dev-A",)].children] == ["run-3"]
        @test updated.index[("dev-A", "run-3")].items[1].id == "a3"
        # Untouched subtree nodes are shared, and their stats survive the edit.
        @test updated.index[("dev-B", "run-1")] === original.index[("dev-B", "run-1")]
        @test updated.index[("dev-B",)].analysis == Dict(:count => 1)

        # Removing the last record of a device prunes the whole branch up to the root.
        second = HE_INDEX.edit_hierarchy(updated)
        HE_INDEX.remove_records!(second, [replacement])
        pruned = HE_INDEX.finish_edit!(second, source)
        @test !haskey(pruned.index, ("dev-A",))
        @test [child.name for child in pruned.root.children] == ["dev-B"]
        @test haskey(updated.index, ("dev-A", "run-3"))
    end
end
