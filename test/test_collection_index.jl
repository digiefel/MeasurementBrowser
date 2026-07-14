using DataBrowser
using Test

const COLLECTION_RECORD_INDEX = DataBrowserAPI.ItemIndex

struct IndexedCollection <: AbstractCollection
    key::Int
    shown::String
    area::Float64
end

DataBrowser.id(collection::IndexedCollection) = collection.key
DataBrowser.label(collection::IndexedCollection) = collection.shown
DataBrowser.metadata(collection::IndexedCollection) = Dict(:area => collection.area)

function _index_collection!(index, item_id, path)
    inputs = COLLECTION_RECORD_INDEX.collection_inputs(AbstractCollection[path...])
    key = COLLECTION_RECORD_INDEX.resolve_collection_path!(index, inputs)
    COLLECTION_RECORD_INDEX.insert_item!(index, item_id, key)
    return key
end

@testset "collection record index" begin
    index = COLLECTION_RECORD_INDEX.CollectionIndex("source")
    first = _index_collection!(index, "a", [
        IndexedCollection(1, "parent", 10.0),
        IndexedCollection(2, "same label", 20.0),
    ])
    second = _index_collection!(index, "b", [
        IndexedCollection(1, "renamed parent", 11.0),
        IndexedCollection(3, "same label", 30.0),
    ])

    @test first != second
    @test COLLECTION_RECORD_INDEX.collection_location(index, first) ==
        ["renamed parent", "same label"]
    @test COLLECTION_RECORD_INDEX.collection_item_ids(index, first) == ["a"]
    @test COLLECTION_RECORD_INDEX.collection_item_ids(index, second) == ["b"]
    @test COLLECTION_RECORD_INDEX.collection_metadata(index, first) ==
        Dict(:area => 20.0)

    copied = copy(index)
    COLLECTION_RECORD_INDEX.remove_item!(copied, "a", first)
    @test !haskey(copied.records, first)
    @test haskey(index.records, first)

    COLLECTION_RECORD_INDEX.set_collection_analysis!(index, second, Dict(:count => 1))
    @test index.records[second].analysis == Dict(:count => 1)
    COLLECTION_RECORD_INDEX.clear_collection_analysis!(index, second)
    @test isempty(index.records[second].analysis)
end
