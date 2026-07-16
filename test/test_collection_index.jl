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
        IndexedCollection(1, "parent", 10.0),
        IndexedCollection(3, "same label", 30.0),
    ])

    @test first != second
    @test COLLECTION_RECORD_INDEX.collection_location(index, first) ==
        ["parent", "same label"]
    @test COLLECTION_RECORD_INDEX.collection_item_ids(index, first) == ["a"]
    @test COLLECTION_RECORD_INDEX.collection_item_ids(index, second) == ["b"]
    @test COLLECTION_RECORD_INDEX.collection_metadata(index, first) ==
        Dict(:area => 20.0)

    copied = copy(index)
    COLLECTION_RECORD_INDEX.remove_item!(copied, "a", first)
    @test !haskey(copied.records, first)
    @test haskey(index.records, first)

    registered = COLLECTION_RECORD_INDEX.CollectionIndex("registered")
    registered_inputs = COLLECTION_RECORD_INDEX.collection_inputs(AbstractCollection[
        COLLECTION_RECORD_INDEX.RegisteredCollection("wafer"),
        COLLECTION_RECORD_INDEX.RegisteredCollection("device"),
    ])
    registered_key = COLLECTION_RECORD_INDEX.resolve_collection_path!(
        registered, registered_inputs)
    @test COLLECTION_RECORD_INDEX.registration_names(registered, registered_key) ==
        ["wafer", "device"]

    conflicting = COLLECTION_RECORD_INDEX.CollectionIndex("conflicting")
    original_inputs = COLLECTION_RECORD_INDEX.collection_inputs(AbstractCollection[
        IndexedCollection(1, "parent", 10.0),
    ])
    conflicting_inputs = COLLECTION_RECORD_INDEX.collection_inputs(AbstractCollection[
        IndexedCollection(1, "renamed parent", 11.0),
    ])
    key = COLLECTION_RECORD_INDEX.resolve_collection_path!(conflicting, original_inputs)
    @test_throws ArgumentError COLLECTION_RECORD_INDEX.resolve_collection_path!(
        conflicting, conflicting_inputs)
    @test conflicting.records[key].label == "parent"
    @test conflicting.records[key].own_metadata == Dict(:area => 10.0)

    # Authoritative source metadata refreshes are allowed to replace an existing projection.
    @test COLLECTION_RECORD_INDEX.resolve_collection_path!(
        conflicting, conflicting_inputs; update_existing=true) == key
    @test conflicting.records[key].label == "renamed parent"
    @test conflicting.records[key].own_metadata == Dict(:area => 11.0)

    batched = COLLECTION_RECORD_INDEX.CollectionIndex("batched")
    @test_throws ArgumentError COLLECTION_RECORD_INDEX.resolve_collection_paths!(
        batched,
        [original_inputs, conflicting_inputs],
    )
    @test isempty(batched.records)
end
