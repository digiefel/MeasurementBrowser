using DataBrowserAPI: AbstractCollection
using DataBrowserAPI
using Test

struct DefaultIdentityCollection <: AbstractCollection
    value::Int
end

DataBrowserAPI.id(collection::DefaultIdentityCollection) = string(collection.value)

struct ExplicitIdentityCollection <: AbstractCollection
    key::Int
    shown::String
end

DataBrowserAPI.id(collection::ExplicitIdentityCollection) = string(collection.key)
DataBrowserAPI.label(collection::ExplicitIdentityCollection) = collection.shown

struct OtherIdentityCollection <: AbstractCollection
    value::Int
end

DataBrowserAPI.id(collection::OtherIdentityCollection) = string(collection.value)

@testset "collection identity contract" begin
    original = only(DataBrowserAPI.collection_id_path([
        ExplicitIdentityCollection(42, "before"),
    ]))
    @test original == only(DataBrowserAPI.collection_id_path([
        ExplicitIdentityCollection(42, "after"),
    ]))
    @test original != only(DataBrowserAPI.collection_id_path([
        ExplicitIdentityCollection(43, "before"),
    ]))

    @test only(DataBrowserAPI.collection_id_path([DefaultIdentityCollection(1)])) !=
        only(DataBrowserAPI.collection_id_path([OtherIdentityCollection(1)]))

    parent_a = DataBrowserAPI.collection_id_path([
        DefaultIdentityCollection(1),
        DefaultIdentityCollection(3),
    ])
    parent_b = DataBrowserAPI.collection_id_path([
        DefaultIdentityCollection(2),
        DefaultIdentityCollection(3),
    ])
    @test last(parent_a) != last(parent_b)
end
