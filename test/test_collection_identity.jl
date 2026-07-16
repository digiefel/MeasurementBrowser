using DataBrowser
using Test

struct DefaultIdentityCollection <: AbstractCollection
    value::Int
end

struct ExplicitIdentityCollection <: AbstractCollection
    key::Int
    shown::String
end

DataBrowser.id(collection::ExplicitIdentityCollection) = collection.key
DataBrowser.label(collection::ExplicitIdentityCollection) = collection.shown

struct OtherIdentityCollection <: AbstractCollection
    value::Int
end

@testset "collection identity contract" begin
    default = DefaultIdentityCollection(2)
    @test DataBrowser.id(default) == default
    @test DataBrowser.label(default) == string(default)
    @test DataBrowser.metadata(default) == Dict()

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
