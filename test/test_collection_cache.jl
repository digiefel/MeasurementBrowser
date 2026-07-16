using DataBrowser
using DataBrowserAPI
using DataBrowserCache
using DataBrowserSources
using DBInterface
using DuckDB
using Test

struct CacheCollectionLevel <: AbstractCollection
    name::String
    value::Int
end

DataBrowserAPI.label(collection::CacheCollectionLevel) = collection.name
DataBrowserAPI.metadata(collection::CacheCollectionLevel) = Dict(:value => collection.value)

@testset "cache restores package-owned collection records" begin
    mktempdir() do dir
        filepath = joinpath(dir, "item.dat")
        write(filepath, "data")
        source_item = DataBrowserSources.index_source_file(filepath)
        cache_identity = ProjectCacheIdentity(
            "CollectionRoundTrip", dir, basename(dir), joinpath(dir, "cache.duckdb"))
        collections = DataBrowserAPI.ItemIndex.CollectionIndex(dir)
        path = AbstractCollection[
            CacheCollectionLevel("parent", 1),
            CacheCollectionLevel("leaf", 2),
        ]
        leaf_key = DataBrowserAPI.ItemIndex.resolve_collection_path!(
            collections, DataBrowserAPI.ItemIndex.collection_inputs(path))
        record = DataBrowserAPI.ItemIndex.ItemRecord(;
            id="item-1",
            source_item_id=filepath,
            source_item_path=filepath,
            label="item",
            kind=:test,
            collection_key=leaf_key,
        )
        DataBrowserAPI.ItemIndex.insert_item!(collections, record.id, leaf_key)

        cache = open_cache_db(cache_identity; rebuild=true)
        try
            write_meta_header!(cache)
            item = DataBrowserAPI.ItemIndex.RegisteredDataItem(record, nothing)
            DataBrowserCache.store_interpreted_records!(cache, source_item, [record], [item])
            store_collection_index!(cache, collections, [record])
            store_collection_metadata!(cache, leaf_key, Dict(:mean => 2.5))
            store_collection_process_result!(cache, leaf_key)
        finally
            close_cache_db!(cache)
        end

        reopened = open_cache_db(cache_identity)
        try
            index = load_cache_index(reopened)
            restored = only(index.source.items)
            @test restored.collection_key == leaf_key
            restored_collection = index.source.collections.records[leaf_key]
            @test restored_collection.id == collections.records[leaf_key].id
            @test restored_collection.label == "leaf"
            @test restored_collection.own_metadata == Dict(:value => 2)
            @test restored_collection.analysis[:mean] == 2.5
            @test haskey(index.result_states,
                CacheResultKey(COLLECTION_PROCESS_RESULT, leaf_key))
        finally
            close_cache_db!(reopened)
        end

        db = DBInterface.connect(DuckDB.DB, cache_identity.cache_path)
        connection = DBInterface.connect(db)
        try
            columns = Set(String(row.name) for row in
                DBInterface.execute(connection, "PRAGMA table_info('collections')"))
            @test "id" in columns
            @test !("value_hex" in columns)
            @test !("id_preimage_hex" in columns)
        finally
            DBInterface.close!(connection)
            DBInterface.close!(db)
        end
    end
end
