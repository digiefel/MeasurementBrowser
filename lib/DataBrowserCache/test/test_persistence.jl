using DataBrowserAPI: AbstractDataSourceItem, AbstractDataItem
using Test, DataBrowserCache, DataBrowserAPI, Tables
using DataBrowserAPI: SOURCE_READ, SOURCE_INTERPRET, ITEM_PROCESS, ITEM_ANALYZE
using DataBrowserAPI.ItemIndex: ItemRecord

struct CacheSourceItem <: AbstractDataSourceItem end
DataBrowserAPI.id(::CacheSourceItem) = "source"
DataBrowserAPI.label(::CacheSourceItem) = "source"
DataBrowserAPI.fingerprint(::CacheSourceItem) = 1
struct CacheItem <: AbstractDataItem end

@testset "payload replacement and deletion survive reopening" begin
    mktempdir() do dir
        identity = ProjectCacheIdentity("test", "source", "source", joinpath(dir, "cache.duckdb"))
        cache = open_cache_db(identity)
        record = ItemRecord(id="item", label="item", type=CacheItem,
            source_item_key=source_item_key!(cache, "source"; mint=true))
        try
            write_meta_header!(cache)
            store_interpreted!(cache, CacheSourceItem(), "source", [record], [(x=[1, 2],)])
            store_processed!(cache, record, (x=[3, 4],))
            store_processed!(cache, record, (x=[5],))
            store_item_metadata!(cache, record, Dict(:total => 5))
        finally
            close_cache_db!(cache)
        end
        cache = open_cache_db(identity)
        try
            @test Tables.columntable(something(only(read_payload(cache, [record]; stage=ITEM_PROCESS)))).x == [5]
            @test load_cache_index(cache).item_metadata["item"][:total] == 5
            delete_source_output!(cache, record.source_item_key, [record])
        finally
            close_cache_db!(cache)
        end
        cache = open_cache_db(identity)
        try
            @test !has_payload(cache, "item"; stage=ITEM_PROCESS)
            @test isempty(load_cache_index(cache).source.items)
        finally
            close_cache_db!(cache)
        end
    end
end

@testset "read values distinguish nothing from a cache miss" begin
    mktempdir() do dir
        identity = ProjectCacheIdentity("read", "source", "source", joinpath(dir, "cache.duckdb"))
        cache = open_cache_db(identity)
        try
            key = source_item_key!(cache, "source"; mint=true)
            store_source_read!(cache, key, nothing)
            @test read_payload(cache, key) === Some(nothing)
            delete_source_item!(cache, key, ItemRecord[])
            @test read_payload(cache, key) === nothing
        finally
            close_cache_db!(cache)
        end
    end
end
