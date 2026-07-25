using DataBrowser
using Test
using Dates: Date, DateTime

const WIDE = DataBrowserCache
const WIDE_INDEX = DataBrowserAPI.ItemIndex

"""Open a fresh disk cache for one wide-table test and clean it up afterwards."""
function _with_wide_cache(work::Function, name::AbstractString)
    mktempdir() do dir
        source = DataBrowser.DirectorySource(dir)
        identity = WIDE.project_cache_identity("$(name)_$(basename(dir))", source)
        cache = WIDE.open_cache_db(identity)
        try
            work(cache)
        finally
            WIDE.close_cache_db!(cache)
            rm(dirname(identity.cache_path); force=true, recursive=true)
        end
    end
end

@testset "wide cache round-trips typed values" begin
    _with_wide_cache("WideRoundTrip") do cache
        store = cache.analyzed_collection_metadata
        WIDE.edit!(store, Int64(1), WIDE_INDEX.MetadataDict(
            :polarity => :positive,
            :day => Date(2026, 7, 5),
            :moment => DateTime(2026, 7, 5, 12, 0, 0),
            :flags => [true, false, true],
            :counts => Int64[1, 2, 3],
            :area => 3.5,
            :ok => true,
        ))
        loaded = read(store)[Int64(1)]
        @test loaded[:polarity] === :positive
        @test loaded[:day] == Date(2026, 7, 5)
        @test loaded[:day] isa Date
        @test loaded[:moment] == DateTime(2026, 7, 5, 12, 0, 0)
        @test loaded[:moment] isa DateTime
        @test loaded[:flags] == [true, false, true]
        @test loaded[:counts] == Int64[1, 2, 3]
        @test loaded[:area] == 3.5
        @test loaded[:ok] === true
    end
end

@testset "wide cache drops missing keys on reload" begin
    _with_wide_cache("WideMissing") do cache
        store = cache.analyzed_collection_metadata
        WIDE.edit!(store, Int64(1), WIDE_INDEX.MetadataDict(:present => 1, :absent => missing))
        loaded = read(store)[Int64(1)]
        @test loaded[:present] == 1
        @test !haskey(loaded, :absent)
    end
end

@testset "wide cache conflict surfaces as a dropped key" begin
    _with_wide_cache("WideConflict") do cache
        store = cache.analyzed_collection_metadata
        # Kind A registers :polarity as Float64.
        dropped_a = WIDE.edit!(store, Int64(1), WIDE_INDEX.MetadataDict(:polarity => 1.0, :extra => 2))
        @test isempty(dropped_a)
        # Kind B writes :polarity as a Symbol: the conflicting value drops, others survive.
        dropped_b = WIDE.edit!(
            store, Int64(2), WIDE_INDEX.MetadataDict(:polarity => :up, :other => 9))
        @test dropped_b == [(:polarity, WIDE.VT_FLOAT, WIDE.VT_SYMBOL)]
        # The surfaced message names the key and both types.
        messages = WIDE.store_collection_metadata!(
            cache, Int64(3), WIDE_INDEX.MetadataDict(:polarity => :up))
        @test messages == ["metadata :polarity expected Float64, got Symbol; value dropped"]
        loaded_b = read(store)[Int64(2)]
        @test !haskey(loaded_b, :polarity)
        @test loaded_b[:other] == 9
        # A's original Float64 value is untouched.
        @test read(store)[Int64(1)][:polarity] == 1.0
    end
end

@testset "wide cache read uses only committed columns from disk" begin
    _with_wide_cache("WidePendingColumns") do cache
        store = cache.source_item_metadata
        WIDE.edit!(store, Int64(1), WIDE_INDEX.MetadataDict(
            :timestamp => DateTime(2026, 1, 1), :V_base => 1.0))
        # `:cycle` is registered in memory but not ALTER TABLE'd until the next flush.
        WIDE.edit!(store, Int64(2), WIDE_INDEX.MetadataDict(:cycle => Int64(3)))
        loaded = read(store)
        @test loaded[Int64(1)][:V_base] == 1.0
        @test !haskey(loaded[Int64(1)], :cycle)
        @test loaded[Int64(2)][:cycle] == Int64(3)
    end
end
