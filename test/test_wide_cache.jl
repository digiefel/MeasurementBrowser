using MeasurementBrowser
using Test
using Dates: Date, DateTime

const WIDE = MeasurementBrowser.Cache
const WIDE_INDEX = MeasurementBrowser.ItemIndex

"""Open a fresh disk cache for one wide-table test and clean it up afterwards."""
function _with_wide_cache(work::Function, name::AbstractString)
    mktempdir() do dir
        source = MeasurementBrowser.DirectorySource(dir)
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
        store = cache.collection_metadata
        WIDE.edit!(store, "c1", WIDE_INDEX.MetadataDict(
            :polarity => :positive,
            :day => Date(2026, 7, 5),
            :moment => DateTime(2026, 7, 5, 12, 0, 0),
            :flags => [true, false, true],
            :counts => Int64[1, 2, 3],
            :area => 3.5,
            :ok => true,
        ))
        loaded = read(store)["c1"]
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
        store = cache.collection_metadata
        WIDE.edit!(store, "c1", WIDE_INDEX.MetadataDict(:present => 1, :absent => missing))
        loaded = read(store)["c1"]
        @test loaded[:present] == 1
        @test !haskey(loaded, :absent)
    end
end

@testset "wide cache conflict surfaces as a dropped key" begin
    _with_wide_cache("WideConflict") do cache
        store = cache.collection_metadata
        # Kind A registers :polarity as Float64.
        dropped_a = WIDE.edit!(store, "a", WIDE_INDEX.MetadataDict(:polarity => 1.0, :extra => 2))
        @test isempty(dropped_a)
        # Kind B writes :polarity as a Symbol: the conflicting value drops, others survive.
        dropped_b = WIDE.edit!(
            store, "b", WIDE_INDEX.MetadataDict(:polarity => :up, :other => 9))
        @test dropped_b == [(:polarity, WIDE.VT_FLOAT, WIDE.VT_SYMBOL)]
        # The surfaced message names the key and both types.
        messages = WIDE.store_collection_metadata!(
            cache, "c", WIDE_INDEX.MetadataDict(:polarity => :up))
        @test messages == ["metadata :polarity expected Float64, got Symbol; value dropped"]
        loaded_b = read(store)["b"]
        @test !haskey(loaded_b, :polarity)
        @test loaded_b[:other] == 9
        # A's original Float64 value is untouched.
        @test read(store)["a"][:polarity] == 1.0
    end
end

@testset "query_items filters on committed effective metadata" begin
    mktempdir() do dir
        write(joinpath(dir, "hot.csv"), "x\n5\n")
        write(joinpath(dir, "cold.csv"), "x\n1\n")
        write(joinpath(dir, "metadata.txt"), "collection_path,wafer\ndev,A\n")

        project = MeasurementBrowser.define_project("QueryItems_$(basename(dir))")
        MeasurementBrowser.register_item!(project, :table;
            detect=file -> endswith(file.filename, ".csv"),
            read=file -> read(file.filepath, String),
            entries=(file, _data) -> [MeasurementBrowser.DataItem(
                kind=:table,
                collection=["dev", splitext(file.filename)[1]],
                data=nothing)],
            analyze=item -> Dict{Symbol,Any}(
                :peak => splitext(basename(item.id))[1] == "hot" ? 5.0 : 1.0),
        )
        workspace = MeasurementBrowser.open_workspace(
            project, MeasurementBrowser.DirectorySource(dir); background_processing=true)
        try
            wait_workspace_idle!(workspace)
            # Give the buffered wide-column writes their flush window before querying committed DB;
            # until the ALTER commits, the predicate column is not yet present.
            settled() = try
                length(MeasurementBrowser.query_items(workspace, "peak > 4")) == 1
            catch
                false
            end
            @test Base.timedwait(settled, 10) === :ok
            matched = MeasurementBrowser.query_items(workspace, "peak > 4 AND wafer = 'A'")
            @test length(matched) == 1
            @test occursin("hot", only(matched))
        finally
            MeasurementBrowser.close_workspace!(workspace)
        end
    end
end

@testset "memory cache answers query_items with an empty vector" begin
    mktempdir() do dir
        write(joinpath(dir, "a.csv"), "x\n1\n")
        workspace = MeasurementBrowser.open_workspace(
            TEST_PROJECT, MeasurementBrowser.DirectorySource(dir); cache=false)
        try
            wait_workspace_idle!(workspace)
            @test MeasurementBrowser.query_items(workspace, "1 = 1") == String[]
        finally
            MeasurementBrowser.close_workspace!(workspace)
        end
    end
end

@testset "wide cache widens under concurrent reads" begin
    _with_wide_cache("WideConcurrent") do cache
        store = cache.collection_metadata
        stop = Base.Threads.Atomic{Bool}(false)
        reader = Base.Threads.@spawn while !stop[]
            read(store)
            yield()
        end
        for index in 1:40
            WIDE.edit!(store, "c$index", WIDE_INDEX.MetadataDict(
                Symbol("col$index") => index))
        end
        stop[] = true
        wait(reader)
        loaded = read(store)
        @test loaded["c40"][Symbol("col40")] == 40
    end
end
