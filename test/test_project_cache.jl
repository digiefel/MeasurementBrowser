using MeasurementBrowser
using MeasurementBrowser: ItemFailure
using Dates
using Test
import DuckDB, DBInterface

const ProjectCache = MeasurementBrowser.Cache

"""Persist a finished SourceScan through the live CacheDB streaming API, the way scanning does."""
function cache_source_scan!(cachedb, scan; data_by_id=Dict{String,Any}())
    ProjectCache.write_scan_identity!(cachedb)
    batches = Dict{String,Vector{MeasurementBrowser.ItemIndex.ItemRecord}}()
    for record in scan.hierarchy.all_items
        push!(
            get!(() -> MeasurementBrowser.ItemIndex.ItemRecord[], batches, record.source_item_id),
            record)
    end
    record_batches = collect(values(batches))
    data_batches = [Any[get(data_by_id, r.id, nothing) for r in batch] for batch in record_batches]
    written = isempty(record_batches) ? String[] :
        ProjectCache.reconcile_source_items!(cachedb, record_batches, data_batches)
    ProjectCache.finalize_scan!(cachedb, scan, Set(written))
    return nothing
end

@testset "project DuckDB cache" begin
    mktempdir() do dir
        source = test_source(TEST_PROJECT, dir)
        identity = ProjectCache.project_cache_identity("Cache path test", source)
        @test identity.cache_path == joinpath(
            DEPOT_PATH[1], "measurementbrowser", "Cache path test", "cache.duckdb")
        @test identity.project_name == "Cache path test"
    end

    for invalid_name in ("", ".", "..", "nested/project", "nested\\project")
        mktempdir() do dir
            source = test_source(TEST_PROJECT, dir)
            @test_throws ArgumentError ProjectCache.project_cache_identity(invalid_name, source)
        end
    end

    mktempdir() do first_root
        mktempdir() do second_root
            first = ProjectCache.ProjectCacheIdentity(
                "One project", first_root, "First", joinpath(first_root, "cache.duckdb"))
            cachedb = ProjectCache.open_cache_db(first)
            ProjectCache.write_scan_identity!(cachedb)
            ProjectCache.close_cache_db!(cachedb)

            second = ProjectCache.ProjectCacheIdentity(
                "One project", second_root, "Second", first.cache_path)
            other = ProjectCache.open_cache_db(second)
            try
                error = try
                    ProjectCache.load_cache_index(other)
                    nothing
                catch caught
                    caught
                end
                @test error isa ProjectCache.ProjectCacheError
                @test occursin("use Rebuild Cache", sprint(showerror, error))
                @test occursin(first_root, sprint(showerror, error))
                @test occursin(second_root, sprint(showerror, error))
            finally
                ProjectCache.close_cache_db!(other)
            end
        end
    end

    mktempdir() do dir
        write_test_source(joinpath(dir, "first.csv"))
        write_test_source(joinpath(dir, "second.csv"), 10)
        write(joinpath(dir, "metadata.txt"), "collection_path,wafer\ntest,A\n")
        source = scan_test_source(TEST_PROJECT, dir)
        adapter = test_source(TEST_PROJECT, dir)
        identity = ProjectCache.project_cache_identity("Cache round trip", adapter)
        ProjectCache._remove_cache_files(identity.cache_path)

        try
            cachedb = ProjectCache.open_cache_db(identity)
            try
                cache_source_scan!(cachedb, source)
                loaded = ProjectCache.load_cache_index(cachedb)
                @test length(loaded.source.hierarchy.all_items) == 2
                @test loaded.source.hierarchy.index[("test", "first")].parameters[:wafer] == "A"
                @test Set(item.id for item in loaded.source.hierarchy.all_items) ==
                      Set(item.id for item in source.hierarchy.all_items)
                @test all(item -> item.parameters[:wafer] == "A", loaded.source.hierarchy.all_items)
                # No item data was written, so an interactive read finds nothing cached.
                @test ProjectCache.read_cached_item_data(
                    cachedb, loaded.source.hierarchy.all_items) == Any[nothing, nothing]
            finally
                ProjectCache.close_cache_db!(cachedb)
            end

            # Freshness counts compare a fresh scan against the previously cached index.
            sleep(0.05)
            touch(first(source.hierarchy.all_items).source_item_path)
            rm(joinpath(dir, "second.csv"))
            write_test_source(joinpath(dir, "third.csv"), 20)
            updated_source = scan_test_source(TEST_PROJECT, dir)
            status = MeasurementBrowser.Workspace._post_write_status(
                identity, updated_source, ProjectCache.ProjectCacheIndex(identity, source))
            @test status.stale_source_items == 1
            @test status.new_source_items == 1
            @test status.deleted_source_items == 1

            settled = MeasurementBrowser.Workspace._post_write_status(
                identity, updated_source,
                ProjectCache.ProjectCacheIndex(identity, updated_source))
            @test settled.fresh_source_items == 2
            @test settled.stale_source_items == 0
            @test settled.new_source_items == 0
            @test settled.deleted_source_items == 0

            failed_source = MeasurementBrowser.SourceScan(
                updated_source.source_id, updated_source.source_label,
                updated_source.source_item_fingerprints, updated_source.hierarchy,
                [ItemFailure(
                    first(updated_source.hierarchy.all_items).source_item_id,
                    first(updated_source.hierarchy.all_items).id,
                    "synthetic analysis failure")])
            failed = ProjectCache.ProjectCacheIndex(identity, failed_source)
            @test only(values(failed.analysis_errors)) == "synthetic analysis failure"
        finally
            ProjectCache._remove_cache_files(identity.cache_path)
        end
    end

    @testset "metadata value round-trip and query" begin
        rec(id, params, stats) = MeasurementBrowser.ItemIndex.ItemRecord(;
            id=id, source_item_id="si_" * id, source_item_fingerprint=(id, 1),
            source_item_path="/tmp/" * id, item_label="L" * id, kind=:iv,
            collection=["grp"], parameters=params, stats=stats, item_fingerprint=(id, 2))
        r1 = rec("a",
            Dict(:flag => true, :n => Int64(7), :x => 3.5, :name => "hello",
                 :tag => :sym, :day => Date(2026, 6, 20)),
            Dict(:when => DateTime(2026, 6, 20, 12), :absent => missing,
                 :mask => Bool[true, false, true], :ints => Int64[1, 2],
                 :floats => [1.0, 2.5], :words => ["p", "q"]))
        r2 = rec("b", Dict(:n => Int64(99)), Dict(:x => 10.0))
        hierarchy = MeasurementBrowser.ItemIndex.Hierarchy("round_trip_src", false, 0)
        MeasurementBrowser.ItemIndex.insert_item!(hierarchy, r1)
        MeasurementBrowser.ItemIndex.insert_item!(hierarchy, r2)
        scan = MeasurementBrowser.SourceScan(
            "round_trip_src", "RoundTripSrc",
            Dict{String,Any}("si_a" => ("a", 1), "si_b" => ("b", 1)),
            hierarchy, ItemFailure[])

        mktempdir() do dir
            identity = ProjectCache.ProjectCacheIdentity(
                "rt", "round_trip_src", "RoundTripSrc", joinpath(dir, "rt.duckdb"))
            cachedb = ProjectCache.open_cache_db(identity)
            cache_source_scan!(cachedb, scan)
            loaded = ProjectCache.load_cache_index(cachedb)
            ProjectCache.close_cache_db!(cachedb)
            by_id = Dict(r.id => r for r in loaded.source.hierarchy.all_items)
            a = by_id["a"]

            @test a.parameters[:flag] === true
            @test a.parameters[:n] === Int64(7)
            @test a.parameters[:x] === 3.5
            @test a.parameters[:name] == "hello"
            @test a.parameters[:tag] === :sym
            @test a.parameters[:day] === Date(2026, 6, 20)
            @test a.stats[:when] === DateTime(2026, 6, 20, 12)
            @test a.stats[:absent] === missing
            @test a.stats[:mask] isa Vector{Bool} && a.stats[:mask] == Bool[true, false, true]
            @test a.stats[:ints] isa Vector{Int64} && a.stats[:ints] == Int64[1, 2]
            @test a.stats[:floats] isa Vector{Float64} && a.stats[:floats] == [1.0, 2.5]
            @test a.stats[:words] isa Vector{String} && a.stats[:words] == ["p", "q"]

            # The typed EAV columns are genuinely queryable.
            db = DBInterface.connect(DuckDB.DB, identity.cache_path)
            try
                statement = DBInterface.prepare(
                    db, "SELECT entity FROM metadata WHERE scope = 1 AND key = ? AND d > ?")
                hits = [row.entity for row in DBInterface.execute(statement, ("x", 5.0))]
                @test hits == ["b"]
            finally
                DBInterface.close!(db)
            end
        end
    end

    @testset "opening a project updates a stale cache" begin
        mktempdir() do dir
            write_test_source(joinpath(dir, "first.csv"))
            source = scan_test_source(TEST_PROJECT, dir)
            adapter = test_source(TEST_PROJECT, dir)
            identity = ProjectCache.project_cache_identity(
                MeasurementBrowser.project_name(TEST_PROJECT), adapter)
            ProjectCache._remove_cache_files(identity.cache_path)
            try
                stale_cache = ProjectCache.open_cache_db(identity)
                cache_source_scan!(stale_cache, source)
                ProjectCache.close_cache_db!(stale_cache)
                write_test_source(joinpath(dir, "second.csv"), 10)

                state = MeasurementBrowser.Browser.BrowserState(project_locked=true)
                MeasurementBrowser.Browser._open_project_path!(
                    state,
                    dir;
                    project=TEST_PROJECT,
                    persist=false,
                )
                saw_writing = false
                deadline = time() + 10
                while time() < deadline
                    workspace = state.workspace
                    MeasurementBrowser.Workspace.poll_workspace!(workspace)
                    saw_writing |= workspace.cache_state == :writing
                    status = workspace.cache.status
                    if saw_writing &&
                       workspace.cache_state == :ready &&
                       status isa ProjectCache.ProjectCacheStatus &&
                       status.cached_source_items == 2 &&
                       status.new_source_items == 1
                        break
                    end
                    sleep(0.02)
                end

                workspace = state.workspace
                status = workspace.cache.status
                @test saw_writing
                @test workspace.scan.state == :done
                @test workspace.cache_state == :ready
                @test status isa ProjectCache.ProjectCacheStatus
                @test status.cached_source_items == 2
                # second.csv is genuinely new relative to the stale one-item cache; first.csv was
                # reused without a re-read.
                @test status.new_source_items == 1
                @test status.stale_source_items == 0
                close_workspace!(workspace)
            finally
                ProjectCache._remove_cache_files(identity.cache_path)
            end
        end
    end

    @testset "incremental writes are durable, reconcilable, and stats persist" begin
        mkrec(id, sid) = MeasurementBrowser.ItemIndex.ItemRecord(;
            id=id, source_item_id=sid, source_item_fingerprint=(sid, 1),
            source_item_path="/tmp/" * sid, item_label="L" * id, kind=:iv,
            collection=["grp"], item_fingerprint=(id, 2))
        function make_scan(records)
            hierarchy = MeasurementBrowser.ItemIndex.Hierarchy("inc_src", false, 0)
            for record in records
                MeasurementBrowser.ItemIndex.insert_item!(hierarchy, record)
            end
            fingerprints = Dict{String,Any}(
                record.source_item_id => record.source_item_fingerprint for record in records)
            return MeasurementBrowser.SourceScan(
                "inc_src", "IncSrc", fingerprints, hierarchy, ItemFailure[])
        end
        ids(index) = Set(record.id for record in index.source.hierarchy.all_items)
        # Write one source item's records as scanning does (one bounded batch, no loaded data).
        reconcile_one!(cachedb, records) = ProjectCache.reconcile_source_items!(
            cachedb, [records], [Any[nothing for _ in records]])

        mktempdir() do dir
            identity = ProjectCache.ProjectCacheIdentity(
                "inc", "inc_src", "IncSrc", joinpath(dir, "inc.duckdb"))
            cachedb = ProjectCache.open_cache_db(identity)
            try
                recs_a = [mkrec("a1", "si_a"), mkrec("a2", "si_a")]
                recs_b = [mkrec("b1", "si_b")]

                # An interrupted scan (identity stamped, one source item written, no finalize)
                # is still loadable and shows the progress made so far.
                ProjectCache.write_scan_identity!(cachedb)
                reconcile_one!(cachedb, recs_a)
                partial = ProjectCache.load_cache_index(cachedb)
                @test ids(partial) == Set(["a1", "a2"])

                reconcile_one!(cachedb, recs_b)
                ProjectCache.finalize_scan!(
                    cachedb, make_scan([recs_a; recs_b]), Set(["si_a", "si_b"]))
                @test ids(ProjectCache.load_cache_index(cachedb)) == Set(["a1", "a2", "b1"])

                # Stats computed after the scan persist without any further index write.
                ProjectCache.persist_stats!(
                    cachedb,
                    Dict("a1" => Dict{Symbol,Any}(:vmax => 3.5)),
                    Dict{Tuple{Vararg{String}},Dict{Symbol,Any}}(),
                )
                reloaded = ProjectCache.load_cache_index(cachedb)
                a1 = only(r for r in reloaded.source.hierarchy.all_items if r.id == "a1")
                @test a1.stats[:vmax] === 3.5

                # Dropping si_a from the next scan cascades its items out of the cache.
                reconcile_one!(cachedb, recs_b)
                ProjectCache.finalize_scan!(cachedb, make_scan(recs_b), Set(["si_b"]))
                @test ids(ProjectCache.load_cache_index(cachedb)) == Set(["b1"])
            finally
                ProjectCache.close_cache_db!(cachedb)
            end
        end
    end

    @testset "a corrupt cache file is discarded and rebuilt on open" begin
        mktempdir() do dir
            path = joinpath(dir, "corrupt.duckdb")
            write(path, "this is not a duckdb file")
            identity = ProjectCache.ProjectCacheIdentity("c", "c_src", "CSrc", path)
            cachedb = ProjectCache.open_cache_db(identity)
            try
                @test cachedb isa ProjectCache.CacheDB
                ProjectCache.write_scan_identity!(cachedb)
                @test ProjectCache.load_cache_index(cachedb) isa ProjectCache.ProjectCacheIndex
            finally
                ProjectCache.close_cache_db!(cachedb)
            end

            connection = DBInterface.connect(DuckDB.DB, path)
            try
                DBInterface.execute(connection, """
                    UPDATE meta SET value = '1' WHERE key = 'schema_version'
                """)
            finally
                DBInterface.close!(connection)
            end
            rebuilt = ProjectCache.open_cache_db(identity)
            try
                meta_rows = ProjectCache.with_persistent_reader(rebuilt) do reader
                    only(DBInterface.execute(reader, "SELECT count(*) AS count FROM meta")).count
                end
                @test meta_rows == 0
            finally
                ProjectCache.close_cache_db!(rebuilt)
            end
        end
    end

    @testset "native DataFrame round-trip, validation, and skipping" begin
        prec(id, sif, iff) = MeasurementBrowser.ItemIndex.ItemRecord(;
            id=id, source_item_id="si_" * id, source_item_fingerprint=sif,
            source_item_path="/tmp/" * id, item_label="L" * id, kind=:iv,
            collection=["grp"], item_fingerprint=iff)
        mktempdir() do dir
            identity = ProjectCache.ProjectCacheIdentity(
                "pl", "pl_src", "PlSrc", joinpath(dir, "pl.duckdb"))
            cachedb = ProjectCache.open_cache_db(identity)
            try
                cacheable = prec("a", ("si_a", 1), ("a", 1))
                uncacheable = prec("b", nothing, nothing)
                source_frame = DataFrame(
                    voltage=[1.0, 2.0, 3.0],
                    current=Union{Missing,Float64}[0.1, missing, 0.3],
                )
                cached_frame = view(source_frame, 1:3, :)
                cached_item = MeasurementBrowser.DataItem(cacheable, cached_frame)
                @test MeasurementBrowser.cacheable(cached_item)
                ignored_item = MeasurementBrowser.DataItem(uncacheable, DataFrame(x=[1]))
                ProjectCache.write_cached_item_data!(
                    cachedb, [cacheable, uncacheable], Any[cached_item, ignored_item])

                # The cacheable item round-trips; the fingerprint-less one is never stored.
                got = ProjectCache.read_cached_item_data(cachedb, [cacheable, uncacheable])
                @test got[1] isa MeasurementBrowser.DataItem
                @test isequal(MeasurementBrowser.item_data(got[1]), DataFrame(cached_frame))
                @test MeasurementBrowser.item_data(got[1]) isa DataFrame
                @test got[2] === nothing

                rows = ProjectCache.with_persistent_reader(cachedb) do connection
                    collect(DBInterface.execute(connection, """
                        SELECT storage_id, row_count FROM item_data
                    """))
                end
                @test length(rows) == 1
                @test only(rows).row_count == 3

                # A changed item fingerprint invalidates the stored data.
                restamped = prec("a", ("si_a", 1), ("a", 2))
                @test only(ProjectCache.read_cached_item_data(cachedb, [restamped])) === nothing
            finally
                ProjectCache.close_cache_db!(cachedb)
            end
        end
    end

    @testset "bulk native DataFrame read" begin
        mktempdir() do dir
            identity = ProjectCache.ProjectCacheIdentity(
                "bulk", "bulk_src", "BulkSrc", joinpath(dir, "bulk.duckdb"))
            cachedb = ProjectCache.open_cache_db(identity)
            try
                records = MeasurementBrowser.ItemIndex.ItemRecord[
                    MeasurementBrowser.ItemIndex.ItemRecord(;
                        id="item_$index",
                        source_item_id="source_$index",
                        source_item_fingerprint=(index, 1),
                        source_item_path=nothing,
                        item_label="Item $index",
                        kind=:table,
                        collection=["bulk"],
                    )
                    for index in 1:17
                ]
                frames = [DataFrame(x=[Float64(index), index + 0.5]) for index in 1:17]
                items = Any[
                    MeasurementBrowser.DataItem(record, frame)
                    for (record, frame) in zip(records, frames)
                ]
                ProjectCache.write_cached_item_data!(cachedb, records, items)
                loaded = ProjectCache.read_cached_item_data(cachedb, records)

                @test length(loaded) == length(records)
                @test all(
                    isequal(MeasurementBrowser.item_data(item), frame)
                    for (item, frame) in zip(loaded, frames)
                )
            finally
                ProjectCache.close_cache_db!(cachedb)
            end
        end
    end

    @testset "item LRU evicts least-recently-used entries" begin
        Cache = MeasurementBrowser.Workspace
        lru = Cache.ItemCache(2)
        Cache.cache_item!(lru, "a", :fp, "A")
        Cache.cache_item!(lru, "b", :fp, "B")
        @test Cache.lookup_item(lru, "a", nothing) == (:fp, "A")   # refreshes "a" as most-recent
        Cache.cache_item!(lru, "c", :fp, "C")                      # capacity 2 → evict "b" (oldest)
        @test Cache.lookup_item(lru, "b", nothing) === nothing
        @test Cache.lookup_item(lru, "a", nothing) == (:fp, "A")
        @test Cache.lookup_item(lru, "c", nothing) == (:fp, "C")
    end

    @testset "materialization is served from the item-data cache without the origin" begin
        mktempdir() do dir
            write_test_source(joinpath(dir, "first.csv"))
            records = scan_test_source(TEST_PROJECT, dir).hierarchy.all_items
            workspace = MeasurementBrowser.Workspace.Workspace(
                TEST_PROJECT, test_source(TEST_PROJECT, dir))
            try
                first_pass = read_item_data(workspace, records)
                @test all(!isnothing, first_pass)

                # Populate the durable item-data cache the way source scanning does (the interactive
                # read path itself never writes), then drop the origin file and the in-memory cache.
                materialized = MeasurementBrowser.Workspace.materialize_items(workspace, records)
                ProjectCache.write_cached_item_data!(workspace.cache.db, records, materialized)
                rm(joinpath(dir, "first.csv"))
                workspace.loaded_items = MeasurementBrowser.Workspace.ItemCache()

                second_pass = read_item_data(workspace, records)
                @test second_pass == first_pass
            finally
                close_workspace!(workspace)
            end
        end
    end
end
