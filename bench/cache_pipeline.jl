#!/usr/bin/env julia
# Headless benchmark of the core pipeline: read -> process -> item-analyze -> cache write/read.
#
# Run:  julia --project=bench --threads=auto bench/cache_pipeline.jl [n_files] [rows_per_file]
#
# It generates a synthetic dataset, then reports:
#   1. micro-benchmarks that settle "is serialization the bottleneck?" (parse vs serialize vs
#      deserialize vs DuckDB blob round-trip vs DuckDB columnar round-trip),
#   2. macro timings for the real scan (parallel) and a single-threaded analysis loop,
#   3. the engine's own per-region TimerOutputs report for a representative cache exercise.
#
# Everything here measures real functions on real data — no estimates.

using MeasurementBrowser
const MB = MeasurementBrowser
using DataFrames, CSV, Random, Printf, Statistics, Serialization
using BenchmarkTools
import DuckDB, DBInterface

const N_FILES = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 500
const ROWS = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 1000

section(title) = (println("\n", "="^78); println(title); println("="^78))
ms(seconds) = @sprintf("%.3f ms", seconds * 1e3)
kb(bytes) = @sprintf("%.1f KB", bytes / 1024)

# --------------------------------------------------------------------------------------------------
# Synthetic dataset: many small CSVs shaped like an I-V sweep.
# --------------------------------------------------------------------------------------------------
function generate_dataset(dir::AbstractString; n_files::Int, rows::Int)
    rng = Random.MersenneTwister(0x5eed)
    for i in 1:n_files
        v = collect(range(-1.0, 1.0; length=rows))
        df = DataFrame(v=v, i=(v .^ 3) .+ 0.01 .* randn(rng, rows), t=cumsum(rand(rng, rows)))
        CSV.write(joinpath(dir, "meas_$(lpad(i, 5, '0')).csv"), df)
    end
    return nothing
end

function bench_project()
    project = MB.define_project("Bench")
    MB.register_item!(
        project, :iv;
        detect=file -> endswith(file.filename, ".csv"),
        read=file -> CSV.read(file.filepath, DataFrame),
        entries=(file, data) -> [MB.DataItem(
            kind=:iv,
            collection=[splitext(file.filename)[1]],
            parameters=Dict{Symbol,Any}(),
            data=data,
        )],
        stats=item -> Dict{Symbol,Any}(
            :vmax => maximum(item.data.v),
            :imax => maximum(item.data.i),
            :imean => Statistics.mean(item.data.i),
        ),
    )
    return project
end

# --------------------------------------------------------------------------------------------------
# 1. Micro-benchmarks on one representative file/DataFrame.
# --------------------------------------------------------------------------------------------------
function micro_benchmarks(sample_csv::AbstractString)
    section("1. PER-ITEM MICRO-BENCHMARKS (one $(ROWS)-row file)")
    df = CSV.read(sample_csv, DataFrame)
    bytes = serialize(IOBuffer(), df)  # warmup
    blob = (io = IOBuffer(); serialize(io, df); take!(io))

    t_parse = @belapsed CSV.read($sample_csv, DataFrame)
    t_ser = @belapsed (io = IOBuffer(); serialize(io, $df); take!(io))
    t_deser = @belapsed deserialize(IOBuffer($blob))

    # DuckDB blob round-trip (what the cache does today: serialize -> BLOB -> deserialize).
    con = DBInterface.connect(DuckDB.DB)
    DBInterface.execute(con, "CREATE TABLE b(id INTEGER, blob BLOB)")
    ins = DBInterface.prepare(con, "INSERT INTO b VALUES (?, ?)")
    sel = DBInterface.prepare(con, "SELECT blob FROM b WHERE id = ?")
    t_blob_write = @belapsed DBInterface.execute($ins, (1, (io = IOBuffer(); serialize(io, $df); take!(io))))
    # The write benchmark left many rows behind; keep exactly one so the read benchmark reads one blob.
    DBInterface.execute(con, "DELETE FROM b")
    DBInterface.execute(ins, (1, blob))
    t_blob_read = @belapsed begin
        local out = nothing
        for row in DBInterface.execute($sel, (1,))
            out = deserialize(IOBuffer(Vector{UInt8}(row.blob)))
        end
        out
    end

    # DuckDB columnar round-trip (store the DataFrame's columns natively; no Julia serialization).
    DuckDB.register_data_frame(con, df, "src")
    t_col_write = @belapsed begin
        DBInterface.execute($con, "DROP TABLE IF EXISTS cols")
        DBInterface.execute($con, "CREATE TABLE cols AS SELECT * FROM src")
    end
    DBInterface.execute(con, "DROP TABLE IF EXISTS cols")
    DBInterface.execute(con, "CREATE TABLE cols AS SELECT * FROM src")
    t_col_read = @belapsed DataFrame(DBInterface.execute($con, "SELECT * FROM cols"))
    DBInterface.close!(con)

    @printf("  blob size (serialized DataFrame):     %s\n", kb(length(blob)))
    println()
    @printf("  parse CSV from disk (read once):      %s\n", ms(t_parse))
    println("  --- if we cache the parsed payload, a later read costs one of: ---")
    @printf("  serialize DataFrame -> bytes:         %s\n", ms(t_ser))
    @printf("  deserialize bytes -> DataFrame:       %s\n", ms(t_deser))
    @printf("  DuckDB blob write (serialize+insert): %s\n", ms(t_blob_write))
    @printf("  DuckDB blob read  (select+deserialize): %s\n", ms(t_blob_read))
    @printf("  DuckDB columnar write (CREATE TABLE): %s\n", ms(t_col_write))
    @printf("  DuckDB columnar read  (SELECT *):     %s\n", ms(t_col_read))
    println()
    @printf("  VERDICT: cached read via blob = %s vs re-parse = %s  (speedup %.2fx)\n",
        ms(t_blob_read), ms(t_parse), t_parse / t_blob_read)
    @printf("           cached read via columnar = %s  (speedup %.2fx)\n",
        ms(t_col_read), t_parse / t_col_read)
    return nothing
end

# --------------------------------------------------------------------------------------------------
# 2. Macro timings on the whole dataset.
# --------------------------------------------------------------------------------------------------
function macro_timings(project, dir::AbstractString)
    section("2. WHOLE-DATASET MACRO TIMINGS ($(N_FILES) files x $(ROWS) rows)")
    source = MB.DirectorySource(dir)

    MB.reset_scan_profile!(project)
    MB.scan_source(project, source)              # warmup / compile
    MB.reset_scan_profile!(project)
    t_scan = @elapsed scan = MB.scan_source(project, source)
    records = scan.hierarchy.all_items
    @printf("  scan (parallel, %d threads): %s for %d items  (%s/item)\n",
        Threads.nthreads(), ms(t_scan), length(records), ms(t_scan / length(records)))
    for row in MB.scan_profile_summary(project)
        @printf("    kind=%s  read=%s  stats=%s  items=%d\n",
            row.kind, ms(row.read_seconds), ms(row.stats_seconds), row.items)
    end

    # Single-threaded analysis loop: load each item from origin + compute stats (today's behavior).
    function analysis_loop()
        total = 0.0
        for record in records
            item = MB.load_data_item(project, source, record.source_item_id, record.id)
            s = MB.Projects.compute_item_stats(project, source, item)
            total += s[:imean]
        end
        return total
    end
    analysis_loop()                              # warmup
    t_analysis = @elapsed analysis_loop()
    @printf("  analysis (single-thread): %s for %d items  (%s/item)\n",
        ms(t_analysis), length(records), ms(t_analysis / length(records)))
    return scan, source, records
end

# --------------------------------------------------------------------------------------------------
# 3. Engine per-region timings (TimerOutputs) over a representative cache exercise.
# --------------------------------------------------------------------------------------------------
function region_timings(project, source, records)
    section("3. ENGINE PER-REGION TIMINGS (cache write-all then read-all)")
    identity = MB.Cache.project_cache_identity(MB.Cache.project_cache_id(source), source)
    MB.Cache._remove_cache_files(identity.cache_path)
    cachedb = MB.Cache.open_cache_db(identity)
    try
        # Materialize every item from origin once, then persist all payloads, then read them back.
        items = [MB.load_data_item(project, source, r.source_item_id, r.id) for r in records]

        MB.Profiling.reset!()   # timings already enabled at top level (avoids a world-age miss)
        MB.Cache.write_item_payloads!(cachedb, records, items)   # serialize + insert (timed)
        got = MB.Cache.read_item_payloads(cachedb, records)      # select + deserialize (timed)
        MB.Profiling.report()
        @printf("  payloads read back: %d/%d non-nothing\n", count(!isnothing, got), length(records))
    finally
        MB.Cache.close_cache_db!(cachedb)
        MB.Cache._remove_cache_files(identity.cache_path)
    end
    return nothing
end

function main()
    println("MeasurementBrowser pipeline benchmark")
    @printf("threads=%d  n_files=%d  rows=%d\n", Threads.nthreads(), N_FILES, ROWS)
    mktempdir() do dir
        print("generating dataset... ")
        @printf("%s\n", ms(@elapsed generate_dataset(dir; n_files=N_FILES, rows=ROWS)))
        project = bench_project()
        sample = joinpath(dir, "meas_00001.csv")
        micro_benchmarks(sample)
        scan, source, records = macro_timings(project, dir)
        region_timings(project, source, records)
    end
    return nothing
end

# Enable engine timings here at top level so the new method definitions are in a newer world age than
# the workload calls below (enabling inside the workload function would not take effect — world age).
MB.Profiling.enable!()
try
    main()
finally
    MB.Profiling.disable!()
end
