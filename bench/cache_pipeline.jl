#!/usr/bin/env julia
# Headless benchmark of the core pipeline: read -> process -> item-analyze -> cache write/read.
#
# Run:  julia --project=bench --threads=auto bench/cache_pipeline.jl [n_files] [rows_per_file] [workspace]
#
# It reuses a synthetic dataset under bench/data, then reports:
#   1. micro-benchmarks that settle "is serialization the bottleneck?" (parse vs serialize vs
#      deserialize vs DuckDB blob round-trip vs DuckDB columnar round-trip),
#   2. macro timings for the real scan (parallel) and a single-threaded analysis loop,
#   3. the engine's own per-region TimerOutputs report for a representative cache exercise.
#
# Everything here measures real functions on real data — no estimates. Reports are saved under
# bench/results/<datetime> <commit> results.txt.

using MeasurementBrowser
const MB = MeasurementBrowser
using DataFrames, CSV, Dates, Random, Printf, Statistics, Serialization
using BenchmarkTools
import DuckDB, DBInterface

const N_FILES = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 500
const ROWS = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 1000
const WORKSPACE_ONLY = length(ARGS) >= 3 && ARGS[3] == "workspace"
const EXPANDED_SOURCE_ITEMS = max(1, N_FILES ÷ 10)
const ITEMS_PER_EXPANDED_SOURCE = 100
const BENCH_DIR = @__DIR__
const ROOT_DIR = normpath(joinpath(BENCH_DIR, ".."))
const DATA_DIR = joinpath(BENCH_DIR, "data", "cache_pipeline_$(N_FILES)x$(ROWS)")
const RESULTS_DIR = joinpath(BENCH_DIR, "results")
const OUT = Ref{IO}(stdout)

ms(seconds) = @sprintf("%.3f ms", seconds * 1e3)
kb(bytes) = @sprintf("%.1f KB", bytes / 1024)
mib(bytes) = @sprintf("%.1f MiB", bytes / 1024^2)

struct TeeIO <: IO
    a::IO
    b::IO
end

Base.iswritable(::TeeIO) = true
Base.flush(io::TeeIO) = (flush(io.a); flush(io.b))
Base.write(io::TeeIO, byte::UInt8) = (write(io.a, byte); write(io.b, byte); 1)
Base.write(io::TeeIO, bytes::AbstractVector{UInt8}) =
    (n = write(io.a, bytes); write(io.b, bytes); n)
Base.write(io::TeeIO, text::Union{String,SubString{String}}) =
    (n = write(io.a, text); write(io.b, text); n)

out() = OUT[]
section(title) = (println(out(), "\n", "="^78); println(out(), title); println(out(), "="^78))
say(args...) = println(out(), args...)
sayf(fmt::AbstractString, args...) = print(out(), Printf.format(Printf.Format(fmt), args...))

function git_output(args::AbstractString...)::String
    try
        return chomp(read(Cmd(["git", "-C", ROOT_DIR, args...]), String))
    catch
        return "unknown"
    end
end

function commit_label()::String
    commit = git_output("rev-parse", "--short", "HEAD")
    status = git_output("status", "--short")
    return isempty(status) || status == "unknown" ? commit : commit * "-dirty"
end

function result_path()::String
    mkpath(RESULTS_DIR)
    stamp = Dates.format(now(), dateformat"yyyy-mm-dd HHMMSS")
    return joinpath(RESULTS_DIR, "$stamp $(commit_label()) results.txt")
end

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

function dataset_ready(dir::AbstractString; n_files::Int, rows::Int)::Bool
    marker = joinpath(dir, "COMPLETE.txt")
    first_file = joinpath(dir, "meas_00001.csv")
    last_file = joinpath(dir, "meas_$(lpad(n_files, 5, '0')).csv")
    return isfile(marker) && isfile(first_file) && isfile(last_file)
end

function ensure_dataset(dir::AbstractString; n_files::Int, rows::Int)::Nothing
    if dataset_ready(dir; n_files, rows)
        sayf("using dataset: %s\n", dir)
        return nothing
    end
    rm(dir; force=true, recursive=true)
    mkpath(dir)
    sayf("generating dataset: %s\n", dir)
    sayf("  %s\n", ms(@elapsed generate_dataset(dir; n_files, rows)))
    write(
        joinpath(dir, "COMPLETE.txt"),
        "n_files=$n_files\nrows=$rows\nseed=0x5eed\n",
    )
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

"""
Build a project where ten percent of source items expand into 100 data items.

Expanded items copy their columns into separate `DataFrame`s, matching ordinary project code that
constructs one table per logical item. The read counter makes repeated origin reads visible instead
of inferring them from elapsed time.
"""
function expansion_project()
    expanded_filenames = Set(
        "meas_$(lpad(index, 5, '0')).csv"
        for index in 1:EXPANDED_SOURCE_ITEMS
    )
    read_count = Threads.Atomic{Int}(0)
    project = MB.define_project("ExpansionBench")
    MB.register_item!(
        project, :iv;
        detect=file -> endswith(file.filename, ".csv"),
        read=function (file)
            Threads.atomic_add!(read_count, 1)
            return CSV.read(file.filepath, DataFrame)
        end,
        entries=function (file, data)
            file.filename in expanded_filenames || return [MB.DataItem(
                kind=:iv,
                collection=[splitext(file.filename)[1]],
                parameters=Dict{Symbol,Any}(),
                data=data,
            )]

            row_count = nrow(data)
            items = MB.AbstractDataItem[]
            sizehint!(items, ITEMS_PER_EXPANDED_SOURCE)
            for part in 1:ITEMS_PER_EXPANDED_SOURCE
                first_row = fld((part - 1) * row_count, ITEMS_PER_EXPANDED_SOURCE) + 1
                last_row = fld(part * row_count, ITEMS_PER_EXPANDED_SOURCE)
                rows = first_row:last_row
                push!(items, MB.DataItem(
                    kind=:iv,
                    collection=[splitext(file.filename)[1]],
                    parameters=Dict{Symbol,Any}(:part => part),
                    data=DataFrame(
                        v=data.v[rows],
                        i=data.i[rows],
                        t=data.t[rows],
                    ),
                ))
            end
            return items
        end,
        stats=item -> Dict{Symbol,Any}(
            :rows => nrow(item.data),
            :imean => Statistics.mean(item.data.i),
        ),
    )
    return project, read_count
end

# --------------------------------------------------------------------------------------------------
# 1. Micro-benchmarks on one representative file/DataFrame.
# --------------------------------------------------------------------------------------------------
function micro_benchmarks(sample_csv::AbstractString)
    section("1. PER-ITEM MICRO-BENCHMARKS (one $(ROWS)-row file)")
    df = CSV.read(sample_csv, DataFrame)
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

    sayf("  blob size (serialized DataFrame):     %s\n", kb(length(blob)))
    say()
    sayf("  parse CSV from disk (read once):      %s\n", ms(t_parse))
    say("  --- if we cache the parsed payload, a later read costs one of: ---")
    sayf("  serialize DataFrame -> bytes:         %s\n", ms(t_ser))
    sayf("  deserialize bytes -> DataFrame:       %s\n", ms(t_deser))
    sayf("  DuckDB blob write (serialize+insert): %s\n", ms(t_blob_write))
    sayf("  DuckDB blob read  (select+deserialize): %s\n", ms(t_blob_read))
    sayf("  DuckDB columnar write (CREATE TABLE): %s\n", ms(t_col_write))
    sayf("  DuckDB columnar read  (SELECT *):     %s\n", ms(t_col_read))
    say()
    sayf("  VERDICT: cached read via blob = %s vs re-parse = %s  (speedup %.2fx)\n",
        ms(t_blob_read), ms(t_parse), t_parse / t_blob_read)
    sayf("           cached read via columnar = %s  (speedup %.2fx)\n",
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
    sayf("  scan (parallel, %d threads): %s for %d items  (%s/item)\n",
        Threads.nthreads(), ms(t_scan), length(records), ms(t_scan / length(records)))
    for row in MB.scan_profile_summary(project)
        sayf("    kind=%s  read=%s  stats=%s  items=%d\n",
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
    sayf("  analysis (single-thread): %s for %d items  (%s/item)\n",
        ms(t_analysis), length(records), ms(t_analysis / length(records)))
    return scan, source, records
end

# --------------------------------------------------------------------------------------------------
# 3. Source-item expansion: current per-item reloads versus a read-once source pass.
# --------------------------------------------------------------------------------------------------
function expansion_timings(dir::AbstractString)
    section(
        "3. SOURCE-ITEM EXPANSION " *
        "($(EXPANDED_SOURCE_ITEMS)/$(N_FILES) sources expand to " *
        "$(ITEMS_PER_EXPANDED_SOURCE) data items)",
    )
    project, read_count = expansion_project()
    source = MB.DirectorySource(dir)

    # Compile the scan path before measuring it.
    MB.scan_source(project, source)
    Threads.atomic_xchg!(read_count, 0)
    scan_timing = @timed scan = MB.scan_source(project, source)
    scan_reads = read_count[]
    records = scan.hierarchy.all_items

    function per_item_analysis()
        total = 0.0
        for record in records
            item = MB.load_data_item(project, source, record.source_item_id, record.id)
            total += MB.Projects.compute_item_stats(project, source, item)[:imean]
        end
        return total
    end

    per_item_analysis()
    Threads.atomic_xchg!(read_count, 0)
    per_item_timing = @timed per_item_analysis()
    per_item_reads = read_count[]

    discovered = MB.source_items(source)
    function source_pass_analysis()
        total = 0.0
        for source_item in discovered
            for item in MB.data_items(project, source, source_item)
                total += MB.Projects.compute_item_stats(project, source, item)[:imean]
            end
        end
        return total
    end

    source_pass_analysis()
    Threads.atomic_xchg!(read_count, 0)
    source_pass_timing = @timed source_pass_analysis()
    source_pass_reads = read_count[]

    expected_items =
        EXPANDED_SOURCE_ITEMS * ITEMS_PER_EXPANDED_SOURCE +
        (N_FILES - EXPANDED_SOURCE_ITEMS)
    length(records) == expected_items || error(
        "Expansion benchmark produced $(length(records)) data items; expected $expected_items",
    )
    scan_reads == N_FILES || error(
        "Expansion scan read $scan_reads source items; expected exactly $N_FILES",
    )
    source_pass_reads == N_FILES || error(
        "Read-once analysis read $source_pass_reads source items; expected exactly $N_FILES",
    )

    sayf("  scan:                       %s  alloc=%s  gc=%s  source reads=%d  data items=%d\n",
        ms(scan_timing.time), mib(scan_timing.bytes), ms(scan_timing.gctime),
        scan_reads, length(records))
    sayf("  current per-item analysis:  %s  alloc=%s  gc=%s  source reads=%d\n",
        ms(per_item_timing.time), mib(per_item_timing.bytes), ms(per_item_timing.gctime),
        per_item_reads)
    sayf("  read-once source pass:      %s  alloc=%s  gc=%s  source reads=%d\n",
        ms(source_pass_timing.time), mib(source_pass_timing.bytes), ms(source_pass_timing.gctime),
        source_pass_reads)
    sayf("  legacy two-pass total:      %s  alloc=%s  source reads=%d\n",
        ms(scan_timing.time + per_item_timing.time),
        mib(scan_timing.bytes + per_item_timing.bytes),
        scan_reads + per_item_reads)
    sayf("  read-once pipeline gain:    %.2fx time, %.2fx allocation\n",
        (scan_timing.time + per_item_timing.time) / scan_timing.time,
        (scan_timing.bytes + per_item_timing.bytes) / scan_timing.bytes)
    sayf("  isolated analysis gain:     %.2fx time, %.2fx allocation\n",
        per_item_timing.time / source_pass_timing.time,
        per_item_timing.bytes / source_pass_timing.bytes)
    return nothing
end

# --------------------------------------------------------------------------------------------------
# 4. Full workspace cold build: scan, stats, cache writes, then cache-only materialization.
# --------------------------------------------------------------------------------------------------
function workspace_cold_build(dir::AbstractString)
    section("4. WORKSPACE COLD BUILD (expanded dataset, including DuckDB writes)")
    project, read_count = expansion_project()
    source = MB.DirectorySource(dir)
    identity = MB.Cache.project_cache_identity(MB.Cache.project_cache_id(source), source)
    MB.Cache._remove_cache_files(identity.cache_path)
    workspace = MB.Workspace.Workspace(project, source)

    try
        Threads.atomic_xchg!(read_count, 0)
        MB.Profiling.reset!()
        started = time_ns()
        MB.Workspace.scan_source!(workspace; rebuild=true)
        deadline = time() + 120
        while time() < deadline
            MB.Workspace.poll_workspace!(workspace)
            workspace.scan.state == :error &&
                error("Workspace benchmark scan failed: $(workspace.scan.error)")
            workspace.analysis.state == :error &&
                error("Workspace benchmark analysis failed: $(workspace.analysis.error)")
            workspace.scan.state == :done && workspace.analysis.state == :done && break
            sleep(0.005)
        end
        workspace.scan.state == :done || error(
            "Workspace benchmark scan did not finish within 120 seconds",
        )
        workspace.analysis.state == :done || error(
            "Workspace benchmark analysis did not finish within 120 seconds",
        )
        build_seconds = (time_ns() - started) / 1e9
        build_reads = read_count[]
        records = workspace.index.hierarchy.all_items

        materialize_started = time_ns()
        loaded = MB.Workspace.materialize_items(workspace, records)
        materialize_seconds = (time_ns() - materialize_started) / 1e9
        materialize_reads = read_count[] - build_reads
        cache_size = isfile(identity.cache_path) ? filesize(identity.cache_path) : 0

        build_reads == N_FILES || error(
            "Workspace build read $build_reads source items; expected exactly $N_FILES",
        )
        materialize_reads == 0 || error(
            "Cache materialization performed $materialize_reads unexpected source reads",
        )
        length(loaded) == length(records) || error(
            "Materialized $(length(loaded)) items; expected $(length(records))",
        )

        sayf("  cold build:                 %s  source reads=%d  data items=%d\n",
            ms(build_seconds), build_reads, length(records))
        sayf("  cache-only materialization: %s  source reads=%d\n",
            ms(materialize_seconds), materialize_reads)
        sayf("  DuckDB file size:           %s\n", mib(cache_size))
        MB.Profiling.report(out())
    finally
        MB.close_workspace!(workspace)
        MB.Cache._remove_cache_files(identity.cache_path)
    end
    return nothing
end

# --------------------------------------------------------------------------------------------------
# 5. Engine per-region timings (TimerOutputs) over a representative cache exercise.
# --------------------------------------------------------------------------------------------------
function region_timings(project, source, records)
    section("5. ENGINE PER-REGION TIMINGS (cache write-all then read-all)")
    identity = MB.Cache.project_cache_identity(MB.Cache.project_cache_id(source), source)
    MB.Cache._remove_cache_files(identity.cache_path)
    cachedb = MB.Cache.open_cache_db(identity)
    try
        # Materialize every item from origin once, then persist all payloads, then read them back.
        items = [MB.load_data_item(project, source, r.source_item_id, r.id) for r in records]

        MB.Profiling.reset!()   # timings already enabled at top level (avoids a world-age miss)
        MB.Cache.write_item_payloads!(cachedb, records, items)   # serialize + insert (timed)
        got = MB.Cache.read_item_payloads(cachedb, records)      # select + deserialize (timed)
        MB.Profiling.report(out())
        sayf("  payloads read back: %d/%d non-nothing\n", count(!isnothing, got), length(records))
    finally
        MB.Cache.close_cache_db!(cachedb)
        MB.Cache._remove_cache_files(identity.cache_path)
    end
    return nothing
end

function main(report_path::AbstractString)
    say("MeasurementBrowser pipeline benchmark")
    sayf("threads=%d  n_files=%d  rows=%d\n", Threads.nthreads(), N_FILES, ROWS)
    sayf("commit=%s\n", commit_label())
    sayf("results=%s\n", report_path)
    ensure_dataset(DATA_DIR; n_files=N_FILES, rows=ROWS)
    if WORKSPACE_ONLY
        workspace_cold_build(DATA_DIR)
        return nothing
    end
    project = bench_project()
    sample = joinpath(DATA_DIR, "meas_00001.csv")
    micro_benchmarks(sample)
    scan, source, records = macro_timings(project, DATA_DIR)
    expansion_timings(DATA_DIR)
    workspace_cold_build(DATA_DIR)
    region_timings(project, source, records)
    return nothing
end

function run_with_report(path::AbstractString)::Nothing
    open(path, "w") do file
        previous = OUT[]
        OUT[] = TeeIO(stdout, file)
        try
            main(path)
        finally
            flush(OUT[])
            OUT[] = previous
        end
    end
    println("saved benchmark report: ", path)
    return nothing
end

# Enable engine timings here at top level so the new method definitions are in a newer world age than
# the workload calls below (enabling inside the workload function would not take effect — world age).
MB.Profiling.enable!()
try
    run_with_report(result_path())
finally
    MB.Profiling.disable!()
end
