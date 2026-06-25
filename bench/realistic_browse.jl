# Realistic browse-while-building benchmark.
#
# Models the shape of a real project (the RuO2 v2 project): a few thousand source items across three
# kinds, with the same diversity — many tiny files, some medium files, and a handful of big files that
# each expand into many items (the fatigue-style one-file-many-cycles case that builds the big payload
# tables). It then does what a user does: while the scan/analysis is still running it *selects and
# reads* items (the expensive half of plotting), timing each read. That is the number that has to stay
# small for the app to feel responsive — the whole point of the cache.
#
#   julia --project=bench --threads=auto bench/realistic_browse.jl [scale]
#
# `scale` (default 1.0) multiplies all file/item counts; use it to hit ~1 min per build on your
# machine. Synthetic data and the cache live in a temp dir that is deleted on exit, so the drive is
# not filled; only the small result files (CSVs + PNG) are kept under bench/results/.
#
# Tunables via ENV (counts are per the documented diversity; see DEFAULTS below):
#   MB_BENCH_KIND1_FILES, MB_BENCH_KIND2_FILES, MB_BENCH_KIND3_FILES, MB_BENCH_KIND3_CYCLES
# Set MB_PROFILE_INTERNAL=1 for a full structured capture; add MB_PROFILE_CPU=1 for Julia sampling.
# Set MB_BENCH_START_PROFILE=0 to measure enabled-but-idle overhead.

using MeasurementBrowser
const MB = MeasurementBrowser
using CSV
using DBInterface
using DataFrames
using Random
using Printf
using Statistics: mean, median, quantile

# Optional: render a PNG when CairoMakie is in the bench env. Imported at top level so the plotting
# call below runs in a new-enough world age (importing inside the function fails with "method too new").
const HAS_CAIRO = try
    @eval import CairoMakie
    true
catch
    false
end

# --------------------------------------------------------------------------------------------------
# Sizing (downscaled from the real 3034-file / 20 GB project, same diversity)
# --------------------------------------------------------------------------------------------------

scale = length(ARGS) >= 1 ? parse(Float64, ARGS[1]) : 1.0
_env_int(key, default) = parse(Int, get(ENV, key, string(default)))
_scaled(n) = max(1, round(Int, n * scale))
const PROFILE_INTERNAL = MB.Profiling.environment_flag("MB_PROFILE_INTERNAL")
const PROFILE_CPU = MB.Profiling.environment_flag("MB_PROFILE_CPU")
const START_PROFILE = MB.Profiling.environment_flag(
    "MB_BENCH_START_PROFILE", PROFILE_INTERNAL)
START_PROFILE && !PROFILE_INTERNAL && error(
    "MB_BENCH_START_PROFILE=1 requires MB_PROFILE_INTERNAL=1",
)

# kind1: tiny files, one item each (IV-style, ~120 rows).      many → the long tail of small files
# kind2: medium files, one item each (CV-style, ~5000 rows).   some → the mid bucket
# kind3: big files, one item PER CYCLE (fatigue-style).        few files, many items, big payload
# Defaults tuned for ~1 min per build on an 8-thread laptop (~10k items, the fatigue-style kind3
# dominating the payload). Pass a scale arg or set the ENV counts to retune for your machine.
const KIND1_FILES  = _scaled(_env_int("MB_BENCH_KIND1_FILES", 3800))
const KIND2_FILES  = _scaled(_env_int("MB_BENCH_KIND2_FILES", 600))
const KIND3_FILES  = _scaled(_env_int("MB_BENCH_KIND3_FILES", 100))
const KIND3_CYCLES = _env_int("MB_BENCH_KIND3_CYCLES", 60)   # items per big file
const KIND1_ROWS   = 120
const KIND2_ROWS   = 5_000
const KIND3_ROWS   = 280     # rows per cycle

const MAX_BUILD_SECONDS = 600   # safety cap

# --------------------------------------------------------------------------------------------------
# Synthetic data generation (deterministic; written as raw CSV for speed)
# --------------------------------------------------------------------------------------------------

"""Write `n` rows of `cols` (name => generator(i)) as CSV to `path`, fast."""
function _write_csv(path::String, header::String, body::String)
    open(path, "w") do io
        write(io, header)
        write(io, body)
    end
    return nothing
end

function _kind1_body(rng)
    io = IOBuffer()
    for i in 1:KIND1_ROWS
        @printf(io, "%.4f,%.6e\n", -2 + 4i / KIND1_ROWS, (1e-9) * sinpi(i / 40) + 1e-11 * randn(rng))
    end
    return String(take!(io))
end

function _kind2_body(rng)
    io = IOBuffer()
    for i in 1:KIND2_ROWS
        @printf(io, "%.4f,%.6e,%.1f\n", -3 + 6i / KIND2_ROWS,
            1e-12 * (1 + cospi(i / 500)) + 1e-14 * randn(rng), 1.0e6)
    end
    return String(take!(io))
end

function _kind3_body(rng)
    io = IOBuffer()
    for c in 1:KIND3_CYCLES, r in 1:KIND3_ROWS
        @printf(io, "%d,%.6e,%.4f,%.6e\n", c, r * 1e-7,
            3 * sinpi(r / (KIND3_ROWS / 2)), (1e-6) * cospi(r / (KIND3_ROWS / 2)) + 1e-8 * randn(rng))
    end
    return String(take!(io))
end

"""Generate the whole synthetic tree under `root`; return (file_count, approx_item_count)."""
function generate_data(root::String)
    rng = MersenneTwister(20260624)
    mkpath(root)
    files = 0
    # Spread files across a few wafers/sites so the hierarchy has real structure.
    wafers = ["W$(i)" for i in 1:6]
    for n in 1:KIND1_FILES
        w = wafers[mod1(n, length(wafers))]
        dir = joinpath(root, w, "kind1"); mkpath(dir)
        _write_csv(joinpath(dir, @sprintf("dev%04d_kind1.csv", n)),
            "voltage,current\n", _kind1_body(rng)); files += 1
    end
    for n in 1:KIND2_FILES
        w = wafers[mod1(n, length(wafers))]
        dir = joinpath(root, w, "kind2"); mkpath(dir)
        _write_csv(joinpath(dir, @sprintf("dev%04d_kind2.csv", n)),
            "voltage,cap,freq\n", _kind2_body(rng)); files += 1
    end
    for n in 1:KIND3_FILES
        w = wafers[mod1(n, length(wafers))]
        dir = joinpath(root, w, "kind3"); mkpath(dir)
        _write_csv(joinpath(dir, @sprintf("dev%04d_kind3.csv", n)),
            "cycle,time,voltage,current\n", _kind3_body(rng)); files += 1
    end
    items = KIND1_FILES + KIND2_FILES + KIND3_FILES * KIND3_CYCLES
    return files, items
end

# --------------------------------------------------------------------------------------------------
# Project: three kinds mirroring the real read/entries/process/stats/plot shape
# --------------------------------------------------------------------------------------------------

_collection(file) = [splitpath(dirname(file.filepath))[end-1], splitpath(dirname(file.filepath))[end]]

function build_project()
    project = MB.define_project("BenchRealistic"; description="Realistic browse-while-build benchmark")

    MB.register_item!(project, :kind1;
        detect  = file -> endswith(file.filename, "_kind1.csv"),
        read    = file -> DataFrame(CSV.File(file.filepath; ntasks=1)),
        entries = (file, data) -> [MB.DataItem(; kind=:kind1, collection=_collection(file),
            label=file.filename, data=data, id=file.filepath * "#kind1")],
        process = item -> MB.DataItem(item, transform(item.data, [:voltage, :current] =>
            ByRow((v, i) -> iszero(v) ? 0.0 : i / v) => :conductance)),
        stats   = item -> Dict{Symbol,Any}(:imax => maximum(abs, item.data.current)),
        label   = item -> "K1 $(item.label)")

    MB.register_item!(project, :kind2;
        detect  = file -> endswith(file.filename, "_kind2.csv"),
        read    = file -> DataFrame(CSV.File(file.filepath; ntasks=1)),
        entries = (file, data) -> [MB.DataItem(; kind=:kind2, collection=_collection(file),
            label=file.filename, data=data, id=file.filepath * "#kind2")],
        stats   = item -> Dict{Symbol,Any}(:cmean => mean(item.data.cap)),
        label   = item -> "K2 $(item.label)")

    # The fatigue-style kind: one file → one item per cycle (where item count explodes).
    MB.register_item!(project, :kind3;
        detect  = file -> endswith(file.filename, "_kind3.csv"),
        read    = file -> DataFrame(CSV.File(file.filepath; ntasks=1)),
        entries = (file, data) -> [
            MB.DataItem(; kind=:kind3, collection=_collection(file),
                label=string(file.filename, " cycle ", c), parameters=Dict{Symbol,Any}(:cycle => c),
                data=data[data.cycle .== c, :], id=string(file.filepath, "#kind3,cycle=", c))
            for c in sort(unique(data.cycle))],
        process = item -> MB.DataItem(item, transform(item.data, [:voltage, :current] =>
            ByRow((v, i) -> v * i) => :power)),
        stats   = item -> Dict{Symbol,Any}(:pmax => maximum(abs, item.data.current)),
        label   = item -> "K3 $(item.label)")

    for k in (:kind1, :kind2, :kind3)
        MB.register_plot!(project, k; label="$k", setup=(ws, items) -> nothing,
            draw=(ws, items, fig) -> nothing)
    end
    return project
end

# --------------------------------------------------------------------------------------------------
# Driver: poll the workspace while timing interactive reads (selection + plot data load)
# --------------------------------------------------------------------------------------------------

build_idle(ws) =
    !MB.Workspace.source_scan_running(ws) &&
    !MB.Workspace.processing_work_running(ws) &&
    !MB.Workspace.analysis_work_running(ws) &&
    ws.analysis.state != :pending &&
    ws.scan.state in (:done, :unchanged, :error, :canceled)

"""Time one interactive read of `k` already-processed items of `kind`."""
function timed_read!(ws, kind::Symbol, k::Int)
    ready = String[]
    for id in keys(ws.index.item_stats)
        rec = get(ws.index.items, id, nothing)
        rec === nothing && continue
        rec.kind === kind && push!(ready, id)
        length(ready) >= 64 && break
    end
    length(ready) < k && return (nothing, 0, 0, length(ready))
    pick = ready[1:k]
    records = MB.ItemIndex.ItemRecord[ws.index.items[id] for id in pick]
    MB.select_items!(ws, records)                 # mirror the GUI selecting them
    timed = @timed MB.Workspace.materialize_items(ws, records)
    return (timed.time * 1e3, timed.bytes, k, length(ready))
end

mutable struct Sample
    elapsed_s::Float64
    phase::Symbol
    kind::Symbol
    n::Int
    read_ms::Float64
    allocated_bytes::Int
    ready::Int
end

mutable struct MemorySample
    elapsed_s::Float64
    rss_bytes::Int64
    gc_live_bytes::Int64
    gc_allocated_bytes::Int64
    rss_minus_gc_live_bytes::Int64
    index_items::Int64
    index_collections::Int64
    item_stats::Int64
    analysis_errors::Int64
    processing_jobs::Int64
    pending_writes::Int64
    pending_write_rows::Int64
    selected_queue::Int64
    background_waiting::Int64
end

function MemorySample(elapsed_s::Float64, snapshot)::MemorySample
    return MemorySample(
        elapsed_s,
        snapshot.rss_bytes,
        snapshot.gc_live_bytes,
        snapshot.gc_allocated_bytes,
        snapshot.rss_minus_gc_live_bytes,
        snapshot.index_items,
        snapshot.index_collections,
        snapshot.item_stats,
        snapshot.analysis_errors,
        snapshot.processing_jobs,
        snapshot.pending_writes,
        snapshot.pending_write_rows,
        snapshot.selected_queue,
        snapshot.background_waiting,
    )
end

"""Measure one full-table aggregate for every processed payload schema."""
function global_query_samples(ws)
    return MB.Cache.with_reader(ws.cache.db) do connection
        tables = collect(DBInterface.execute(connection, """
            SELECT min(i.kind) AS kind, d.storage_id
            FROM item_data d
            JOIN items i ON i.id = d.item_id
            WHERE d.stage = 'processed'
            GROUP BY d.storage_id
            ORDER BY kind
        """))
        results = NamedTuple[]
        for row in tables
            kind = Symbol(String(row.kind))
            table = "\"dataframe_$(String(row.storage_id))\""
            times = Float64[]
            for _ in 1:6
                timed = @timed only(DBInterface.execute(
                    connection, "SELECT sum(c1), avg(c2) FROM $table"))
                push!(times, timed.time * 1e3)
            end
            push!(results, (kind=kind, median_ms=median(times[2:end])))
        end
        return results
    end
end

function run_benchmark()
    tmp = mktempdir()
    pushfirst!(DEPOT_PATH, tmp)          # cache lands in temp, deleted with everything else
    data_root = joinpath(tmp, "data")

    outdir = joinpath(@__DIR__, "results",
        "realistic-" * replace(string(round(Int, time())), r"\D" => ""))
    mkpath(outdir)

    println("Generating synthetic data … (scale=$scale)")
    gen_t = @elapsed (n_files, n_items) = generate_data(data_root)
    data_bytes = sum(filesize(joinpath(r, f))
                     for (r, _, fs) in walkdir(data_root) for f in fs)
    @printf("  %d files, ~%d items, %.1f MB on disk, generated in %.1fs\n",
        n_files, n_items, data_bytes / 1024^2, gen_t)

    project = build_project()
    kinds = (:kind1, :kind2, :kind3)
    samples = Sample[]

    println("Building cache + browsing during the scan …")
    profile_output = PROFILE_INTERNAL ? something(
        MB.Profiling.environment_path("MB_PROFILE_OUTPUT"),
        joinpath(outdir, "profile.json"),
    ) : nothing
    ws = MB.open_workspace(
        project,
        data_root;
        profile_internal=PROFILE_INTERNAL,
        profile_cpu=PROFILE_CPU,
        profile_output=profile_output,
    )
    if START_PROFILE
        MB.Workspace.start_internal_profile!(ws)
        deadline = time() + 60
        while ws.profiler.state !== :recording
            MB.Workspace.poll_workspace!(ws)
            time() < deadline || error(
                "Profiler did not reach recording state within 60 seconds; state=$(ws.profiler.state)",
            )
            sleep(0.004)
        end
    end
    t_start = time()
    rss_start_bytes = MB.Profiling.process_rss_bytes()
    rss_peak_bytes = rss_start_bytes
    rss_end_bytes = rss_start_bytes
    build_seconds = 0.0
    scan_seconds = 0.0
    processing_started = nothing
    processing_seconds = 0.0
    analysis_started = nothing
    analysis_seconds = 0.0
    global_queries = NamedTuple[]
    build_stats = nothing
    profile_report = nothing
    try
        last_probe = 0.0
        last_rss_sample = 0.0
        last_memory_sample = -Inf
        memory_samples = MemorySample[]
        kind_cursor = 1
        while true
            MB.Workspace.poll_workspace!(ws)
            now = time() - t_start
            if now - last_rss_sample >= 0.1
                last_rss_sample = now
                rss_peak_bytes = max(rss_peak_bytes, MB.Profiling.process_rss_bytes())
            end
            if now - last_memory_sample >= 0.5
                last_memory_sample = now
                snapshot = MB.Workspace.workspace_memory_snapshot(ws)
                rss_peak_bytes = max(rss_peak_bytes, snapshot.rss_bytes)
                push!(memory_samples, MemorySample(now, snapshot))
            end
            completed, queued = lock(ws.processing.lock) do
                (ws.processing.completed, ws.processing.total)
            end
            processing_started === nothing && queued > 0 && (processing_started = now)
            scan_seconds == 0 && !MB.Workspace.source_scan_running(ws) &&
                ws.scan.state in (:done, :unchanged, :error, :canceled) &&
                (scan_seconds = now)
            if processing_started !== nothing && processing_seconds == 0 && scan_seconds > 0 &&
               completed == queued && !MB.Workspace.processing_work_running(ws)
                processing_seconds = now - processing_started
                analysis_started === nothing && (analysis_started = now)
            end
            analysis_started === nothing && ws.analysis.state == :analyzing &&
                (analysis_started = now)
            analysis_started !== nothing && analysis_seconds == 0 &&
                ws.analysis.state in (:done, :error, :canceled) &&
                (analysis_seconds = now - analysis_started)
            # Probe responsiveness ~6×/s, rotating across kinds, once items exist.
            if now - last_probe >= 0.16
                last_probe = now
                kind = kinds[kind_cursor]; kind_cursor = mod1(kind_cursor + 1, length(kinds))
                ms, bytes, n, ready = timed_read!(ws, kind, 3)
                ms === nothing || push!(samples,
                    Sample(now, :during_build, kind, n, ms, bytes, ready))
            end
            if build_idle(ws)
                build_seconds = now
                break
            end
            (now > MAX_BUILD_SECONDS) && (build_seconds = now;
                @warn("hit MAX_BUILD_SECONDS"); break)
            sleep(0.004)
        end
        MB.Workspace.poll_workspace!(ws)
        final_memory = MB.Workspace.workspace_memory_snapshot(ws)
        rss_end_bytes = final_memory.rss_bytes
        push!(memory_samples, MemorySample(time() - t_start, final_memory))

        # Steady-state sweep: many random selections per kind on the finished cache.
        rng = MersenneTwister(1)
        for kind in kinds, _ in 1:40
            ids = [id for id in keys(ws.index.item_stats)
                   if (r = get(ws.index.items, id, nothing); r !== nothing && r.kind === kind)]
            isempty(ids) && continue
            k = rand(rng, 1:min(4, length(ids)))
            records = MB.ItemIndex.ItemRecord[ws.index.items[id] for id in rand(rng, ids, k)]
            timed = @timed MB.Workspace.materialize_items(ws, records)
            push!(samples, Sample(time() - t_start, :after_build, kind, k,
                timed.time * 1e3, timed.bytes, length(ids)))
        end
        global_queries = global_query_samples(ws)
        completed, queued = lock(ws.processing.lock) do
            (ws.processing.completed, ws.processing.total)
        end
        metrics = ws.metrics
        db = ws.cache.db
        build_stats = (
            scan_seconds,
            processing_seconds,
            analysis_seconds,
            processed_items=length(ws.index.items),
            completed_jobs=completed,
            queued_items=queued,
            collection_nodes=ws.analysis.progress.total_source_items,
            interpreted_write_ns=metrics.interpreted_write_ns[],
            interpreted_writes=metrics.interpreted_writes[],
            processed_write_ns=metrics.processed_write_ns[],
            processed_writes=metrics.processed_writes[],
            stats_write_ns=metrics.stats_write_ns[],
            stats_writes=metrics.stats_writes[],
            writer_busy_ns=db.writer_busy_ns[],
            writer_wait_ns=db.writer_wait_ns[],
            rss_start_bytes,
            rss_peak_bytes,
            rss_end_bytes,
            memory_samples,
        )
        START_PROFILE &&
            (profile_report = MB.Workspace.stop_internal_profile!(ws))
    finally
        MB.close_workspace!(ws)
    end

    report(samples, global_queries, build_stats, profile_report, outdir,
        n_files, n_items, data_bytes, build_seconds)

    rm(tmp; force=true, recursive=true)   # synthetic data + cache gone; results kept
    println("\nResults kept in: $outdir")
    return outdir
end

# --------------------------------------------------------------------------------------------------
# Reporting
# --------------------------------------------------------------------------------------------------

_stat(v, f) = isempty(v) ? NaN : f(v)

function report(samples, global_queries, stats, profile_report, outdir,
    n_files, n_items, data_bytes, build_seconds)
    # responsiveness CSV
    open(joinpath(outdir, "responsiveness.csv"), "w") do io
        println(io, "elapsed_s,phase,kind,n_items,read_ms,allocated_bytes,ready_items")
        for s in samples
            @printf(io, "%.3f,%s,%s,%d,%.3f,%d,%d\n",
                s.elapsed_s, s.phase, s.kind, s.n, s.read_ms, s.allocated_bytes, s.ready)
        end
    end

    open(joinpath(outdir, "memory_samples.csv"), "w") do io
        println(io, "elapsed_s,rss_bytes,gc_live_bytes,gc_allocated_bytes," *
            "rss_minus_gc_live_bytes,index_items,index_collections,item_stats,analysis_errors," *
            "processing_jobs,pending_writes,pending_write_rows,selected_queue,background_waiting")
        for sample in stats.memory_samples
            @printf(io, "%.3f,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d\n",
                sample.elapsed_s,
                sample.rss_bytes,
                sample.gc_live_bytes,
                sample.gc_allocated_bytes,
                sample.rss_minus_gc_live_bytes,
                sample.index_items,
                sample.index_collections,
                sample.item_stats,
                sample.analysis_errors,
                sample.processing_jobs,
                sample.pending_writes,
                sample.pending_write_rows,
                sample.selected_queue,
                sample.background_waiting)
        end
    end

    open(joinpath(outdir, "global_queries.csv"), "w") do io
        println(io, "kind,median_ms")
        for sample in global_queries
            @printf(io, "%s,%.3f\n", sample.kind, sample.median_ms)
        end
    end

    if profile_report isa MB.Profiling.ProfileReport
        open(joinpath(outdir, "profile_summary.csv"), "w") do io
            println(io, "category,operation,count,total_ms,p50_ms,p90_ms,p99_ms,max_ms," *
                "wait_ms,service_ms,mean_batch,p50_batch,p90_batch,max_batch")
            for row in profile_report.summary
                @printf(io, "%s,%s,%d,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f\n",
                    row.category, row.operation, row.count, row.total_ms, row.median_ms,
                    row.p90_ms, row.p99_ms, row.max_ms, row.wait_ms, row.service_ms,
                    row.mean_batch, row.median_batch, row.p90_batch, row.max_batch)
            end
        end
        if profile_report.cpu !== nothing
            open(joinpath(outdir, "cpu_hotspots.csv"), "w") do io
                println(io, "samples,self_samples,function,file,line")
                for row in profile_report.cpu.rows
                    function_name = replace(row.function_name, '"' => "\"\"")
                    file = replace(row.file, '"' => "\"\"")
                    @printf(io, "%d,%d,\"%s\",\"%s\",%d\n",
                        row.samples, row.self_samples, function_name, file, row.line)
                end
            end
        end
    end

    write_calls = stats.interpreted_writes + stats.processed_writes + stats.stats_writes
    write_ns = stats.interpreted_write_ns + stats.processed_write_ns + stats.stats_write_ns
    mean_write_ms = write_calls == 0 ? 0.0 : write_ns / write_calls / 1e6
    during = [s.read_ms for s in samples if s.phase === :during_build]
    after = [s.read_ms for s in samples if s.phase === :after_build]
    read_stat(values, statistic) = isempty(values) ? NaN : statistic(values)
    open(joinpath(outdir, "scorecard.csv"), "w") do io
        println(io, "build_s,scan_s,scan_per_s,processing_s,items_per_s,analysis_s," *
            "end_to_end_items_per_s,mean_write_ms,writer_busy_s,writer_wait_s," *
            "during_median_ms,during_p90_ms,during_p99_ms,during_max_ms," *
            "after_median_ms,after_p90_ms,after_p99_ms,after_max_ms," *
            "rss_start_mib,rss_peak_mib,rss_end_mib,profile_events,profile_dropped")
        event_count = profile_report isa MB.Profiling.ProfileReport ?
            length(profile_report.events) : 0
        dropped = profile_report isa MB.Profiling.ProfileReport ?
            profile_report.dropped_events : 0
        @printf(io, "%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.1f,%.1f,%.1f,%d,%d\n",
            build_seconds, stats.scan_seconds, n_files / stats.scan_seconds,
            stats.processing_seconds, stats.processed_items / stats.processing_seconds,
            stats.analysis_seconds, stats.processed_items / build_seconds, mean_write_ms,
            stats.writer_busy_ns / 1e9, stats.writer_wait_ns / 1e9,
            read_stat(during, median), read_stat(during, values -> quantile(values, 0.9)),
            read_stat(during, values -> quantile(values, 0.99)), read_stat(during, maximum),
            read_stat(after, median), read_stat(after, values -> quantile(values, 0.9)),
            read_stat(after, values -> quantile(values, 0.99)), read_stat(after, maximum),
            stats.rss_start_bytes / 1024^2,
            stats.rss_peak_bytes / 1024^2,
            stats.rss_end_bytes / 1024^2,
            event_count, dropped)
    end

    println("\n==================== REALISTIC BROWSE BENCHMARK ====================")
    @printf("dataset:  %d files · ~%d items · %.1f MB\n", n_files, n_items, data_bytes / 1024^2)
    @printf("build:    %.1f s wall (scan + processing + collection analysis)\n", build_seconds)
    n_during = count(s -> s.phase === :during_build, samples)
    @printf("probes:   %d during build · %d after\n", n_during,
        count(s -> s.phase === :after_build, samples))

    println("\nThroughput:")
    @printf("  scan                %8.1f source items/s  (%6.1f s)\n",
        n_files / stats.scan_seconds, stats.scan_seconds)
    @printf("  item processing     %8.1f items/s         (%6.1f s, %d unique items)\n",
        stats.processed_items / stats.processing_seconds,
        stats.processing_seconds, stats.processed_items)
    stats.completed_jobs == stats.processed_items || @printf(
        "  duplicate queue work %8d cache-hit jobs\n",
        stats.completed_jobs - stats.processed_items,
    )
    @printf("  collection analysis %8.1f nodes/s         (%6.1f s, %d nodes)\n",
        stats.collection_nodes / stats.analysis_seconds,
        stats.analysis_seconds, stats.collection_nodes)
    @printf("  end to end          %8.1f items/s\n", stats.processed_items / build_seconds)

    println("\nWrites:")
    @printf("  interpreted %6d calls  mean %7.2f ms\n", stats.interpreted_writes,
        stats.interpreted_write_ns / max(stats.interpreted_writes, 1) / 1e6)
    @printf("  processed   %6d calls  mean %7.2f ms  mean batch %5.1f items\n",
        stats.processed_writes,
        stats.processed_write_ns / max(stats.processed_writes, 1) / 1e6,
        stats.processed_items / max(stats.processed_writes, 1))
    @printf("  stats       %6d calls  mean %7.2f ms\n", stats.stats_writes,
        stats.stats_write_ns / max(stats.stats_writes, 1) / 1e6)
    @printf("  overall mean %7.2f ms  writer busy %.1f s  queued wait %.1f s\n",
        mean_write_ms, stats.writer_busy_ns / 1e9, stats.writer_wait_ns / 1e9)

    println("\nProcess memory:")
    @printf("  RSS start %.1f MiB  peak %.1f MiB  end %.1f MiB\n",
        stats.rss_start_bytes / 1024^2,
        stats.rss_peak_bytes / 1024^2,
        stats.rss_end_bytes / 1024^2)
    if !isempty(stats.memory_samples)
        peak_sample = stats.memory_samples[argmax(
            [sample.rss_bytes for sample in stats.memory_samples])]
        @printf("  peak sample: GC live %.1f MiB  RSS-GC-live %.1f MiB\n",
            peak_sample.gc_live_bytes / 1024^2,
            peak_sample.rss_minus_gc_live_bytes / 1024^2)
        @printf("  index counts: items %d  stats %d  collections %d  errors %d\n",
            peak_sample.index_items,
            peak_sample.item_stats,
            peak_sample.index_collections,
            peak_sample.analysis_errors)
        @printf("  queue counts: jobs %d  pending writes %d (%d rows)  selected %d  background waiting %d\n",
            peak_sample.processing_jobs,
            peak_sample.pending_writes,
            peak_sample.pending_write_rows,
            peak_sample.selected_queue,
            peak_sample.background_waiting)
    end

    println("\nAggregate interactive read latency, ms:")
    @printf("  %-13s %8s %8s %8s %8s %6s\n",
        "phase", "median", "p90", "p99", "max", "n")
    for phase in (:during_build, :after_build)
        values = [s.read_ms for s in samples if s.phase === phase]
        isempty(values) && continue
        @printf("  %-13s %8.1f %8.1f %8.1f %8.1f %6d\n",
            phase, median(values), quantile(values, 0.9), quantile(values, 0.99),
            maximum(values), length(values))
    end

    println("\nRead latency by data shape, ms:")
    @printf("  %-13s %-8s %8s %8s %8s %8s %6s\n",
        "phase", "kind", "median", "p90", "p99", "max", "n")
    for phase in (:during_build, :after_build), kind in (:kind1, :kind2, :kind3)
        v = [s.read_ms for s in samples if s.phase === phase && s.kind === kind]
        isempty(v) && continue
        @printf("  %-13s %-8s %8.1f %8.1f %8.1f %8.1f %6d\n",
            phase, kind, median(v), quantile(v, 0.9), quantile(v, 0.99),
            maximum(v), length(v))
    end

    println("\nAllocated per materialization, MiB:")
    for phase in (:during_build, :after_build), kind in (:kind1, :kind2, :kind3)
        bytes = [s.allocated_bytes for s in samples if s.phase === phase && s.kind === kind]
        isempty(bytes) || @printf("  %-13s %-8s median %6.2f  p90 %6.2f\n",
            phase, kind, median(bytes) / 1024^2, quantile(bytes, 0.9) / 1024^2)
    end

    println("\nGlobal aggregate latency, ms:")
    for sample in global_queries
        @printf("  %-8s %8.1f\n", sample.kind, sample.median_ms)
    end

    if profile_report isa MB.Profiling.ProfileReport
        @printf("\nStructured profile: %d events, %d counters, %d dropped, %d CPU samples\n",
            length(profile_report.events), length(profile_report.counters),
            profile_report.dropped_events,
            profile_report.cpu === nothing ? 0 : profile_report.cpu.total_samples)
    end

    allduring = [s.read_ms for s in samples if s.phase === :during_build]
    if !isempty(allduring)
        slow = count(>(200), allduring)
        @printf("\n  during-build reads over 200 ms: %d / %d (%.1f%%)  ·  worst %.0f ms\n",
            slow, length(allduring), 100 * slow / length(allduring), maximum(allduring))
    end

    _maybe_plot(samples, outdir)
    return nothing
end

"""Render a PNG if CairoMakie is available in the bench env; otherwise point at the CSVs."""
function _maybe_plot(samples, outdir)
    if !HAS_CAIRO
        println("\n(No CairoMakie in the bench env — skipping PNG.)")
        return nothing
    end
    CM = CairoMakie
    fig = CM.Figure(size=(1100, 750), fontsize=14)
    CM.Label(fig[0, 1:2], "Browse-while-build responsiveness", fontsize=19, font=:bold)
    colors = Dict(:kind1 => :steelblue, :kind2 => :darkorange, :kind3 => :crimson)

    ax1 = CM.Axis(fig[1, 1:2]; xlabel="build elapsed (s)", ylabel="read latency (ms)",
        title="Interactive read latency over the build (each point = one selection)")
    for kind in (:kind1, :kind2, :kind3)
        pts = [(s.elapsed_s, s.read_ms) for s in samples
               if s.phase === :during_build && s.kind === kind]
        isempty(pts) && continue
        CM.scatter!(ax1, first.(pts), last.(pts); color=colors[kind], markersize=7, label=string(kind))
    end
    CM.axislegend(ax1; position=:lt)

    ax2 = CM.Axis(fig[2, 1]; xlabel="read latency (ms)", ylabel="count",
        title="During build")
    ax3 = CM.Axis(fig[2, 2]; xlabel="read latency (ms)", ylabel="count",
        title="After build")
    for kind in (:kind1, :kind2, :kind3)
        d = [s.read_ms for s in samples if s.phase === :during_build && s.kind === kind]
        a = [s.read_ms for s in samples if s.phase === :after_build && s.kind === kind]
        isempty(d) || CM.hist!(ax2, d; bins=20, color=(colors[kind], 0.5), label=string(kind))
        isempty(a) || CM.hist!(ax3, a; bins=20, color=(colors[kind], 0.5), label=string(kind))
    end
    CM.axislegend(ax2; position=:rt)

    path = joinpath(outdir, "responsiveness.png")
    CM.save(path, fig)
    println("\nplot: $path")
    return nothing
end

run_benchmark()
