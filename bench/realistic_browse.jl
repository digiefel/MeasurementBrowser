# Realistic browse-while-building benchmark.
#
# Models the pressure profile of a real project (the RuO2 v2 project): fewer source files than the
# full data set, but realistic fatigue-style row volume. The important scale is staged DataFrame rows:
# the real run crosses the cache buffer row ceiling and forces the writer/backpressure path.
#
# Source bytes are three static CSVs under bench/templates/. BenchSource presents each file as many
# source items (scale multiplies alias counts, not rows per file).
#
# Peak RAM: bench/run.sh realistic_browse.jl [scale]
#
#   julia --project=bench --threads=auto bench/realistic_browse.jl [scale]
#
# `scale` (default 1.0) multiplies alias counts. The DuckDB cache lives in a temp depot deleted on
# exit; only the small result files are kept under bench/results/.
#
# Tunables via ENV (alias counts; see DEFAULTS below):
#   MB_BENCH_KIND1_FILES, MB_BENCH_KIND2_FILES, MB_BENCH_KIND3_FILES
#   MB_BENCH_AFTER_BUILD_PLOTS, MB_BENCH_PROCESSED_STRESS_ROWS, MB_BENCH_REPEATS
# Loading DataBrowserProfiling turns on `@timed_dbg`; the run then `reset_debug_timings!` /
# `take_debug_timings!` and writes the TimerOutput.

using DataBrowserAPI: item_data, label
using DataBrowserRecipes: define_project, register_item!
using DataBrowserAPI.ItemIndex: ItemRecord, collection_item_ids
import DataBrowserCore.Workspace as Workspace
using DataBrowserCore.Workspace: close_workspace!, open_workspace, select_items!
using DataBrowserCache
using DataBrowserPlots:
    register_plot!,
    registered_plot_kinds,
    setup_plot,
    plot_data!
import DataBrowserProfiling as Profiling
using CSV
using DataFrames
using Dates
using Random
using Printf
using Statistics: mean, median, quantile
import GLMakie: Figure, Axis, lines!, contents

include(joinpath(@__DIR__, "custom_data_source.jl"))

# --------------------------------------------------------------------------------------------------
# Sizing (alias counts; row layout is the templates)
# --------------------------------------------------------------------------------------------------

scale = length(ARGS) >= 1 ? parse(Float64, ARGS[1]) : 1.0
_env_int(key, default) = parse(Int, get(ENV, key, string(default)))
_scaled(n) = max(1, round(Int, n * scale))

# kind1: tiny files, one item each (IV-style).      many small records
# kind2: medium files, one item each (CV-style).    mid-size table writes
# kind3: big files, one item PER CYCLE (fatigue).   few aliases, many rows/items
const KIND1_FILES = _scaled(_env_int("MB_BENCH_KIND1_FILES", 500))
const KIND2_FILES = _scaled(_env_int("MB_BENCH_KIND2_FILES", 120))
const KIND3_FILES = _scaled(_env_int("MB_BENCH_KIND3_FILES", 16))
# Rows / cycles in bench/templates/*.csv (kind3: one item per cycle).
const KIND1_ROWS = 120
const KIND2_ROWS = 5_000
const KIND3_CYCLES = 96
const KIND3_ROWS = 6_000

const CACHE_ROW_CEILING = DataBrowserCache.CACHE_BUFFER_ROW_LIMIT
const PROCESSED_STRESS_ROWS = _env_int(
    "MB_BENCH_PROCESSED_STRESS_ROWS",
    2 * CACHE_ROW_CEILING,
)
const AFTER_BUILD_PLOTS = _env_int("MB_BENCH_AFTER_BUILD_PLOTS", 40)
const BENCH_REPEATS = _env_int("MB_BENCH_REPEATS", 1)
const MAX_BUILD_SECONDS = 600   # safety cap
const REQUIRE_SATURATION = scale >= 1.0
const ESTIMATED_PAYLOAD_ROWS = Int64(KIND1_FILES) * KIND1_ROWS +
                               Int64(KIND2_FILES) * KIND2_ROWS +
                               Int64(KIND3_FILES) * KIND3_CYCLES * KIND3_ROWS

const RUN_LOG = Ref{Union{Nothing,IO}}(nothing)
const BENCH_ENV_KEYS = (
    "MB_BENCH_KIND1_FILES",
    "MB_BENCH_KIND2_FILES",
    "MB_BENCH_KIND3_FILES",
    "MB_BENCH_AFTER_BUILD_PLOTS",
    "MB_BENCH_PROCESSED_STRESS_ROWS",
    "MB_BENCH_REPEATS",
    "MB_BENCH_OUTDIR",
)

function tee_println(args...)::Nothing
    println(stdout, args...)
    io = RUN_LOG[]
    io === nothing || println(io, args...)
    return nothing
end

function tee_printf(format::AbstractString, args...)::Nothing
    Printf.format(stdout, Printf.Format(format), args...)
    io = RUN_LOG[]
    io === nothing || Printf.format(io, Printf.Format(format), args...)
    return nothing
end

function _repo_command(args::Vector{String})::String
    try
        return strip(read(Cmd(Cmd(args); dir=joinpath(@__DIR__, "..")), String))
    catch error
        return "unavailable ($(typeof(error)))"
    end
end

function _print_run_header(log_path::String, outdir::String)::Nothing
    tee_println("DataBrowser realistic benchmark")
    tee_println("started_at: ", Dates.format(now(), dateformat"yyyy-mm-dd HH:MM:SS"))
    tee_println("output_dir: ", outdir)
    tee_println("log_file:   ", log_path)
    tee_println("branch:     ", _repo_command(["git", "rev-parse", "--abbrev-ref", "HEAD"]))
    tee_println("commit:     ", _repo_command(["git", "rev-parse", "HEAD"]))
    status = _repo_command(["git", "status", "--short", "--",
        "bench/realistic_browse.jl", "bench/custom_data_source.jl", "bench/README.md",
        "bench/run.sh"])
    tee_println("benchmark_file_status:")
    if isempty(status)
        tee_println("  <clean>")
    else
        for line in split(status, '\n')
            tee_println("  ", line)
        end
    end
    tee_println("julia:      ", string(VERSION))
    tee_println("threads:    ", string(Base.Threads.nthreads()))
    tee_println("args:       ", isempty(ARGS) ? "<none>" : join(ARGS, " "))
    tee_println("environment:")
    for key in BENCH_ENV_KEYS
        tee_println("  ", key, "=", get(ENV, key, "<unset>"))
    end
    tee_println()
    return nothing
end

_template(name) = joinpath(@__DIR__, "templates", name)

function _aliases(name, n, kind)
    alias_file(_template(name), n, i -> joinpath(
        "W$(mod1(i, 6))", kind, string("dev", lpad(i, 4, '0'), "_", kind, ".csv")))
end

_bench_source() = BenchSource(abspath(_template("")), vcat(
    _aliases("kind1.csv", KIND1_FILES, "kind1"),
    _aliases("kind2.csv", KIND2_FILES, "kind2"),
    _aliases("kind3.csv", KIND3_FILES, "kind3"),
))

# --------------------------------------------------------------------------------------------------
# Project: three kinds mirroring the real read/entries/process/analyze/plot shape
# --------------------------------------------------------------------------------------------------

# Load the shared template at `filepath`. Wafer and kind come from the alias path
# (`W2/kind1/dev0042_kind1.csv`), not from the template's real folder.
function _read_table(file)
    parts = splitpath(dirname(file.relative_path))
    return (
        data=DataFrame(CSV.File(file.filepath; ntasks=1)),
        metadata=Dict{Symbol,Any}(
            :wafer => parts[end-1],
            :measurement_kind => parts[end],
        ),
    )
end

_collection(_data, metadata) =
    String[metadata[:wafer], metadata[:measurement_kind]]

function build_project(; plots::Bool=true)
    project = define_project("BenchRealistic"; description="Realistic browse-while-build benchmark")

    register_item!(project, :kind1;
        detect  = file -> endswith(file.filename, "_kind1.csv"),
        read    = _read_table,
        collection = _collection,
        process = (data, _metadata) -> transform(data, [:voltage, :current] =>
            ByRow((v, i) -> iszero(v) ? 0.0 : i / v) => :conductance),
        analyze = (data, _metadata) -> Dict{Symbol,Any}(:imax => maximum(abs, data.current)),
        label   = (_data, metadata) -> "K1 $(metadata[:filename])")

    register_item!(project, :kind2;
        detect  = file -> endswith(file.filename, "_kind2.csv"),
        read    = _read_table,
        collection = _collection,
        analyze = (data, _metadata) -> Dict{Symbol,Any}(:cmean => mean(data.cap)),
        label   = (_data, metadata) -> "K2 $(metadata[:filename])")

    # The fatigue-style kind: one file → one item per cycle (where item count explodes).
    register_item!(project, :kind3;
        detect  = file -> endswith(file.filename, "_kind3.csv"),
        read    = _read_table,
        entries = (data, _metadata) -> [
            (data=(@view data[data.cycle .== cycle, :]),
                metadata=Dict{Symbol,Any}(:cycle => cycle))
            for cycle in sort(unique(data.cycle))],
        collection = _collection,
        id = (_data, metadata) -> metadata[:cycle],
        process = (data, _metadata) -> transform(data, [:voltage, :current] =>
            ByRow((v, i) -> v * i) => :power),
        analyze = (data, _metadata) -> Dict{Symbol,Any}(:pmax => maximum(abs, data.current)),
        label   = (_data, metadata) ->
            "K3 $(metadata[:filename]) cycle $(metadata[:cycle])")

    plots || return project

    # Real plot recipes: one axis, one line per selected item, reading the processed columns — the
    # same scene-graph work the GUI does when a user selects items and the plot panel renders.
    _axis(fig, xlabel, ylabel) = Axis(fig[1, 1]; xlabel, ylabel)
    register_plot!(project, :kind1; label="IV",
        setup=(ws, items) -> (fig = Figure(); _axis(fig, "voltage", "conductance"); fig),
        draw=(ws, items, fig) -> for item in items
            d = item_data(item)
            lines!(contents(fig[1, 1])[1], d.voltage, d.conductance)
        end)
    register_plot!(project, :kind2; label="CV",
        setup=(ws, items) -> (fig = Figure(); _axis(fig, "voltage", "cap"); fig),
        draw=(ws, items, fig) -> for item in items
            d = item_data(item)
            lines!(contents(fig[1, 1])[1], d.voltage, d.cap)
        end)
    register_plot!(project, :kind3; label="power",
        setup=(ws, items) -> (fig = Figure(); _axis(fig, "time", "power"); fig),
        draw=(ws, items, fig) -> for item in items
            d = item_data(item)
            lines!(contents(fig[1, 1])[1], d.time, d.power)
        end)
    return project
end

# --------------------------------------------------------------------------------------------------
# Driver: poll the workspace while timing interactive reads (selection + plot data load)
# --------------------------------------------------------------------------------------------------

build_idle(ws) =
    !Workspace.workspace_busy(ws) &&
        ws.scan.state in (:done, :unchanged, :error, :canceled)

function _active_work_count(ws, kinds::Tuple)::Int
    return lock(ws.work.lock) do
        count(
            node -> node.key.kind in kinds && node.state in (:queued, :running),
            values(ws.work.nodes),
        )
    end
end

_processing_active(ws)::Bool =
    _active_work_count(ws, (Workspace.ITEM_PROCESS, Workspace.ITEM_ANALYZE)) > 0 ||
    Workspace.cache_has_pending_writes(ws.cache.db)

_analysis_active(ws)::Bool =
    _active_work_count(ws, (Workspace.COLLECTION_ANALYZE,)) > 0

"""Collect up to 64 item ids of `kind`."""
function _ready_ids(ws, kind::Symbol)
    ready = String[]
    for (id, rec) in ws.index.items
        label(rec.type) === kind && push!(ready, id)
        length(ready) >= 64 && break
    end
    return ready
end

"""
Time one plot probe on `k` items of `kind`: select, materialize, setup, and draw.
`records` may be supplied to plot a fixed selection (used for reopen probes); otherwise the first
`k` ready items are used.
"""
function timed_plot!(ws, plot_kinds, kind::Symbol, k::Int; records=nothing)
    if records === nothing
        ready = _ready_ids(ws, kind)
        length(ready) < k && return nothing
        records = ItemRecord[ws.index.items[id] for id in ready[1:k]]
    end
    n_ready = records === nothing ? 0 : length(records)
    result = @timed begin
        select_items!(ws, records)             # mirror the GUI selecting them
        items = Workspace.materialize_items(ws, records)
        plot_kind = plot_kinds[kind]
        figure = setup_plot(ws, plot_kind, items)
        plot_data!(ws, plot_kind, items, figure)
    end
    return (plot_ms=result.time * 1e3, bytes=result.bytes, n=length(records), ready=n_ready)
end

mutable struct Sample
    elapsed_s::Float64
    phase::Symbol
    kind::Symbol
    n::Int
    plot_ms::Float64       # secondary plot probe
    allocated_bytes::Int
    ready::Int
end

struct SaturationSample
    kind::Symbol
    requested_items::Int
    materialized_items::Int
    estimated_rows::Int64
    load_ms::Float64
    flush_ms::Float64
    peak_pending_rows::Int64
    processed_writes::Int
end

function _records_of_kind(ws, kind::Symbol)::Vector{ItemRecord}
    records = ItemRecord[
        record for record in values(ws.index.items)
        if label(record.type) === kind
    ]
    sort!(records; by=record -> record.id)
    return records
end

"""
Materialize one large selected batch and report whether that created processed writes.

If background work already wrote those processed payloads, this is a cache-read/materialization
probe, not a writer-saturation probe; `processed_writes` makes that visible.
"""
function saturate_processed_writes!(ws, kind::Symbol)::SaturationSample
    records = _records_of_kind(ws, kind)
    stress_items = ceil(Int, PROCESSED_STRESS_ROWS / KIND3_ROWS)
    stress_items = min(length(records), stress_items)
    if stress_items == 0
        @warn "Skipping processed-writer saturation: no completed $kind items exist"
        return SaturationSample(kind, 0, 0, 0, 0.0, 0.0, 0, 0)
    end
    selected = records[end-stress_items+1:end]
    estimated_rows = Int64(stress_items) * Int64(KIND3_ROWS)
    if REQUIRE_SATURATION && estimated_rows < CACHE_ROW_CEILING
        error(
            "Processed-writer stress is too small: selected $stress_items $kind item(s) " *
            "for about $estimated_rows rows, below the cache row ceiling $CACHE_ROW_CEILING. " *
            "Increase MB_BENCH_KIND3_FILES, or rewrite the kind3 template with more rows.",
        )
    end

    select_items!(ws, selected)
    processed_writes_before = ws.metrics.processed_writes[]
    load = @timed Workspace.materialize_items(ws, selected)
    peak_pending_rows = Int64(0)
    flush_started = time()
    while Workspace.cache_has_pending_writes(ws.cache.db)
        counts = Workspace.cache_pending_counts(ws.cache.db)
        peak_pending_rows = max(peak_pending_rows, Int64(counts.rows))
        sleep(0.004)
    end
    flush_ms = (time() - flush_started) * 1e3
    return SaturationSample(
        kind,
        stress_items,
        length(load.value),
        estimated_rows,
        load.time * 1e3,
        flush_ms,
        peak_pending_rows,
        ws.metrics.processed_writes[] - processed_writes_before,
    )
end

"""Write a taken `TimerOutput` as a readable tree and a Tables.jl CSV."""
function write_debug_timings(outdir, timings)
    open(joinpath(outdir, "debug_timings.txt"), "w") do io
        show(IOContext(io, :displaysize => (typemax(Int), typemax(Int))),
            MIME("text/plain"), timings)
    end
    CSV.write(joinpath(outdir, "debug_timings.csv"), timings)
    return nothing
end

function run_benchmark()
    tmp = mktempdir()
    pushfirst!(DEPOT_PATH, tmp)          # cache lands in temp, deleted with everything else
    source = _bench_source()
    n_files = length(source.files)
    n_items = KIND1_FILES + KIND2_FILES + KIND3_FILES * KIND3_CYCLES
    data_bytes = filesize(_template("kind1.csv")) +
        filesize(_template("kind2.csv")) +
        filesize(_template("kind3.csv"))

    default_outdir = joinpath(@__DIR__, "results",
        "realistic-" * replace(string(round(Int, time())), r"\D" => ""))
    outdir = get(ENV, "MB_BENCH_OUTDIR", default_outdir)
    mkpath(outdir)
    log_path = joinpath(outdir, "benchmark.log")
    log_io = open(log_path, "w")
    RUN_LOG[] = log_io
    Profiling.reset_debug_timings!()

    try
        _print_run_header(log_path, outdir)

        tee_printf("Templates: 3 files, %.1f MB on disk; %d aliases, ~%d items (scale=%s)\n",
            data_bytes / 1024^2, n_files, n_items, string(scale))

        project = build_project(; plots=true)
        kinds = (:kind1, :kind2, :kind3)
        plot_kinds = Dict(k => first(registered_plot_kinds(project, k)) for k in kinds)
        samples = Sample[]

        tee_println("Building cache + browsing during the scan ...")
        ws = open_workspace(project, source)
        t_start = time()
        build_seconds = 0.0
        scan_seconds = 0.0
        processing_started = nothing
        processing_seconds = 0.0
        analysis_started = nothing
        analysis_seconds = 0.0
        build_stats = nothing
        saturation_stats = nothing
        try
            last_probe = 0.0
            kind_cursor = 1
            while true
                now = time() - t_start
                processing_active = _processing_active(ws)
                analysis_active = _analysis_active(ws)
                processing_started === nothing && processing_active && (processing_started = now)
                scan_seconds == 0 && !Workspace.source_scan_running(ws) &&
                    ws.scan.state in (:done, :unchanged, :error, :canceled) &&
                    (scan_seconds = now)
                if processing_started !== nothing && processing_seconds == 0 && !processing_active
                    processing_seconds = now - processing_started
                    analysis_started === nothing && (analysis_started = now)
                end
                analysis_started === nothing && analysis_active && (analysis_started = now)
                analysis_started !== nothing && analysis_seconds == 0 &&
                    !analysis_active &&
                    (analysis_seconds = now - analysis_started)
                # Probe responsiveness ~6×/s, rotating across kinds, once items exist. Each probe is the
                # full select → load → plot probe a user performs while the build is still running.
                if now - last_probe >= 0.16
                    last_probe = now
                    kind = kinds[kind_cursor]; kind_cursor = mod1(kind_cursor + 1, length(kinds))
                    probe = timed_plot!(ws, plot_kinds, kind, 3)
                    probe === nothing || push!(samples, Sample(now, :during_build, kind, probe.n,
                        probe.plot_ms, probe.bytes, probe.ready))
                end
                if build_idle(ws)
                    build_seconds = now
                    break
                end
                (now > MAX_BUILD_SECONDS) && (build_seconds = now;
                    @warn("hit MAX_BUILD_SECONDS"); break)
                sleep(0.004)
            end

            tee_println("Saturating processed-payload writer ...")
            saturation_stats = saturate_processed_writes!(ws, :kind3)

            # Steady-state sweep: random plot probes per kind on the finished cache.
            rng = MersenneTwister(1)
            for kind in kinds, _ in 1:AFTER_BUILD_PLOTS
                ids = [
                    id for (id, record) in ws.index.items
                    if label(record.type) === kind
                ]
                isempty(ids) && continue
                k = rand(rng, 1:min(4, length(ids)))
                records = ItemRecord[ws.index.items[id] for id in rand(rng, ids, k)]
                probe = timed_plot!(ws, plot_kinds, kind, k; records)
                probe === nothing || push!(samples, Sample(
                    time() - t_start,
                    :after_build,
                    kind,
                    probe.n,
                    probe.plot_ms,
                    probe.bytes,
                    probe.ready,
                ))
            end
            completed, total, active = Workspace.work_counts(ws)
            collections = ws.index.collections
            collection_nodes = count(keys(collections.records)) do collection_key
                member_ids = collection_item_ids(collections, collection_key)
                isempty(member_ids) && return false
                key = Workspace.WorkKey(Workspace.COLLECTION_ANALYZE, collection_key)
                return Workspace.cache_work_status(ws, key) === :ready
            end
            metrics = ws.metrics
            build_stats = (
                scan_seconds,
                processing_seconds,
                analysis_seconds,
                processed_items=length(ws.index.items),
                completed_jobs=completed,
                total_jobs=total,
                active_jobs=active,
                collection_nodes=collection_nodes,
                interpreted_write_ns=metrics.interpreted_write_ns[],
                interpreted_writes=metrics.interpreted_writes[],
                processed_write_ns=metrics.processed_write_ns[],
                processed_writes=metrics.processed_writes[],
                metadata_write_ns=metrics.metadata_write_ns[],
                metadata_writes=metrics.metadata_writes[],
            )
        finally
            close_workspace!(ws)
        end

        # Warm reopen on the same cache: the incremental rescan finds every fingerprint unchanged and
        # reuses the cached index. Surfaces the true warm-reopen cost (rescan + cached-index handling +
        # any re-processing the post-scan readiness probe triggers).
        reopen_stats = measure_reopen(project, source, plot_kinds, kinds)

        report(samples, build_stats, saturation_stats, reopen_stats, outdir,
            n_files, n_items, data_bytes, build_seconds)

        tee_println("\nResults kept in: $outdir")
        tee_println("Log kept in: $log_path")
        timings = Profiling.take_debug_timings!()
        write_debug_timings(outdir, timings)
        return outdir
    finally
        RUN_LOG[] = nothing
        close(log_io)
        first(DEPOT_PATH) == tmp && popfirst!(DEPOT_PATH)
        ispath(tmp) && rm(tmp; force=true, recursive=true)
    end
end

# --------------------------------------------------------------------------------------------------
# Warm reopen
# --------------------------------------------------------------------------------------------------

"""One close-and-reopen on the warm cache: time to first view and idle, allocation, first plots."""
function _reopen_once(project, source, plot_kinds, kinds)
    GC.gc()
    t0 = time()
    bytes0 = Base.gc_bytes()
    ws = open_workspace(project, source)
    first_view_s = 0.0
    idle_s = 0.0
    deadline = time() + MAX_BUILD_SECONDS
    try
        while true
            now = time() - t0
            first_view_s == 0 && !isempty(ws.index.items) && (first_view_s = now)
            if build_idle(ws)
                idle_s = now
                break
            end
            time() > deadline && (idle_s = now; @warn("reopen hit MAX_BUILD_SECONDS"); break)
            sleep(0.002)
        end
        alloc_bytes = Base.gc_bytes() - bytes0
        first_plots = NamedTuple[]
        for kind in kinds
            probe = timed_plot!(ws, plot_kinds, kind, 1)
            probe === nothing || push!(first_plots,
                (kind=kind, plot_ms=probe.plot_ms))
        end
        return (first_view_s=first_view_s, idle_s=idle_s, alloc_bytes=alloc_bytes,
            items=length(ws.index.items), unchanged=ws.scan.state === :unchanged,
            first_plots=first_plots)
    finally
        close_workspace!(ws)
    end
end

"""
Close-and-reopen on the same warm cache and time the incremental rescan and first plot.

A first discarded pass warms the reopen-specific code paths (cache-index load, incremental reuse,
cached-index handling) so the reported allocation reflects work, not first-call compilation. Reports
the wall time to first cached view and to idle, the total bytes allocated getting to idle (the whole
warm-reopen cost — rescan, cached-index handling, and any re-processing the readiness probe triggers),
and the first warm plot per kind (its data is read from disk, never the staged buffer).
"""
function measure_reopen(project, source, plot_kinds, kinds)
    _reopen_once(project, source, plot_kinds, kinds)   # warm up JIT, discard
    return _reopen_once(project, source, plot_kinds, kinds)
end

# --------------------------------------------------------------------------------------------------
# Reporting
# --------------------------------------------------------------------------------------------------

function report(samples, stats, saturation, reopen, outdir,
    n_files, n_items, data_bytes, build_seconds)
    saturation === nothing && error("Missing processed-writer saturation sample")
    open(joinpath(outdir, "responsiveness.csv"), "w") do io
        println(io, "elapsed_s,phase,kind,n_items,plot_ms,allocated_bytes,ready_items")
        for s in samples
            @printf(io, "%.3f,%s,%s,%d,%.3f,%d,%d\n",
                s.elapsed_s, s.phase, s.kind, s.n, s.plot_ms,
                s.allocated_bytes, s.ready)
        end
    end

    open(joinpath(outdir, "saturation.csv"), "w") do io
        println(io, "kind,requested_items,materialized_items,estimated_rows,load_ms,flush_ms,peak_pending_rows,processed_writes,row_ceiling")
        @printf(io, "%s,%d,%d,%d,%.3f,%.3f,%d,%d,%d\n",
            saturation.kind,
            saturation.requested_items,
            saturation.materialized_items,
            saturation.estimated_rows,
            saturation.load_ms,
            saturation.flush_ms,
            saturation.peak_pending_rows,
            saturation.processed_writes,
            CACHE_ROW_CEILING)
    end

    open(joinpath(outdir, "reopen.csv"), "w") do io
        println(io, "first_view_s,idle_s,alloc_mib,items,unchanged")
        @printf(io, "%.3f,%.3f,%.1f,%d,%s\n", reopen.first_view_s, reopen.idle_s,
            reopen.alloc_bytes / 1024^2, reopen.items, reopen.unchanged)
        println(io, "kind,first_plot_ms")
        for p in reopen.first_plots
            @printf(io, "%s,%.3f\n", p.kind, p.plot_ms)
        end
    end

    write_calls = stats.interpreted_writes + stats.processed_writes + stats.metadata_writes
    write_ns = stats.interpreted_write_ns + stats.processed_write_ns + stats.metadata_write_ns
    mean_write_ms = write_calls == 0 ? 0.0 : write_ns / write_calls / 1e6
    if REQUIRE_SATURATION
        stats.interpreted_writes > 0 || error("Benchmark did not exercise interpreted writes")
        stats.processed_writes > 0 || error("Benchmark did not exercise processed writes")
        stats.metadata_writes > 0 || error("Benchmark did not exercise metadata writes")
        saturation.estimated_rows >= CACHE_ROW_CEILING || error(
            "Processed-writer saturation selected only $(saturation.estimated_rows) rows, " *
            "below the cache row ceiling $CACHE_ROW_CEILING",
        )
        saturation.processed_writes > 0 || error(
            "Processed-writer saturation created no processed writes; selected items were " *
            "already cached or memory-resident",
        )
    end
    during = [s.plot_ms for s in samples if s.phase === :during_build]
    after = [s.plot_ms for s in samples if s.phase === :after_build]
    read_stat(values, statistic) = isempty(values) ? NaN : statistic(values)
    source_files = max(n_files, 1)
    indexed_items = max(stats.processed_items, 1)
    payload_rows = max(ESTIMATED_PAYLOAD_ROWS, 1)
    per_second(count, seconds) = seconds > 0 ? count / seconds : NaN
    rows_per_file = ESTIMATED_PAYLOAD_ROWS / source_files
    rows_per_item = ESTIMATED_PAYLOAD_ROWS / max(n_items, 1)
    build_ms_per_file = build_seconds * 1e3 / source_files
    build_ms_per_item = build_seconds * 1e3 / indexed_items
    scan_ms_per_file = stats.scan_seconds * 1e3 / source_files
    processing_ms_per_item = stats.processing_seconds * 1e3 / indexed_items
    write_ms_per_file = write_ns / 1e6 / source_files
    write_ms_per_item = write_ns / 1e6 / indexed_items
    write_ns_per_payload_row = write_ns / payload_rows
    open(joinpath(outdir, "scorecard.csv"), "w") do io
        println(io, "source_files,items,estimated_payload_rows,data_mib," *
            "rows_per_file,rows_per_item,build_s,scan_s,processing_s,analysis_s," *
            "scan_files_per_s,processing_items_per_s,build_items_per_s," *
            "build_ms_per_file,build_ms_per_item,scan_ms_per_file,processing_ms_per_item," *
            "write_ms_per_call,write_ms_per_file,write_ms_per_item,write_ns_per_payload_row," *
            "saturation_items,saturation_rows,saturation_load_ms,saturation_flush_ms," *
            "saturation_peak_pending_rows,saturation_processed_writes," *
            "during_plot_median_ms,during_plot_p90_ms,during_plot_p99_ms,during_plot_max_ms," *
            "after_plot_median_ms,after_plot_p90_ms,after_plot_p99_ms,after_plot_max_ms")
        @printf(io, "%d,%d,%d,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%d,%d,%.3f,%.3f,%d,%d,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f\n",
            n_files,
            stats.processed_items,
            ESTIMATED_PAYLOAD_ROWS,
            data_bytes / 1024^2,
            rows_per_file,
            rows_per_item,
            build_seconds,
            stats.scan_seconds,
            stats.processing_seconds,
            stats.analysis_seconds,
            per_second(n_files, stats.scan_seconds),
            per_second(stats.processed_items, stats.processing_seconds),
            per_second(stats.processed_items, build_seconds),
            build_ms_per_file,
            build_ms_per_item,
            scan_ms_per_file,
            processing_ms_per_item,
            mean_write_ms,
            write_ms_per_file,
            write_ms_per_item,
            write_ns_per_payload_row,
            saturation.requested_items, saturation.estimated_rows, saturation.load_ms,
            saturation.flush_ms, saturation.peak_pending_rows, saturation.processed_writes,
            read_stat(during, median), read_stat(during, values -> quantile(values, 0.9)),
            read_stat(during, values -> quantile(values, 0.99)), read_stat(during, maximum),
            read_stat(after, median), read_stat(after, values -> quantile(values, 0.9)),
            read_stat(after, values -> quantile(values, 0.99)), read_stat(after, maximum))
    end

    tee_println("\n==================== REALISTIC BROWSE BENCHMARK ====================")
    tee_printf("dataset:  %d aliases · ~%d items · %d estimated rows · %.1f MB templates\n",
        n_files, n_items, ESTIMATED_PAYLOAD_ROWS, data_bytes / 1024^2)
    tee_printf("build:    %.1f s wall (scan + processing + collection analysis)\n", build_seconds)
    n_during = count(s -> s.phase === :during_build, samples)
    tee_printf("plot probes: %d during build · %d after (secondary CSV only)\n", n_during,
        count(s -> s.phase === :after_build, samples))

    tee_println("\nThroughput:")
    tee_printf("  scan                %8.1f source items/s  (%6.1f s)\n",
        per_second(n_files, stats.scan_seconds), stats.scan_seconds)
    tee_printf("  item processing     %8.1f items/s         (%6.1f s, %d unique items)\n",
        per_second(stats.processed_items, stats.processing_seconds),
        stats.processing_seconds, stats.processed_items)
    stats.completed_jobs == stats.processed_items || tee_printf(
        "  duplicate queue work %8d cache-hit jobs\n",
        stats.completed_jobs - stats.processed_items,
    )
    tee_printf("  collection analysis %8.1f nodes/s         (%6.1f s, %d nodes)\n",
        per_second(stats.collection_nodes, stats.analysis_seconds),
        stats.analysis_seconds, stats.collection_nodes)
    tee_printf("  build average       %8.1f items/s\n",
        per_second(stats.processed_items, build_seconds))

    tee_println("\nNormalized averages:")
    tee_printf("  payload shape       %8.0f rows/file  %8.0f rows/item\n",
        rows_per_file, rows_per_item)
    tee_printf("  build               %8.2f ms/file    %8.2f ms/item\n",
        build_ms_per_file, build_ms_per_item)
    tee_printf("  scan/process        %8.2f ms/file    %8.2f ms/item\n",
        scan_ms_per_file, processing_ms_per_item)
    tee_printf("  writes              %8.2f ms/file    %8.2f ms/item  %8.1f ns/row\n",
        write_ms_per_file, write_ms_per_item, write_ns_per_payload_row)

    tee_println("\nWrites:")
    tee_printf("  interpreted %6d calls  mean %7.2f ms\n", stats.interpreted_writes,
        stats.interpreted_write_ns / max(stats.interpreted_writes, 1) / 1e6)
    tee_printf("  processed   %6d calls  mean %7.2f ms  mean batch %5.1f items\n",
        stats.processed_writes,
        stats.processed_write_ns / max(stats.processed_writes, 1) / 1e6,
        stats.processed_items / max(stats.processed_writes, 1))
    tee_printf("  stats       %6d calls  mean %7.2f ms\n", stats.metadata_writes,
        stats.metadata_write_ns / max(stats.metadata_writes, 1) / 1e6)
    tee_printf("  combined mean %7.2f ms\n", mean_write_ms)

    tee_println("\nProcessed-writer saturation:")
    tee_printf("  selected %d %s items  estimated rows %d  row ceiling %d\n",
        saturation.requested_items, saturation.kind, saturation.estimated_rows,
        CACHE_ROW_CEILING)
    tee_printf("  materialize %.1f ms  flush %.1f ms  peak pending rows %d  processed writes %d\n",
        saturation.load_ms, saturation.flush_ms, saturation.peak_pending_rows,
        saturation.processed_writes)

    tee_println("\nWarm reopen (same cache, every fingerprint unchanged):")
    tee_printf("  first cached view %.0f ms  ·  idle %.2f s  ·  allocated %.1f MiB  ·  %d items%s\n",
        reopen.first_view_s * 1e3, reopen.idle_s, reopen.alloc_bytes / 1024^2, reopen.items,
        reopen.unchanged ? "  (reused)" : "  (re-scanned!)")
    for p in reopen.first_plots
        tee_printf("  first %-6s plot  %6.1f ms\n", p.kind, p.plot_ms)
    end

    return nothing
end

for repeat in 1:BENCH_REPEATS
    BENCH_REPEATS == 1 || println("\n===== benchmark repeat $repeat / $BENCH_REPEATS =====\n")
    run_benchmark()
end
