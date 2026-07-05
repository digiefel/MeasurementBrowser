# Realistic browse-while-building benchmark.
#
# Models the pressure profile of a real project (the RuO2 v2 project): fewer source files than the
# full data set, but realistic fatigue-style row volume. The important scale is staged DataFrame rows:
# the real run crosses the cache buffer row ceiling and forces the writer/backpressure path. The
# default workload exceeds the real per-source/item burst shape and adds a processed-payload stress
# pass, while staying small enough for routine engine checks.
#
#   julia --project=bench --threads=auto bench/realistic_browse.jl [scale]
#
# `scale` (default 1.0) multiplies file/item counts, not rows per item. Synthetic data and the cache
# live in a temp dir that is deleted on exit, so the drive is not filled; only the small result files
# (CSVs + PNG) are kept under bench/results/.
#
# Tunables via ENV (counts are per the documented diversity; see DEFAULTS below):
#   MB_BENCH_KIND1_FILES, MB_BENCH_KIND2_FILES, MB_BENCH_KIND3_FILES, MB_BENCH_KIND3_CYCLES
#   MB_BENCH_KIND3_ROWS, MB_BENCH_AFTER_BUILD_PLOTS, MB_BENCH_PROCESSED_STRESS_ROWS
#   MB_BENCH_REQUIRE_SATURATION
# The benchmark records a structured engine profile by default. Set MB_PROFILE_INTERNAL=0 to disable
# it; add MB_PROFILE_CPU=1 for Julia sampling. Set MB_BENCH_START_PROFILE=0 to measure
# enabled-but-idle overhead.

using MeasurementBrowser
const MB = MeasurementBrowser
const ENGINE_ONLY = MB.Profiling.environment_flag("MB_BENCH_ENGINE_ONLY", false)
using CSV
using DataFrames
using Dates
using Random
using Printf
using Statistics: mean, median, quantile
# GLMakie ships with the engine; the plot recipes below build real figures from the cached data, the
# same scene-graph work the GUI does on a selection (no window is shown, so no GPU is touched).
ENGINE_ONLY || @eval import GLMakie: Figure, Axis, lines!, contents

# Optional: render a PNG when CairoMakie is in the bench env. Imported at top level so the plotting
# call below runs in a new-enough world age (importing inside the function fails with "method too new").
const HAS_CAIRO = !ENGINE_ONLY && try
    @eval import CairoMakie
    true
catch
    false
end

# --------------------------------------------------------------------------------------------------
# Sizing (downscaled from the real RuO2 project, preserving row-pressure shape)
# --------------------------------------------------------------------------------------------------

scale = length(ARGS) >= 1 ? parse(Float64, ARGS[1]) : 1.0
_env_int(key, default) = parse(Int, get(ENV, key, string(default)))
_scaled(n) = max(1, round(Int, n * scale))
const PROFILE_INTERNAL = MB.Profiling.environment_flag("MB_PROFILE_INTERNAL", true)
const PROFILE_CPU = MB.Profiling.environment_flag("MB_PROFILE_CPU")
const START_PROFILE = MB.Profiling.environment_flag(
    "MB_BENCH_START_PROFILE", PROFILE_INTERNAL)
START_PROFILE && !PROFILE_INTERNAL && error(
    "MB_BENCH_START_PROFILE=1 requires MB_PROFILE_INTERNAL=1",
)

# kind1: tiny files, one item each (IV-style, ~120 rows).      many small records
# kind2: medium files, one item each (CV-style, ~5000 rows).   mid-size table writes
# kind3: big files, one item PER CYCLE (fatigue-style).        few files, many rows/items
# Defaults target a short run while exceeding the real fatigue burst: one synthetic big file expands
# into more cycles and rows per cycle than the typical RuO2 fatigue file, but there are far fewer
# physical files.
const KIND1_FILES  = _scaled(_env_int("MB_BENCH_KIND1_FILES", 500))
const KIND2_FILES  = _scaled(_env_int("MB_BENCH_KIND2_FILES", 120))
const KIND3_FILES  = _scaled(_env_int("MB_BENCH_KIND3_FILES", 16))
const KIND3_CYCLES = _env_int("MB_BENCH_KIND3_CYCLES", 96)   # items per big file
const KIND1_ROWS   = 120
const KIND2_ROWS   = 5_000
const KIND3_ROWS   = _env_int("MB_BENCH_KIND3_ROWS", 6_000) # rows per cycle

const CACHE_ROW_CEILING = MB.Cache.CACHE_BUFFER_ROW_LIMIT
const PROCESSED_STRESS_ROWS = _env_int(
    "MB_BENCH_PROCESSED_STRESS_ROWS",
    2 * CACHE_ROW_CEILING,
)
const AFTER_BUILD_PLOTS = _env_int("MB_BENCH_AFTER_BUILD_PLOTS", 40)
const BENCH_REPEATS = _env_int("MB_BENCH_REPEATS", 1)
const MAX_BUILD_SECONDS = 600   # safety cap
const REQUIRE_SATURATION = MB.Profiling.environment_flag(
    "MB_BENCH_REQUIRE_SATURATION",
    scale >= 1.0,
)
const ESTIMATED_PAYLOAD_ROWS = Int64(KIND1_FILES) * KIND1_ROWS +
                               Int64(KIND2_FILES) * KIND2_ROWS +
                               Int64(KIND3_FILES) * KIND3_CYCLES * KIND3_ROWS

const RUN_LOG = Ref{Union{Nothing,IO}}(nothing)
const BENCH_ENV_KEYS = (
    "MB_BENCH_KIND1_FILES",
    "MB_BENCH_KIND2_FILES",
    "MB_BENCH_KIND3_FILES",
    "MB_BENCH_KIND3_CYCLES",
    "MB_BENCH_KIND3_ROWS",
    "MB_BENCH_AFTER_BUILD_PLOTS",
    "MB_BENCH_ENGINE_ONLY",
    "MB_BENCH_PROCESSED_STRESS_ROWS",
    "MB_BENCH_REPEATS",
    "MB_BENCH_REQUIRE_SATURATION",
    "MB_BENCH_START_PROFILE",
    "MB_PROFILE_INTERNAL",
    "MB_PROFILE_CPU",
    "MB_PROFILE_OUTPUT",
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
    tee_println("MeasurementBrowser realistic benchmark")
    tee_println("started_at: ", Dates.format(now(), dateformat"yyyy-mm-dd HH:MM:SS"))
    tee_println("output_dir: ", outdir)
    tee_println("log_file:   ", log_path)
    tee_println("branch:     ", _repo_command(["git", "rev-parse", "--abbrev-ref", "HEAD"]))
    tee_println("commit:     ", _repo_command(["git", "rev-parse", "HEAD"]))
    status = _repo_command(["git", "status", "--short", "--", "bench/realistic_browse.jl", "bench/README.md"])
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
# Project: three kinds mirroring the real read/entries/process/analyze/plot shape
# --------------------------------------------------------------------------------------------------

_collection(file) = [splitpath(dirname(file.filepath))[end-1], splitpath(dirname(file.filepath))[end]]

function build_project(; plots::Bool=true)
    project = MB.define_project("BenchRealistic"; description="Realistic browse-while-build benchmark")

    MB.register_item!(project, :kind1;
        detect  = file -> endswith(file.filename, "_kind1.csv"),
        read    = file -> DataFrame(CSV.File(file.filepath; ntasks=1)),
        entries = (file, data) -> [MB.DataItem(; kind=:kind1, collection=_collection(file),
            label=file.filename, data=data, id=file.filepath * "#kind1")],
        process = item -> MB.DataItem(item, transform(item.data, [:voltage, :current] =>
            ByRow((v, i) -> iszero(v) ? 0.0 : i / v) => :conductance)),
        analyze = item -> Dict{Symbol,Any}(:imax => maximum(abs, item.data.current)),
        label   = item -> "K1 $(item.label)")

    MB.register_item!(project, :kind2;
        detect  = file -> endswith(file.filename, "_kind2.csv"),
        read    = file -> DataFrame(CSV.File(file.filepath; ntasks=1)),
        entries = (file, data) -> [MB.DataItem(; kind=:kind2, collection=_collection(file),
            label=file.filename, data=data, id=file.filepath * "#kind2")],
        analyze = item -> Dict{Symbol,Any}(:cmean => mean(item.data.cap)),
        label   = item -> "K2 $(item.label)")

    # The fatigue-style kind: one file → one item per cycle (where item count explodes).
    MB.register_item!(project, :kind3;
        detect  = file -> endswith(file.filename, "_kind3.csv"),
        read    = file -> DataFrame(CSV.File(file.filepath; ntasks=1)),
        entries = (file, data) -> [
            MB.DataItem(; kind=:kind3, collection=_collection(file),
                label=string(file.filename, " cycle ", c), metadata=Dict{Symbol,Any}(:cycle => c),
                data=(@view data[data.cycle .== c, :]), id=string(file.filepath, "#kind3,cycle=", c))
            for c in sort(unique(data.cycle))],
        process = item -> MB.DataItem(item, transform(item.data, [:voltage, :current] =>
            ByRow((v, i) -> v * i) => :power)),
        analyze = item -> Dict{Symbol,Any}(:pmax => maximum(abs, item.data.current)),
        label   = item -> "K3 $(item.label)")

    plots || return project

    # Real plot recipes: one axis, one line per selected item, reading the processed columns — the
    # same scene-graph work the GUI does when a user selects items and the plot panel renders.
    _axis(fig, xlabel, ylabel) = Axis(fig[1, 1]; xlabel, ylabel)
    MB.register_plot!(project, :kind1; label="IV",
        setup=(ws, items) -> (fig = Figure(); _axis(fig, "voltage", "conductance"); fig),
        draw=(ws, items, fig) -> for item in items
            d = MB.item_data(item)
            lines!(contents(fig[1, 1])[1], d.voltage, d.conductance)
        end)
    MB.register_plot!(project, :kind2; label="CV",
        setup=(ws, items) -> (fig = Figure(); _axis(fig, "voltage", "cap"); fig),
        draw=(ws, items, fig) -> for item in items
            d = MB.item_data(item)
            lines!(contents(fig[1, 1])[1], d.voltage, d.cap)
        end)
    MB.register_plot!(project, :kind3; label="power",
        setup=(ws, items) -> (fig = Figure(); _axis(fig, "time", "power"); fig),
        draw=(ws, items, fig) -> for item in items
            d = MB.item_data(item)
            lines!(contents(fig[1, 1])[1], d.time, d.power)
        end)
    return project
end

# --------------------------------------------------------------------------------------------------
# Driver: poll the workspace while timing interactive reads (selection + plot data load)
# --------------------------------------------------------------------------------------------------

build_idle(ws) = begin
    _, _, active = MB.Workspace.work_counts(ws)
    return active == 0 &&
        !MB.Workspace.source_scan_running(ws) &&
        !MB.Workspace.processing_work_running(ws) &&
        !MB.Workspace.analysis_work_running(ws) &&
        !MB.Workspace.cache_work_running(ws) &&
        ws.scan.state in (:done, :unchanged, :error, :canceled)
end

function _active_work_count(ws, kinds::Tuple)::Int
    return lock(ws.work.lock) do
        count(
            node -> node.key.kind in kinds && node.state in (:queued, :running),
            values(ws.work.nodes),
        )
    end
end

_processing_active(ws)::Bool =
    _active_work_count(ws, (MB.Workspace.ITEM_PROCESS, MB.Workspace.ITEM_ANALYZE)) > 0 ||
    MB.Workspace.cache_has_pending_writes(ws.cache.db)

_analysis_active(ws)::Bool =
    _active_work_count(ws, (MB.Workspace.COLLECTION_ANALYZE,)) > 0

"""Collect up to 64 already-processed item ids of `kind`."""
function _ready_ids(ws, kind::Symbol)
    ready = String[]
    for id in keys(ws.index.item_metadata)
        rec = get(ws.index.items, id, nothing)
        rec === nothing && continue
        rec.kind === kind && push!(ready, id)
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
        records = MB.ItemIndex.ItemRecord[ws.index.items[id] for id in ready[1:k]]
    end
    n_ready = records === nothing ? 0 : length(records)
    result = @timed begin
        MB.select_items!(ws, records)             # mirror the GUI selecting them
        items = MB.Workspace.materialize_items(ws, records)
        plot_kind = plot_kinds[kind]
        figure = MB.setup_plot(ws, plot_kind, items)
        MB.plot_data!(ws, plot_kind, items, figure)
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

mutable struct MemorySample
    elapsed_s::Float64
    rss_bytes::Int64
    gc_live_bytes::Int64
    gc_allocated_bytes::Int64
    rss_minus_gc_live_bytes::Int64
    index_items::Int64
    index_collections::Int64
    item_metadata::Int64
    analysis_errors::Int64
    processing_jobs::Int64
    pending_writes::Int64
    pending_write_rows::Int64
    selected_queue::Int64
    background_waiting::Int64
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

function MemorySample(elapsed_s::Float64, snapshot)::MemorySample
    return MemorySample(
        elapsed_s,
        snapshot.rss_bytes,
        snapshot.gc_live_bytes,
        snapshot.gc_allocated_bytes,
        snapshot.rss_minus_gc_live_bytes,
        snapshot.index_items,
        snapshot.index_collections,
        snapshot.item_metadata,
        snapshot.analysis_errors,
        snapshot.processing_jobs,
        snapshot.pending_writes,
        snapshot.pending_write_rows,
        snapshot.selected_queue,
        snapshot.background_waiting,
    )
end

function _records_of_kind(ws, kind::Symbol)::Vector{MB.ItemIndex.ItemRecord}
    records = MB.ItemIndex.ItemRecord[
        record for record in values(ws.index.items)
        if record.kind === kind && haskey(ws.index.item_metadata, record.id)
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
            "Increase MB_BENCH_KIND3_FILES, MB_BENCH_KIND3_CYCLES, or MB_BENCH_KIND3_ROWS.",
        )
    end

    MB.select_items!(ws, selected)
    processed_writes_before = ws.metrics.processed_writes[]
    load = @timed MB.Workspace.materialize_items(ws, selected)
    peak_pending_rows = Int64(0)
    flush_started = time()
    while MB.Workspace.cache_has_pending_writes(ws.cache.db)
        counts = MB.Workspace.cache_pending_counts(ws.cache.db)
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

_event_ms(event)::Float64 = event.duration_ns / 1e6
_event_wait_ms(event)::Float64 = event.attributes.wait_ns / 1e6

function _push_event_times!(
    rows::Vector{NamedTuple},
    metric::String,
    events::Vector,
    value::Function,
)::Nothing
    for event in events
        ms = value(event)
        isfinite(ms) || continue
        push!(rows, (metric=metric, ms=Float64(max(ms, 0.0))))
    end
    return nothing
end

function _matching_events(report, category::Symbol, operation::Symbol)::Vector
    return [event for event in report.events
            if event.category === category && event.operation === operation]
end

function _profile_event_groups(report)::Dict{Tuple{Symbol,Symbol},Vector}
    groups = Dict{Tuple{Symbol,Symbol},Vector}()
    for event in report.events
        push!(get!(() -> Any[], groups, (event.category, event.operation)), event)
    end
    return groups
end

function _child_duration_ms(report)::Dict{UInt64,Float64}
    by_parent = Dict{UInt64,Float64}()
    for event in report.events
        event.parent_id == 0 && continue
        by_parent[event.parent_id] = get(by_parent, event.parent_id, 0.0) + _event_ms(event)
    end
    return by_parent
end

function pipeline_event_time_rows(report)::Vector{NamedTuple}
    child_ms = _child_duration_ms(report)
    rows = NamedTuple[]
    _push_event_times!(rows, "processing_queue_wait_ms",
        _matching_events(report, :processing, :item), _event_wait_ms)
    _push_event_times!(rows, "processing_engine_overhead_ms",
        _matching_events(report, :processing, :item), event -> begin
            _event_ms(event) - get(child_ms, event.id, 0.0)
        end)
    _push_event_times!(rows, "interpret_engine_overhead_ms",
        _matching_events(report, :project, :interpret_source_item), event -> begin
            _event_ms(event) - get(child_ms, event.id, 0.0)
        end)
    for operation in sort!(unique(event.operation for event in report.events
                                  if event.category === :cache &&
                                     startswith(String(event.operation), "flush_")))
        _push_event_times!(rows, "cache_$(operation)_ms",
            _matching_events(report, :cache, operation), _event_ms)
    end
    return rows
end

function _write_pipeline_event_summary!(rows::Vector{NamedTuple}, outdir::String)::Nothing
    groups = Dict{String,Vector{Float64}}()
    for row in rows
        push!(get!(() -> Float64[], groups, row.metric), row.ms)
    end
    open(joinpath(outdir, "pipeline_event_summary.csv"), "w") do io
        println(io, "metric,count,p50_ms,p90_ms,p99_ms,max_ms")
        for metric in sort!(collect(keys(groups)))
            values = groups[metric]
            @printf(io, "%s,%d,%.3f,%.3f,%.3f,%.3f\n",
                metric,
                length(values),
                median(values),
                quantile(values, 0.9),
                quantile(values, 0.99),
                maximum(values))
        end
    end
    return nothing
end

function _write_pipeline_event_times!(report, outdir::String)::Vector{NamedTuple}
    rows = pipeline_event_time_rows(report)
    table = isempty(rows) ? DataFrame(metric=String[], ms=Float64[]) : DataFrame(rows)
    CSV.write(joinpath(outdir, "pipeline_event_times.csv"), table)
    _write_pipeline_event_summary!(rows, outdir)
    return rows
end

function _write_pipeline_timeseries!(memory_samples::Vector{MemorySample}, outdir::String)::Nothing
    CSV.write(joinpath(outdir, "pipeline_timeseries.csv"), DataFrame(memory_samples))
    return nothing
end

function run_benchmark()
    tmp = mktempdir()
    pushfirst!(DEPOT_PATH, tmp)          # cache lands in temp, deleted with everything else
    data_root = joinpath(tmp, "data")

    outdir = joinpath(@__DIR__, "results",
        "realistic-" * replace(string(round(Int, time())), r"\D" => ""))
    mkpath(outdir)
    log_path = joinpath(outdir, "benchmark.log")
    log_io = open(log_path, "w")
    RUN_LOG[] = log_io

    try
        _print_run_header(log_path, outdir)

        tee_println("Generating synthetic data ... (scale=$scale)")
        gen_t = @elapsed (n_files, n_items) = generate_data(data_root)
        data_bytes = sum(filesize(joinpath(r, f))
                         for (r, _, fs) in walkdir(data_root) for f in fs)
        tee_printf("  %d files, ~%d items, %.1f MB on disk, generated in %.1fs\n",
            n_files, n_items, data_bytes / 1024^2, gen_t)

    project = build_project(; plots=!ENGINE_ONLY)
    kinds = (:kind1, :kind2, :kind3)
    plot_kinds = ENGINE_ONLY ?
        Dict{Symbol,Any}() :
        Dict(k => first(MB.registered_plot_kinds(project, k)) for k in kinds)
    samples = Sample[]

    tee_println("Building cache + browsing during the scan ...")
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
    build_stats = nothing
    saturation_stats = nothing
    profile_report = nothing
    try
        last_probe = 0.0
        last_rss_sample = 0.0
        last_memory_sample = -Inf
        memory_samples = MemorySample[]
        kind_cursor = 1
        while true
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
            processing_active = _processing_active(ws)
            analysis_active = _analysis_active(ws)
            processing_started === nothing && processing_active && (processing_started = now)
            scan_seconds == 0 && !MB.Workspace.source_scan_running(ws) &&
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
            if !ENGINE_ONLY && now - last_probe >= 0.16
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
        final_memory = MB.Workspace.workspace_memory_snapshot(ws)
        rss_end_bytes = final_memory.rss_bytes
        push!(memory_samples, MemorySample(time() - t_start, final_memory))

        tee_println("Saturating processed-payload writer ...")
        saturation_stats = saturate_processed_writes!(ws, :kind3)

        # Steady-state sweep: random plot probes per kind on the finished cache.
        if !ENGINE_ONLY
            rng = MersenneTwister(1)
            for kind in kinds, _ in 1:AFTER_BUILD_PLOTS
                ids = [id for id in keys(ws.index.item_metadata)
                       if (r = get(ws.index.items, id, nothing); r !== nothing && r.kind === kind)]
                isempty(ids) && continue
                k = rand(rng, 1:min(4, length(ids)))
                records = MB.ItemIndex.ItemRecord[ws.index.items[id] for id in rand(rng, ids, k)]
                probe = timed_plot!(ws, plot_kinds, kind, k; records=records)
                probe === nothing || push!(samples, Sample(time() - t_start, :after_build, kind, probe.n,
                    probe.plot_ms, probe.bytes, probe.ready))
            end
        end
        completed, queued, _active = MB.Workspace.work_counts(ws)
        collection_nodes = lock(ws.work.lock) do
            count(
                node -> node.key.kind === MB.Workspace.COLLECTION_ANALYZE &&
                    node.state in (:ready, :failed),
                values(ws.work.nodes),
            )
        end
        metrics = ws.metrics
        build_stats = (
            scan_seconds,
            processing_seconds,
            analysis_seconds,
            processed_items=length(ws.index.items),
            completed_jobs=completed,
            queued_items=queued,
            collection_nodes=collection_nodes,
            interpreted_write_ns=metrics.interpreted_write_ns[],
            interpreted_writes=metrics.interpreted_writes[],
            processed_write_ns=metrics.processed_write_ns[],
            processed_writes=metrics.processed_writes[],
            metadata_write_ns=metrics.metadata_write_ns[],
            metadata_writes=metrics.metadata_writes[],
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

    # Warm reopen on the same cache: the incremental rescan finds every fingerprint unchanged and
    # reuses the cached index. Surfaces the true warm-reopen cost (rescan + cached-index handling +
    # any re-processing the post-scan readiness probe triggers).
    reopen_stats = measure_reopen(project, data_root, plot_kinds, kinds)

    report(samples, build_stats, saturation_stats, reopen_stats, profile_report, outdir,
        n_files, n_items, data_bytes, build_seconds)

    rm(tmp; force=true, recursive=true)   # synthetic data + cache gone; results kept
    tee_println("\nResults kept in: $outdir")
    tee_println("Log kept in: $log_path")
    return outdir
    finally
        RUN_LOG[] = nothing
        close(log_io)
    end
end

# --------------------------------------------------------------------------------------------------
# Warm reopen
# --------------------------------------------------------------------------------------------------

"""One close-and-reopen on the warm cache: time to first view and idle, allocation, first plots."""
function _reopen_once(project, data_root, plot_kinds, kinds)
    GC.gc()
    t0 = time()
    bytes0 = Base.gc_bytes()
    ws = MB.open_workspace(project, data_root)
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
        if !ENGINE_ONLY
            for kind in kinds
                probe = timed_plot!(ws, plot_kinds, kind, 1)
                probe === nothing || push!(first_plots,
                    (kind=kind, plot_ms=probe.plot_ms))
            end
        end
        return (first_view_s=first_view_s, idle_s=idle_s, alloc_bytes=alloc_bytes,
            items=length(ws.index.items), unchanged=ws.scan.state === :unchanged,
            first_plots=first_plots)
    finally
        MB.close_workspace!(ws)
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
function measure_reopen(project, data_root, plot_kinds, kinds)
    _reopen_once(project, data_root, plot_kinds, kinds)   # warm up JIT, discard
    return _reopen_once(project, data_root, plot_kinds, kinds)
end

# --------------------------------------------------------------------------------------------------
# Reporting
# --------------------------------------------------------------------------------------------------

_stat(v, f) = isempty(v) ? NaN : f(v)

function report(samples, stats, saturation, reopen, profile_report, outdir,
    n_files, n_items, data_bytes, build_seconds)
    saturation === nothing && error("Missing processed-writer saturation sample")
    # responsiveness CSV
    open(joinpath(outdir, "responsiveness.csv"), "w") do io
        println(io, "elapsed_s,phase,kind,n_items,plot_ms,allocated_bytes,ready_items")
        for s in samples
            @printf(io, "%.3f,%s,%s,%d,%.3f,%d,%d\n",
                s.elapsed_s, s.phase, s.kind, s.n, s.plot_ms,
                s.allocated_bytes, s.ready)
        end
    end

    open(joinpath(outdir, "memory_samples.csv"), "w") do io
        println(io, "elapsed_s,rss_bytes,gc_live_bytes,gc_allocated_bytes," *
            "rss_minus_gc_live_bytes,index_items,index_collections,item_metadata,analysis_errors," *
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
                sample.item_metadata,
                sample.analysis_errors,
                sample.processing_jobs,
                sample.pending_writes,
                sample.pending_write_rows,
                sample.selected_queue,
                sample.background_waiting)
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

    if profile_report isa MB.Profiling.ProfileReport
        event_groups = _profile_event_groups(profile_report)
        open(joinpath(outdir, "profile_summary.csv"), "w") do io
            println(io, "category,operation,count,total_ms,p50_ms,p90_ms,p99_ms,max_ms," *
                "wait_ms,service_ms,mean_batch,p50_batch,p90_batch,max_batch," *
                "mean_rows,p50_rows,p90_rows,max_rows")
            for row in profile_report.summary
                events = event_groups[(row.category, row.operation)]
                row_counts = Int64[event.attributes.rows for event in events
                                   if event.attributes.rows > 0]
                mean_rows = isempty(row_counts) ? 0.0 : sum(row_counts) / length(row_counts)
                p50_rows = isempty(row_counts) ? 0.0 : median(row_counts)
                p90_rows = isempty(row_counts) ? 0.0 : quantile(row_counts, 0.9)
                max_rows = isempty(row_counts) ? 0.0 : maximum(row_counts)
                @printf(io, "%s,%s,%d,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f\n",
                    row.category, row.operation, row.count, row.total_ms, row.median_ms,
                    row.p90_ms, row.p99_ms, row.max_ms, row.wait_ms, row.service_ms,
                    row.mean_batch, row.median_batch, row.p90_batch, row.max_batch,
                    mean_rows, p50_rows, p90_rows, max_rows)
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
    peak_rss_kib_per_item = stats.rss_peak_bytes / 1024 / indexed_items
    peak_gc_live_kib_per_item = maximum(
        (sample.gc_live_bytes for sample in stats.memory_samples);
        init=0,
    ) / 1024 / indexed_items
    open(joinpath(outdir, "scorecard.csv"), "w") do io
        println(io, "source_files,items,estimated_payload_rows,data_mib," *
            "rows_per_file,rows_per_item,build_s,scan_s,processing_s,analysis_s," *
            "scan_files_per_s,processing_items_per_s,build_items_per_s," *
            "build_ms_per_file,build_ms_per_item,scan_ms_per_file,processing_ms_per_item," *
            "write_ms_per_call,write_ms_per_file,write_ms_per_item,write_ns_per_payload_row," *
            "saturation_items,saturation_rows,saturation_load_ms,saturation_flush_ms," *
            "saturation_peak_pending_rows,saturation_processed_writes," *
            "during_plot_median_ms,during_plot_p90_ms,during_plot_p99_ms,during_plot_max_ms," *
            "after_plot_median_ms,after_plot_p90_ms,after_plot_p99_ms,after_plot_max_ms," *
            "rss_start_mib,rss_peak_mib,rss_end_mib,peak_rss_kib_per_item," *
            "peak_gc_live_kib_per_item,profile_events,profile_dropped")
        event_count = profile_report isa MB.Profiling.ProfileReport ?
            length(profile_report.events) : 0
        dropped = profile_report isa MB.Profiling.ProfileReport ?
            profile_report.dropped_events : 0
        @printf(io, "%d,%d,%d,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%d,%d,%.3f,%.3f,%d,%d,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.1f,%.1f,%.1f,%.3f,%.3f,%d,%d\n",
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
            read_stat(after, values -> quantile(values, 0.99)), read_stat(after, maximum),
            stats.rss_start_bytes / 1024^2,
            stats.rss_peak_bytes / 1024^2,
            stats.rss_end_bytes / 1024^2,
            peak_rss_kib_per_item,
            peak_gc_live_kib_per_item,
            event_count, dropped)
    end

    tee_println("\n==================== REALISTIC BROWSE BENCHMARK ====================")
    tee_printf("dataset:  %d files · ~%d items · %d estimated rows · %.1f MB\n",
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
    tee_printf("  memory              %8.1f KiB RSS/item  %8.1f KiB GC-live/item\n",
        peak_rss_kib_per_item, peak_gc_live_kib_per_item)

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

    tee_println("\nProcess memory:")
    tee_printf("  RSS start %.1f MiB  peak %.1f MiB  end %.1f MiB\n",
        stats.rss_start_bytes / 1024^2,
        stats.rss_peak_bytes / 1024^2,
        stats.rss_end_bytes / 1024^2)
    if !isempty(stats.memory_samples)
        peak_sample = stats.memory_samples[argmax(
            [sample.rss_bytes for sample in stats.memory_samples])]
        tee_printf("  peak sample: GC live %.1f MiB  RSS-GC-live %.1f MiB\n",
            peak_sample.gc_live_bytes / 1024^2,
            peak_sample.rss_minus_gc_live_bytes / 1024^2)
        tee_printf("  index counts: items %d  metadata %d  collections %d  errors %d\n",
            peak_sample.index_items,
            peak_sample.item_metadata,
            peak_sample.index_collections,
            peak_sample.analysis_errors)
        tee_printf("  queue counts: jobs %d  pending writes %d (%d rows)  selected %d  background waiting %d\n",
            peak_sample.processing_jobs,
            peak_sample.pending_writes,
            peak_sample.pending_write_rows,
            peak_sample.selected_queue,
            peak_sample.background_waiting)
    end

    tee_println("\nWarm reopen (same cache, every fingerprint unchanged):")
    tee_printf("  first cached view %.0f ms  ·  idle %.2f s  ·  allocated %.1f MiB  ·  %d items%s\n",
        reopen.first_view_s * 1e3, reopen.idle_s, reopen.alloc_bytes / 1024^2, reopen.items,
        reopen.unchanged ? "  (reused)" : "  (re-scanned!)")
    for p in reopen.first_plots
        tee_printf("  first %-6s plot  %6.1f ms\n", p.kind, p.plot_ms)
    end

    if profile_report isa MB.Profiling.ProfileReport
        tee_printf("\nStructured profile: %d events, %d counters, %d dropped, %d CPU samples\n",
            length(profile_report.events), length(profile_report.counters),
            profile_report.dropped_events,
            profile_report.cpu === nothing ? 0 : profile_report.cpu.total_samples)
    end
    _maybe_plot_pipeline(profile_report, stats.memory_samples, outdir)
    return nothing
end

function _metric_values(rows::Vector{NamedTuple}, metric::String)::Vector{Float64}
    return [row.ms for row in rows if row.metric == metric]
end

function _plot_time_hist!(CM, fig, position, rows, metric::String, title::String)::Nothing
    values = _metric_values(rows, metric)
    ax = CM.Axis(fig[position...]; xlabel="milliseconds", ylabel="events",
        title="$title (n=$(length(values)))")
    isempty(values) || CM.hist!(ax, values; bins=min(80, max(10, ceil(Int, sqrt(length(values))))))
    return nothing
end

"""Render pipeline timing plots and write the plotted timing tables."""
function _maybe_plot_pipeline(profile_report, memory_samples::Vector{MemorySample}, outdir::String)::Nothing
    if !(profile_report isa MB.Profiling.ProfileReport)
        tee_println("\nNo structured profile captured; skipping pipeline event plots.")
        _write_pipeline_timeseries!(memory_samples, outdir)
        return nothing
    end
    event_rows = _write_pipeline_event_times!(profile_report, outdir)
    _write_pipeline_timeseries!(memory_samples, outdir)

    tee_println("\nPipeline event datapoints:")
    if isempty(event_rows)
        tee_println("  no event timing rows captured")
    else
        summary = combine(groupby(DataFrame(event_rows), :metric), nrow => :count)
        sort!(summary, :metric)
        for row in eachrow(summary)
            tee_printf("  %-34s %6d\n", row.metric, row.count)
        end
    end

    if !HAS_CAIRO
        tee_println("\n(No CairoMakie in the bench env; wrote pipeline CSVs only.)")
        return nothing
    end

    CM = CairoMakie
    fig = CM.Figure(size=(1500, 1200), fontsize=13)
    CM.Label(fig[0, 1:3], "Data Pipeline Timing", fontsize=20, font=:bold)
    _plot_time_hist!(CM, fig, (1, 1), event_rows,
        "processing_queue_wait_ms", "Processing queue wait")
    _plot_time_hist!(CM, fig, (1, 2), event_rows,
        "processing_engine_overhead_ms", "Processing engine overhead")
    _plot_time_hist!(CM, fig, (1, 3), event_rows,
        "interpret_engine_overhead_ms", "Interpretation engine overhead")
    _plot_time_hist!(CM, fig, (2, 1), event_rows,
        "cache_writer_ms", "Cache writer")
    _plot_time_hist!(CM, fig, (2, 2), event_rows,
        "cache_write_dataframe_ms", "DataFrame cache writes")
    _plot_time_hist!(CM, fig, (2, 3), event_rows,
        "cache_read_item_data_ms", "Cache item reads")

    elapsed = [sample.elapsed_s for sample in memory_samples]
    pending_rows = [sample.pending_write_rows for sample in memory_samples]
    rss_mib = [sample.rss_bytes / 1024^2 for sample in memory_samples]
    gc_live_mib = [sample.gc_live_bytes / 1024^2 for sample in memory_samples]

    ax_rows = CM.Axis(fig[3, 1:3]; xlabel="elapsed seconds", ylabel="rows",
        title="Pending write rows over time")
    isempty(elapsed) || CM.lines!(ax_rows, elapsed, pending_rows)

    ax_mem = CM.Axis(fig[4, 1:3]; xlabel="elapsed seconds", ylabel="MiB",
        title="RSS and GC-live memory over time")
    if !isempty(elapsed)
        CM.lines!(ax_mem, elapsed, rss_mib; label="RSS")
        CM.lines!(ax_mem, elapsed, gc_live_mib; label="GC live")
        CM.axislegend(ax_mem; position=:lt)
    end

    path = joinpath(outdir, "pipeline.png")
    CM.save(path, fig)
    tee_println("\npipeline plot: $path")
    return nothing
end

for repeat in 1:BENCH_REPEATS
    BENCH_REPEATS == 1 || println("\n===== benchmark repeat $repeat / $BENCH_REPEATS =====\n")
    run_benchmark()
end
