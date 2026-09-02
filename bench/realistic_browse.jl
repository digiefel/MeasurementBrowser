# Realistic browse-while-building benchmark.
#
#   julia --project=bench --threads=auto bench/realistic_browse.jl [scale]
#
# `scale` (default 0.05, or MB_BENCH_SCALE) multiplies alias counts, not rows per file.
# Standalone prints debug timings. The committed file is written by bench/run.jl.
#
# ENV: MB_BENCH_KIND1_FILES, MB_BENCH_KIND2_FILES, MB_BENCH_KIND3_FILES,
#      MB_BENCH_AFTER_BUILD_PLOTS, MB_BENCH_PROCESSED_STRESS_ROWS, MB_BENCH_SCALE

using DataBrowserAPI: item_data, label
using DataBrowserRecipes: define_project, register_item!
using DataBrowserAPI.ItemIndex: ItemRecord
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
using Random
using Statistics: mean
import GLMakie: Figure, Axis, lines!, contents

if !@isdefined(BenchSource)
    include(joinpath(@__DIR__, "custom_data_source.jl"))
end

# --------------------------------------------------------------------------------------------------
# Sizing (alias counts; row layout is the templates)
# --------------------------------------------------------------------------------------------------

scale = length(ARGS) >= 1 ? parse(Float64, ARGS[1]) :
    parse(Float64, get(ENV, "MB_BENCH_SCALE", "0.05"))
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
const MAX_BUILD_SECONDS = 600
const REQUIRE_SATURATION = scale >= 1.0

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

"""Write a taken `TimerOutput` to stdout."""
function show_debug_timings(timings)
    show(stdout, MIME("text/plain"), timings)
    println()
    return nothing
end

function run_benchmark()
    tmp = mktempdir()
    pushfirst!(DEPOT_PATH, tmp)
    source = _bench_source()
    n_files = length(source.files)
    n_items = KIND1_FILES + KIND2_FILES + KIND3_FILES * KIND3_CYCLES
    Profiling.reset_debug_timings!()

    try
        println("Templates: ", n_files, " aliases, ~", n_items, " items (scale=", scale, ")")

        project = build_project(; plots=true)
        kinds = (:kind1, :kind2, :kind3)
        plot_kinds = Dict(k => first(registered_plot_kinds(project, k)) for k in kinds)
        samples = Sample[]

        println("Building cache + browsing during the scan ...")
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

            println("Saturating processed-payload writer ...")
            saturation_stats = saturate_processed_writes!(ws, :kind3)

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
            metrics = ws.metrics
            build_stats = (
                scan_seconds,
                processing_seconds,
                analysis_seconds,
                processed_items=length(ws.index.items),
                interpreted_writes=metrics.interpreted_writes[],
                processed_writes=metrics.processed_writes[],
                metadata_writes=metrics.metadata_writes[],
            )
        finally
            close_workspace!(ws)
        end

        reopen_stats = measure_reopen(project, source, plot_kinds, kinds)
        return (
            samples=samples,
            stats=build_stats,
            saturation=saturation_stats,
            reopen=reopen_stats,
            n_files=n_files,
            build_seconds=build_seconds,
            scale=scale,
            kind1_files=KIND1_FILES,
            kind2_files=KIND2_FILES,
            kind3_files=KIND3_FILES,
            kind1_rows=KIND1_ROWS,
            kind2_rows=KIND2_ROWS,
            kind3_cycles=KIND3_CYCLES,
            kind3_rows=KIND3_ROWS,
            processed_stress_rows=PROCESSED_STRESS_ROWS,
            after_build_plots=AFTER_BUILD_PLOTS,
            require_saturation=REQUIRE_SATURATION,
            cache_row_ceiling=CACHE_ROW_CEILING,
        )
    finally
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

function measure_reopen(project, source, plot_kinds, kinds)
    _reopen_once(project, source, plot_kinds, kinds)
    return _reopen_once(project, source, plot_kinds, kinds)
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_benchmark()
    show_debug_timings(Profiling.take_debug_timings!())
end
