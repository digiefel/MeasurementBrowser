# Combined performance run. Writes bench/status.txt only on success, and otherwise clears the file.
#   julia --project=bench --threads=auto bench/run.jl     # this file only
#   bench/run.sh                                          # this file + peak RSS at 0.1 scale, prints status.txt
# test/runtests.jl includes this after the unit tests (stdout discarded).
# Tunables: MB_BENCH_SCALE (default 0.05 here, 0.1 in run.sh), MB_BENCH_SCALING_SIZES (250,500,1000).

using Printf
using Statistics: median, var

const STATUS_PATH = joinpath(@__DIR__, "status.txt")

_fmt(x::Integer) = string(x)
_fmt(x::AbstractFloat) = @sprintf("%.4g", x)
_fmt(x::AbstractString) = x

function _line(io, name, value, comment)
    @printf(io, "%-40s  %-24s  # %s\n", name, _fmt(value), comment)
    return nothing
end

function _line(io, name, values::AbstractVector, comment)
    val = join((_fmt(v) for v in values), "  ")
    @printf(io, "%-40s  %-24s  # %s\n", name, val, comment)
    return nothing
end

function _sample(io, name, values::Vector{<:Real}, comment)
    isempty(values) && return _line(io, name, NaN, comment)
    length(values) == 1 && return _line(io, name, Float64(values[1]), comment)
    return _line(io, name, [median(values), var(values)], comment)
end

function write_status(realistic, scaling; startup_s)
    stats = realistic.stats
    saturation = realistic.saturation
    reopen = realistic.reopen
    saturation === nothing && error("Missing processed-writer saturation sample")
    stats === nothing && error("Missing build stats")
    if realistic.require_saturation
        stats.interpreted_writes > 0 || error("Benchmark did not exercise interpreted writes")
        stats.processed_writes > 0 || error("Benchmark did not exercise processed writes")
        stats.metadata_writes > 0 || error("Benchmark did not exercise metadata writes")
        saturation.estimated_rows >= realistic.cache_row_ceiling || error(
            "Processed-writer saturation selected only $(saturation.estimated_rows) rows",
        )
        saturation.processed_writes > 0 || error(
            "Processed-writer saturation created no processed writes",
        )
    end
    during = [s.plot_ms for s in realistic.samples if s.phase === :during_build]
    after = [s.plot_ms for s in realistic.samples if s.phase === :after_build]
    open(STATUS_PATH, "w") do io
        println(io, "# inputs")
        _line(io, "scale", realistic.scale,
            "Multiplier on how many aliases of each template file are advertised. Does not change rows per file.")
        _line(io, "kind1_files", realistic.kind1_files,
            "Kind1 aliases (tiny IV-style tables, one item per file).")
        _line(io, "kind2_files", realistic.kind2_files,
            "Kind2 aliases (medium CV-style tables, one item per file).")
        _line(io, "kind3_files", realistic.kind3_files,
            "Kind3 aliases (large fatigue-style tables; one item per cycle).")
        _line(io, "kind1_rows", realistic.kind1_rows,
            "Row count of bench/templates/kind1.csv.")
        _line(io, "kind2_rows", realistic.kind2_rows,
            "Row count of bench/templates/kind2.csv.")
        _line(io, "kind3_cycles", realistic.kind3_cycles,
            "Distinct cycle values in bench/templates/kind3.csv. Each cycle becomes one item.")
        _line(io, "kind3_rows", realistic.kind3_rows,
            "Row count per cycle in bench/templates/kind3.csv.")
        _line(io, "processed_stress_rows", realistic.processed_stress_rows,
            "Requested row budget for the processed-payload flush probe. The probe takes min(available kind3 items, this budget / kind3_rows).")
        _line(io, "after_build_plots", realistic.after_build_plots,
            "How many plot probes to run per kind after the workspace is idle.")
        _line(io, "scaling_sizes", scaling.sizes,
            "Item counts N for the scaling sweep. Each scaling_*_ms row has one value per size, in this order.")
        _line(io, "threads", Threads.nthreads(),
            "Julia thread count (Threads.nthreads()).")
        _line(io, "julia", string(VERSION),
            "Julia version.")
        println(io)
        println(io, "# outputs")
        _line(io, "startup_s", startup_s,
            "Wall seconds to load packages and include the two harness scripts. On a warm pkgimage cache this is load time; if packages were invalidated it includes precompile.")
        _line(io, "aliases", realistic.n_files,
            "Source items the bench source advertised (kind1_files + kind2_files + kind3_files).")
        _line(io, "items", stats.processed_items,
            "Data items in the index when the cold build finished. Kind3 contributes one item per cycle.")
        _line(io, "build_s", realistic.build_seconds,
            "Wall seconds from open_workspace until the workspace is idle (scan done, workers idle, no pending cache writes). Plot probes run during this wait, so first-plot compilation is inside this number.")
        _line(io, "scan_s", stats.scan_seconds,
            "Wall seconds from open until the source scan reports done, unchanged, error, or canceled.")
        _line(io, "processing_s", stats.processing_seconds,
            "Wall seconds during which item process/analyze workers were busy or the cache still had pending writes. This is a busy-window, not CPU time of process().")
        _line(io, "analysis_s", stats.analysis_seconds,
            "Wall seconds during which collection-analyze workers were busy. Failed collection jobs still count as busy until they settle.")
        _line(io, "interpreted_writes", stats.interpreted_writes,
            "Interpreted-payload cache writes during the cold build.")
        _line(io, "processed_writes", stats.processed_writes,
            "Processed-payload cache writes during the cold build.")
        _line(io, "metadata_writes", stats.metadata_writes,
            "Metadata cache writes during the cold build.")
        println(io, "#                                        median                    variance")
        _sample(io, "plot_during_ms", during,
            "Plot round-trip (select + materialize + setup + draw) sampled every ~160 ms while the cold build is still running, in milliseconds. The first samples include Julia compiling the plot path, which inflates variance.")
        _sample(io, "plot_after_ms", after,
            "The same plot round-trip after the workspace is idle, in milliseconds.")
        _line(io, "saturation_peak_pending_rows", saturation.peak_pending_rows,
            "High-water mark of rows sitting in the cache write buffer while flushing the processed-payload stress selection.")
        _line(io, "saturation_processed_writes", saturation.processed_writes,
            "Processed-payload writes created by that stress selection. Zero means the payloads were already in cache, so this was a read rather than a writer-saturation probe.")
        _line(io, "saturation_load_ms", saturation.load_ms,
            "Milliseconds to materialize the stress selection.")
        _line(io, "saturation_flush_ms", saturation.flush_ms,
            "Milliseconds to wait until the cache write buffer is empty after that materialize.")
        _line(io, "reopen_first_view_s", reopen.first_view_s,
            "Seconds from open until the first item appears in the index on the second close-and-reopen. The first reopen is discarded as JIT warmup.")
        _line(io, "reopen_idle_s", reopen.idle_s,
            "Seconds from open until the workspace is idle on that second reopen.")
        _line(io, "reopen_alloc_mib", reopen.alloc_bytes / 1024^2,
            "MiB allocated according to Base.gc_bytes() during that second reopen. This is Julia GC accounting, not process RSS (see peak_rss_kb).")
        _line(io, "reopen_unchanged", reopen.unchanged ? 1 : 0,
            "1 if that second open reused the cache (scan state :unchanged), 0 if it rebuilt.")
        for p in reopen.first_plots
            _line(io, "reopen_plot_$(p.kind)_ms", p.plot_ms,
                "Milliseconds for one plot round-trip of $(p.kind) after the second reopen.")
        end
        _line(io, "scaling_scan_build_per_item_exponent", scaling.fits["scan_build_per_item"].exponent,
            "Slope of log(scan_build_per_item time) versus log(N). 0 means the cost does not grow with item count, 1 means linear in N, 2 means quadratic.")
        _line(io, "scaling_scan_build_per_item_r2", scaling.fits["scan_build_per_item"].r2,
            "R² of the log(scan_build_per_item time) versus log(N) fit. 1 means the times lie on a clean power law; a low value means the exponent is not a trustworthy summary.")
        _line(io, "scaling_scan_build_per_item_ms", scaling.ms_by_op["scan_build_per_item"],
            "Milliseconds of cold open_workspace + wait-until-idle, divided by N, at each scaling_sizes. One shot per N, not a median of repeats.")
        _line(io, "scaling_status_refresh_exponent", scaling.fits["status_refresh"].exponent,
            "Slope of log(status_refresh time) versus log(N). 0 means the cost does not grow with item count, 1 means linear in N, 2 means quadratic.")
        _line(io, "scaling_status_refresh_r2", scaling.fits["status_refresh"].r2,
            "R² of the log(status_refresh time) versus log(N) fit. 1 means the times lie on a clean power law; a low value means the exponent is not a trustworthy summary.")
        _line(io, "scaling_status_refresh_ms", scaling.ms_by_op["status_refresh"],
            "Median milliseconds per refresh_status! call (the per-frame GUI status rebuild) at each scaling_sizes.")
        _line(io, "scaling_workspace_busy_exponent", scaling.fits["workspace_busy"].exponent,
            "Slope of log(workspace_busy time) versus log(N). 0 means the cost does not grow with item count, 1 means linear in N, 2 means quadratic.")
        _line(io, "scaling_workspace_busy_r2", scaling.fits["workspace_busy"].r2,
            "R² of the log(workspace_busy time) versus log(N) fit. 1 means the times lie on a clean power law; a low value means the exponent is not a trustworthy summary.")
        _line(io, "scaling_workspace_busy_ms", scaling.ms_by_op["workspace_busy"],
            "Median milliseconds per workspace_busy call at each scaling_sizes.")
        _line(io, "scaling_items_panel_exponent", scaling.fits["items_panel"].exponent,
            "Slope of log(items_panel time) versus log(N). 0 means the cost does not grow with item count, 1 means linear in N, 2 means quadratic.")
        _line(io, "scaling_items_panel_r2", scaling.fits["items_panel"].r2,
            "R² of the log(items_panel time) versus log(N) fit. 1 means the times lie on a clean power law; a low value means the exponent is not a trustworthy summary.")
        _line(io, "scaling_items_panel_ms", scaling.ms_by_op["items_panel"],
            "Median milliseconds to gather and sort the items of the selected collection (what the items panel does each frame) at each scaling_sizes.")
        _line(io, "scaling_metadata_publish_exponent", scaling.fits["metadata_publish"].exponent,
            "Slope of log(metadata_publish time) versus log(N). 0 means the cost does not grow with item count, 1 means linear in N, 2 means quadratic.")
        _line(io, "scaling_metadata_publish_r2", scaling.fits["metadata_publish"].r2,
            "R² of the log(metadata_publish time) versus log(N) fit. 1 means the times lie on a clean power law; a low value means the exponent is not a trustworthy summary.")
        _line(io, "scaling_metadata_publish_ms", scaling.ms_by_op["metadata_publish"],
            "Median milliseconds for one reconcile_source_metadata_cache!(refresh_hierarchy=true) at each scaling_sizes.")
    end
    return nothing
end

empty!(ARGS)
get!(ENV, "MB_BENCH_SCALE", "0.05")
get!(ENV, "MB_BENCH_SCALING_SIZES", "250,500,1000")

const STARTUP_S = @elapsed begin
    include(joinpath(@__DIR__, "realistic_browse.jl"))
    include(joinpath(@__DIR__, "scaling.jl"))
end

function run_performance()
    println("performance snapshot: realistic browse, then scaling")
    write(STATUS_PATH, "")
    try
        realistic = run_benchmark()
        scaling = run_scaling()
        write_status(realistic, scaling; startup_s=STARTUP_S)
        return nothing
    catch
        write(STATUS_PATH, "")
        rethrow()
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_performance()
end
