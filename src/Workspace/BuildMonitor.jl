# Live build instrumentation. A real cache build runs for minutes, so the engine reports rolling
# timings while it works rather than only at the end. It covers the cache-I/O phases the project
# scan profile (read/entries/process/stats) omits — interpreted, processed, and stats writes — which
# run on background workers, hence the atomic counters.

mutable struct BuildMonitor
    io::IO
    interval::Float64
    interpreted_write_ns::Base.Threads.Atomic{Int64}
    interpreted_writes::Base.Threads.Atomic{Int64}
    processed_write_ns::Base.Threads.Atomic{Int64}
    processed_writes::Base.Threads.Atomic{Int64}
    stats_write_ns::Base.Threads.Atomic{Int64}
    stats_writes::Base.Threads.Atomic{Int64}
    active::Bool
    operation::Symbol
    started_at::Float64
    last_report_at::Float64
    last_scan_done::Int
    last_processed_done::Int
    peak_rss::Int64
end

function BuildMonitor(; io::IO=stderr, interval::Real=3.0)::BuildMonitor
    zero_counter() = Base.Threads.Atomic{Int64}(0)
    return BuildMonitor(
        io, Float64(interval),
        zero_counter(), zero_counter(), zero_counter(),
        zero_counter(), zero_counter(), zero_counter(),
        false, :idle, 0.0, 0.0, 0, 0, 0)
end

"""Add one timed cache-I/O call to a monitor counter pair (nanoseconds since `t0_ns`)."""
@inline function record_cache_phase!(
    total_ns::Base.Threads.Atomic{Int64},
    calls::Base.Threads.Atomic{Int64},
    t0_ns::UInt64,
)::Nothing
    Base.Threads.atomic_add!(total_ns, Int64(time_ns() - t0_ns))
    Base.Threads.atomic_add!(calls, Int64(1))
    return nothing
end

"""Reset the counters and arm the monitor at the start of a scan/build."""
function begin_build_monitor!(monitor::BuildMonitor, operation::Symbol)::Nothing
    for counter in (
        monitor.interpreted_write_ns, monitor.interpreted_writes,
        monitor.processed_write_ns, monitor.processed_writes,
        monitor.stats_write_ns, monitor.stats_writes,
    )
        counter[] = 0
    end
    now = time()
    monitor.active = true
    monitor.operation = operation
    monitor.started_at = now
    monitor.last_report_at = now
    monitor.last_scan_done = 0
    monitor.last_processed_done = 0
    monitor.peak_rss = 0
    return nothing
end

function _fmt_eta(seconds::Real)::String
    (!isfinite(seconds) || seconds < 0) && return "  —  "
    hours, rest = divrem(round(Int, seconds), 3600)
    minutes, secs = divrem(rest, 60)
    return hours > 0 ? @sprintf("%d:%02d:%02d", hours, minutes, secs) :
           @sprintf("%02d:%02d", minutes, secs)
end

function _fmt_bytes(bytes::Real)::String
    units = ("B", "KiB", "MiB", "GiB", "TiB")
    value, index = Float64(bytes), 1
    while value >= 1024 && index < length(units)
        value /= 1024
        index += 1
    end
    return @sprintf("%.1f %s", value, units[index])
end

_avg_ms(total_ns::Base.Threads.Atomic{Int64}, calls::Base.Threads.Atomic{Int64})::Float64 =
    calls[] == 0 ? 0.0 : total_ns[] / calls[] / 1e6

"""Print a rolling progress group once per interval. Cheap and silent otherwise; safe every poll."""
function maybe_report_build!(workspace)::Nothing
    monitor = workspace.monitor
    monitor.active || return nothing
    now = time()
    (now - monitor.last_report_at) < monitor.interval && return nothing

    progress = workspace.scan.progress
    total_si = progress.total_source_items
    done_si = progress.processed_source_items
    completed, queued = lock(workspace.processing.lock) do
        (workspace.processing.completed, workspace.processing.total)
    end
    elapsed = now - monitor.started_at
    interval = now - monitor.last_report_at
    monitor.peak_rss = max(monitor.peak_rss, Int64(Sys.maxrss()))

    # Instantaneous rate shows the current pace (and stalls); the ETA uses the cumulative average so a
    # momentary stall does not throw it to hours.
    scan_rate = (done_si - monitor.last_scan_done) / interval
    process_rate = (completed - monitor.last_processed_done) / interval
    scan_eta = done_si > 0 ? (total_si - done_si) * elapsed / done_si : Inf
    process_eta = completed > 0 ? (queued - completed) * elapsed / completed : Inf
    scan_pct = total_si > 0 ? 100 * done_si / total_si : 0.0
    process_pct = queued > 0 ? 100 * completed / queued : 0.0

    io = monitor.io
    println(io, "── build $(_fmt_eta(elapsed)) elapsed  ·  $(monitor.operation)  ·  RSS peak $(_fmt_bytes(monitor.peak_rss)) ──")
    @printf(io, "   scan     %5d/%-5d (%4.1f%%)  %6.1f src/s  ETA %s%s\n",
        done_si, total_si, scan_pct, scan_rate, _fmt_eta(scan_eta),
        isempty(progress.current_source_item) ? "" : "  ← " * basename(progress.current_source_item))
    @printf(io, "   analysis %5d/%-5d (%4.1f%%)  %6.1f item/s ETA %s\n",
        completed, queued, process_pct, process_rate, _fmt_eta(process_eta))

    rows = scan_profile_summary(workspace.project)
    if !isempty(rows)
        # read/entries run once per source file (scan); process/stats run once per data item (analysis).
        @printf(io, "   %-13s %6s %8s %9s %9s %8s\n",
            "kind", "items", "read/f", "entries/f", "process/i", "stats/i")
        for row in Iterators.take(rows, 6)
            files = max(row.source_items, 1)
            items = max(row.items, 1)
            @printf(io, "   %-13s %6d %6.1fms %7.1fms %7.1fms %6.1fms\n",
                String(row.kind), row.items,
                1e3 * row.read_seconds / files, 1e3 * row.entries_seconds / files,
                1e3 * row.process_seconds / items, 1e3 * row.stats_seconds / items)
        end
    end

    # The single writer is shared by the scan (interpreted) and analysis (processed/stats). Busy is the
    # real serialized write work; wait is how long threads blocked for their turn. Busy near elapsed
    # with large wait means the one writer is the bottleneck.
    busy_s = workspace.cache.db.writer_busy_ns[] / 1e9
    wait_s = workspace.cache.db.writer_wait_ns[] / 1e9
    @printf(io, "   writer: busy %.1fs (%.0f%% of elapsed) · threads waited %.0fs total · interp %.1fs processed %.1fs stats %.1fs\n",
        busy_s, elapsed > 0 ? 100 * busy_s / elapsed : 0.0, wait_s,
        monitor.interpreted_write_ns[] / 1e9, monitor.processed_write_ns[] / 1e9,
        monitor.stats_write_ns[] / 1e9)
    flush(io)

    monitor.last_report_at = now
    monitor.last_scan_done = done_si
    monitor.last_processed_done = completed
    return nothing
end

"""Print a final one-line summary and disarm. Silent for builds too short to have reported."""
function finish_build_monitor!(workspace)::Nothing
    monitor = workspace.monitor
    monitor.active || return nothing
    monitor.active = false
    elapsed = time() - monitor.started_at
    elapsed < monitor.interval && return nothing
    progress = workspace.scan.progress
    monitor.peak_rss = max(monitor.peak_rss, Int64(Sys.maxrss()))
    rate = elapsed > 0 ? progress.total_source_items / elapsed : 0.0
    busy_s = workspace.cache.db.writer_busy_ns[] / 1e9
    println(monitor.io,
        "✓ build done in $(_fmt_eta(elapsed))  ·  $(progress.total_source_items) source items, " *
        "$(progress.loaded_items) data items ($(@sprintf("%.1f", rate)) src/s)  ·  writer busy " *
        "$(@sprintf("%.1f", busy_s))s ($(elapsed > 0 ? round(Int, 100 * busy_s / elapsed) : 0)% of elapsed) " *
        "· RSS peak $(_fmt_bytes(monitor.peak_rss))")
    flush(monitor.io)
    return nothing
end
