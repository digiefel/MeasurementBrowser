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
    # When recording, one CSV row is appended per sample tick. The file is the deliverable — load it in
    # Julia/pandas/a spreadsheet and plot write-time vs. items to *see* the O(N²), instead of reading spew.
    csv_io::Union{Nothing,IO}
    csv_path::String
    event_log::Bool
end

function BuildMonitor(; io::IO=stderr, interval::Real=1.0)::BuildMonitor
    zero_counter() = Base.Threads.Atomic{Int64}(0)
    return BuildMonitor(
        io, Float64(interval),
        zero_counter(), zero_counter(), zero_counter(),
        zero_counter(), zero_counter(), zero_counter(),
        false, :idle, 0.0, 0.0, 0, 0, 0, nothing, "", false)
end

# Per-build sequence so concurrent/successive builds never write to the same log file.
const BUILD_LOG_SEQ = Base.Threads.Atomic{Int}(0)

"""Insert `suffix` before a path's extension: `/tmp/build.csv` + `-17…-1` → `/tmp/build-17…-1.csv`."""
function _with_log_suffix(base::AbstractString, suffix::AbstractString)::String
    stem, ext = splitext(String(base))
    return string(stem, suffix, isempty(ext) ? ".csv" : ext)
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

"""
Reset the counters and arm the monitor at the start of a scan/build.

`build_profile` (a path or `nothing`) records the per-second CSV of scan/analysis/writer metrics;
`event_log` (a path or `nothing`) records the crash-durable per-operation log. Both come from the
workspace's open_workspace arguments; when both are `nothing` the build runs silent.
"""
function begin_build_monitor!(
    monitor::BuildMonitor,
    operation::Symbol;
    event_log::Union{Nothing,AbstractString}=nothing,
    build_profile::Union{Nothing,AbstractString}=nothing,
)::Nothing
    for counter in (
        monitor.interpreted_write_ns, monitor.interpreted_writes,
        monitor.processed_write_ns, monitor.processed_writes,
        monitor.stats_write_ns, monitor.stats_writes,
    )
        counter[] = 0
    end
    now = time()
    # The first build of a session writes to the exact path you passed; later builds get -2, -3, … so a
    # rescan never clobbers the run you wanted to keep. Predictable and findable — no random names. The
    # actual path is printed on finish.
    seq = Base.Threads.atomic_add!(BUILD_LOG_SEQ, 1) + 1
    suffix = seq == 1 ? "" : "-$seq"
    monitor.csv_io === nothing || (try; close(monitor.csv_io); catch; end)
    monitor.csv_io = nothing
    monitor.csv_path = ""
    if build_profile !== nothing
        path = _with_log_suffix(build_profile, suffix)
        try
            io = open(path, "w")
            println(io,
                "elapsed_s,operation,scan_done,scan_total,analysis_done,analysis_total," *
                "interp_write_s,processed_write_s,stats_write_s,writer_busy_s,writer_wait_s,rss_bytes")
            flush(io)
            monitor.csv_io = io
            monitor.csv_path = path
        catch err
            @warn "Could not open the build-profile CSV; recording disabled" path exception = err
        end
    end
    monitor.event_log = event_log !== nothing
    monitor.event_log && start_event_log!(_with_log_suffix(event_log, suffix))
    monitor.active = monitor.csv_io !== nothing || monitor.event_log
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

"""Append one CSV sample row: a cheap, lock-light snapshot taken from the poll thread."""
function _write_build_csv_row!(workspace, now::Float64)::Nothing
    monitor = workspace.monitor
    io = monitor.csv_io
    io === nothing && return nothing
    progress = workspace.scan.progress
    completed, queued = lock(workspace.processing.lock) do
        (workspace.processing.completed, workspace.processing.total)
    end
    rss = Int64(Sys.maxrss())
    monitor.peak_rss = max(monitor.peak_rss, rss)
    db = workspace.cache.db
    @printf(io, "%.3f,%s,%d,%d,%d,%d,%.6f,%.6f,%.6f,%.6f,%.6f,%d\n",
        now - monitor.started_at, String(monitor.operation),
        progress.processed_source_items, progress.total_source_items,
        completed, queued,
        monitor.interpreted_write_ns[] / 1e9, monitor.processed_write_ns[] / 1e9,
        monitor.stats_write_ns[] / 1e9, db.writer_busy_ns[] / 1e9, db.writer_wait_ns[] / 1e9, rss)
    flush(io)
    return nothing
end

"""Record one CSV sample row once per interval. Cheap and silent otherwise; safe every poll."""
function maybe_report_build!(workspace)::Nothing
    monitor = workspace.monitor
    monitor.active || return nothing
    now = time()
    (now - monitor.last_report_at) < monitor.interval && return nothing
    _write_build_csv_row!(workspace, now)
    monitor.last_report_at = now
    return nothing
end

"""Write the final CSV row, close the file, and point at it on the console. Then disarm."""
function finish_build_monitor!(workspace)::Nothing
    monitor = workspace.monitor
    monitor.active || return nothing
    monitor.active = false
    if monitor.event_log
        path = current_event_log_path()
        stop_event_log!()
        monitor.event_log = false
        println(monitor.io, "✓ event log written to $path")
    end
    _write_build_csv_row!(workspace, time())
    io = monitor.csv_io
    monitor.csv_io = nothing
    io === nothing && (flush(monitor.io); return nothing)
    try; close(io); catch; end
    println(monitor.io, "✓ build profile written to $(monitor.csv_path)")
    flush(monitor.io)
    return nothing
end
