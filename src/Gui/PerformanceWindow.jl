using Printf
using Statistics: mean
import CImGui as ig
using GLMakie: Axis, Figure, Observable, axislegend, lines!
import GLMakie.Makie as Makie
using NativeFileDialog: save_file

import ..Profiling
using ..Projects: scan_profile_summary, scan_source_profile
import ..Workspace
using .MakieImguiIntegration: MakieFigure, makie_context_stats

# ---------------------------------------------------------------------------
# Internal table helpers
# ---------------------------------------------------------------------------

"""
Render a right-to-left aligned string inside the current table column.
Because CImGui.jl does not expose `SetNextItemWidth(-1)` plus right-align text natively
for arbitrary strings, we just render text and let column widths handle alignment.
"""
@inline function _table_text(s::AbstractString)::Nothing
    ig.TextUnformatted(s)
    return nothing
end

"""Begin a compact fixed-layout table with headers, row background, and inner borders."""
function _begin_perf_table(id::String, n_cols::Int, height::Float32=0.0f0)::Bool
    flags = ig.ImGuiTableFlags_RowBg         |
            ig.ImGuiTableFlags_BordersInnerH  |
            ig.ImGuiTableFlags_BordersOuterH  |
            ig.ImGuiTableFlags_BordersOuterV  |
            ig.ImGuiTableFlags_Resizable      |
            ig.ImGuiTableFlags_ScrollY        |
            ig.ImGuiTableFlags_SizingFixedFit
    return ig.BeginTable("##perf_$id", n_cols, flags, (0.0f0, height))
end

# ---------------------------------------------------------------------------
# Tab: Plot Redraw
# ---------------------------------------------------------------------------

"""Render the Plot Redraw tab: per-phase columns mean / max / last / alloc for each phase."""
function _render_plot_redraw_tab(performance::PerformanceState)::Nothing
    timings = performance.timings
    allocs  = performance.allocations

    # The four canonical plot-redraw phase keys, in order
    phases = (
        (:plot_load,  "Load"),
        (:plot_setup, "Setup"),
        (:plot_data,  "Draw data"),
        (:plot_draw,  "Total"),
    )

    has_any = any(haskey(timings, k) && !isempty(timings[k]) for (k, _) in phases)
    if !has_any
        ig.TextDisabled("No plot redraw recorded yet — open a plot to capture timings.")
        return nothing
    end

    if _begin_perf_table("plot_redraw", 6, 160.0f0)
        ig.TableSetupScrollFreeze(0, 1)
        ig.TableSetupColumn("Phase",    ig.ImGuiTableColumnFlags_WidthFixed, 90.0f0)
        ig.TableSetupColumn("n",        ig.ImGuiTableColumnFlags_WidthFixed, 40.0f0)
        ig.TableSetupColumn("Mean ms",  ig.ImGuiTableColumnFlags_WidthFixed, 72.0f0)
        ig.TableSetupColumn("Max ms",   ig.ImGuiTableColumnFlags_WidthFixed, 72.0f0)
        ig.TableSetupColumn("Last ms",  ig.ImGuiTableColumnFlags_WidthFixed, 72.0f0)
        ig.TableSetupColumn("Alloc KB", ig.ImGuiTableColumnFlags_WidthFixed, 72.0f0)
        ig.TableHeadersRow()

        for (key, label) in phases
            samples = get(timings, key, Float64[])
            isempty(samples) && continue
            bytes = get(allocs, key, Int[])

            ig.TableNextRow()
            ig.TableNextColumn(); _table_text(label)
            ig.TableNextColumn(); _table_text(string(length(samples)))
            ig.TableNextColumn(); _table_text(@sprintf("%.2f", mean(samples)))
            ig.TableNextColumn(); _table_text(@sprintf("%.2f", maximum(samples)))
            ig.TableNextColumn(); _table_text(@sprintf("%.2f", samples[end]))
            ig.TableNextColumn(); _table_text(
                isempty(bytes) ? "—" : @sprintf("%.1f", bytes[end] / 1024))
        end
        ig.EndTable()
    end
    return nothing
end

# ---------------------------------------------------------------------------
# Tab: Sampling Hotspots
# ---------------------------------------------------------------------------

"""Render CPU hotspots from the latest internal capture or explicit plot profile."""
function _render_hotspots_tab(
    performance::PerformanceState,
    workspace,
)::Nothing
    internal = workspace isa Workspace.Workspace ? workspace.profiler.report : nothing
    profile = internal isa Profiling.ProfileReport && internal.cpu !== nothing ?
        internal.cpu : performance.plot_sampling_profile

    if profile === nothing
        ig.Spacing()
        ig.TextDisabled("No CPU profile captured yet.")
        return nothing
    end

    ig.Spacing()
    if profile.truncated
        ig.TextColored(
            (0.95f0, 0.65f0, 0.15f0, 1.0f0),
            "Sample buffer filled before the rebuild finished — results are truncated.",
        )
    end
    total = max(profile.total_samples, 1)
    ig.TextDisabled(@sprintf(
        "%d samples at %.0f ms intervals   (%.2f s wall)",
        profile.total_samples, 1000 * profile.delay_seconds,
        profile.delay_seconds * profile.total_samples,
    ))
    ig.Spacing()

    # Table: self% / total% / self ms / function / file:line
    if _begin_perf_table("hotspots", 5, 0.0f0)
        ig.TableSetupScrollFreeze(0, 1)
        ig.TableSetupColumn("Self %",   ig.ImGuiTableColumnFlags_WidthFixed,   64.0f0)
        ig.TableSetupColumn("Total %",  ig.ImGuiTableColumnFlags_WidthFixed,   64.0f0)
        ig.TableSetupColumn("Self ms",  ig.ImGuiTableColumnFlags_WidthFixed,   68.0f0)
        ig.TableSetupColumn("Function", ig.ImGuiTableColumnFlags_WidthStretch, 0.0f0)
        ig.TableSetupColumn("File:line",ig.ImGuiTableColumnFlags_WidthStretch, 0.0f0)
        ig.TableHeadersRow()

        secs_per = profile.delay_seconds
        for row in profile.rows
            ig.TableNextRow()
            ig.TableNextColumn()
            _table_text(@sprintf("%.1f%%", 100 * row.self_samples / total))
            ig.TableNextColumn()
            _table_text(@sprintf("%.1f%%", 100 * row.samples / total))
            ig.TableNextColumn()
            self_ms = 1e3 * secs_per * row.self_samples
            _table_text(self_ms >= 1 ? @sprintf("%.0f", self_ms) : @sprintf("%.2f", self_ms))
            ig.TableNextColumn()
            _table_text(row.function_name)
            ig.TableNextColumn()
            loc = row.line > 0 ?
                "$(basename(row.file)):$(row.line)" : basename(row.file)
            _table_text(loc)
            if ig.BeginItemTooltip()
                ig.TextUnformatted(row.file)
                ig.EndTooltip()
            end
        end
        ig.EndTable()
    end
    return nothing
end

# ---------------------------------------------------------------------------
# Tab: Operations
# ---------------------------------------------------------------------------

"""Render the Operations tab: all perf-sample keys (frame timings, plot phases, etc.)."""
function _render_operations_tab(performance::PerformanceState)::Nothing
    timings = performance.timings
    allocs  = performance.allocations

    if isempty(timings)
        ig.TextDisabled("No operation timings yet.")
        return nothing
    end

    all_keys = sort!(collect(keys(timings)))
    if _begin_perf_table("operations", 6, 0.0f0)
        ig.TableSetupScrollFreeze(0, 1)
        ig.TableSetupColumn("Operation", ig.ImGuiTableColumnFlags_WidthStretch, 0.0f0)
        ig.TableSetupColumn("n",         ig.ImGuiTableColumnFlags_WidthFixed,  40.0f0)
        ig.TableSetupColumn("Mean ms",   ig.ImGuiTableColumnFlags_WidthFixed,  72.0f0)
        ig.TableSetupColumn("Max ms",    ig.ImGuiTableColumnFlags_WidthFixed,  72.0f0)
        ig.TableSetupColumn("Last ms",   ig.ImGuiTableColumnFlags_WidthFixed,  72.0f0)
        ig.TableSetupColumn("AllocLast", ig.ImGuiTableColumnFlags_WidthFixed,  80.0f0)
        ig.TableHeadersRow()

        for key in all_keys
            samples = timings[key]
            isempty(samples) && continue
            bytes = get(allocs, key, Int[])
            last_alloc_kb = isempty(bytes) ? 0 : bytes[end]

            ig.TableNextRow()
            ig.TableNextColumn(); _table_text(String(key))
            ig.TableNextColumn(); _table_text(string(length(samples)))
            ig.TableNextColumn(); _table_text(@sprintf("%.2f", mean(samples)))
            ig.TableNextColumn(); _table_text(@sprintf("%.2f", maximum(samples)))
            ig.TableNextColumn(); _table_text(@sprintf("%.2f", samples[end]))
            ig.TableNextColumn()
            if last_alloc_kb >= 1024^2
                _table_text(@sprintf("%.1f MB", last_alloc_kb / 1024^2))
            elseif last_alloc_kb >= 1024
                _table_text(@sprintf("%.1f KB", last_alloc_kb / 1024))
            else
                _table_text(@sprintf("%d B", last_alloc_kb))
            end
        end
        ig.EndTable()
    end

    ig.Spacing()
    if ig.Button("Clear timings")
        empty!(performance.timings)
        empty!(performance.allocations)
    end
    return nothing
end

# ---------------------------------------------------------------------------
# Tab: Memory & GL
# ---------------------------------------------------------------------------

"""Render the Memory & GL tab: RSS, GC stats, GL strings, and frame/tree counters."""
function _render_memory_gl_tab(state::BrowserState)::Nothing
    performance = state.performance

    # ---- Frame / tree counters ----
    raw_io = ig.GetIO()
    fps = unsafe_load(raw_io.Framerate)
    if fps > 0
        ig.Text(@sprintf("FPS: %.1f   Frame: %.2f ms", fps, 1000 / fps))
    end
    ig.Text("Tree nodes rendered: $(performance.node_count)")
    ig.Text("Items visible / rendered: $(performance.item_rows_visible) / $(performance.item_rows_rendered)")

    # ---- Process memory ----
    mem = _memory_snapshot()
    if mem.vmrss_kb !== nothing
        performance.memory_start_rss_kb === nothing &&
            (performance.memory_start_rss_kb = mem.vmrss_kb)
        start_rss = something(performance.memory_start_rss_kb)
        peak_rss  = max(something(performance.memory_peak_rss_kb, mem.vmrss_kb), mem.vmrss_kb)
        performance.memory_peak_rss_kb = peak_rss

        read_bytes = something(mem.read_bytes, 0)
        performance.memory_start_read_bytes === nothing &&
            (performance.memory_start_read_bytes = read_bytes)
        start_read = something(performance.memory_start_read_bytes)

        ig.Separator()
        ig.TextUnformatted("Process memory")
        if _begin_perf_table("memory", 2, 140.0f0)
            ig.TableSetupColumn("Field",  ig.ImGuiTableColumnFlags_WidthFixed,   110.0f0)
            ig.TableSetupColumn("Value",  ig.ImGuiTableColumnFlags_WidthStretch, 0.0f0)
            ig.TableHeadersRow()

            _mem_row(k, v) = begin
                ig.TableNextRow()
                ig.TableNextColumn(); _table_text(k)
                ig.TableNextColumn(); _table_text(v)
            end
            _mem_row("RSS",         _fmt_kb(mem.vmrss_kb))
            _mem_row("RSS Δstart",  _fmt_kb(mem.vmrss_kb - start_rss))
            _mem_row("RSS peak",    _fmt_kb(peak_rss))
            _mem_row("Anon",        _fmt_kb(mem.rssanon_kb))
            _mem_row("VSize",       _fmt_kb(mem.vmsize_kb))
            _mem_row("VPeak",       _fmt_kb(mem.vmpeak_kb))
            _mem_row("GC live",     _fmt_bytes(mem.gc_live_bytes))
            _mem_row("maxrss",      _fmt_bytes(mem.maxrss_bytes))
            _mem_row("IO Δstart",   _fmt_bytes(read_bytes - start_read))
            ig.EndTable()
        end
    else
        ig.Spacing()
        ig.TextDisabled("Process memory not available (Linux /proc only).")
    end

    # ---- Embedded Makie context ownership ----
    makie_stats = makie_context_stats()
    ig.Separator()
    ig.TextUnformatted("Embedded Makie")
    if _begin_perf_table("makie_contexts", 2, 160.0f0)
        ig.TableSetupColumn("Field", ig.ImGuiTableColumnFlags_WidthFixed, 110.0f0)
        ig.TableSetupColumn("Value", ig.ImGuiTableColumnFlags_WidthStretch, 0.0f0)
        ig.TableHeadersRow()

        _makie_row(k, v) = begin
            ig.TableNextRow()
            ig.TableNextColumn(); _table_text(k)
            ig.TableNextColumn(); _table_text(v)
        end
        _makie_row("Contexts", string(makie_stats.contexts))
        _makie_row("Created", string(makie_stats.created))
        _makie_row("Replaced", string(makie_stats.replaced))
        _makie_row("Destroyed", string(makie_stats.destroyed))
        _makie_row("Rendered", string(makie_stats.rendered))
        _makie_row("IDs", isempty(makie_stats.ids) ? "-" : join(sort(makie_stats.ids), ", "))
        ig.EndTable()
    end

    # ---- GL info ----
    gi = performance.gl_info
    if !isempty(gi)
        ig.Separator()
        ig.TextUnformatted("OpenGL")
        if _begin_perf_table("glinfo", 2, 0.0f0)
            ig.TableSetupColumn("Key",   ig.ImGuiTableColumnFlags_WidthFixed,   80.0f0)
            ig.TableSetupColumn("Value", ig.ImGuiTableColumnFlags_WidthStretch, 0.0f0)
            ig.TableHeadersRow()
            for k in (:vendor, :renderer, :version, :sl)
                haskey(gi, k) || continue
                ig.TableNextRow()
                ig.TableNextColumn(); _table_text(String(k))
                ig.TableNextColumn(); _table_text(gi[k])
            end
            ig.EndTable()
        end
    end
    return nothing
end

# ---------------------------------------------------------------------------
# Tab: Scan Profile
# ---------------------------------------------------------------------------

"""
Refresh the cached scan-profile rows when the project or TTL changes, then render them.

The per-kind and per-source-item tables reuse the existing `render_data_grid!` infrastructure
already wired up in Layout.jl; the grid states live in `PerformanceState`.
"""
function _render_scan_tab(state::BrowserState)::Nothing
    performance = state.performance
    workspace   = state.workspace

    if workspace isa Workspace.Workspace
        if performance.scan_profile_project !== workspace.project ||
           time() >= performance.scan_profile_refresh_at
            performance.scan_profile_project   = workspace.project
            performance.scan_kind_rows         = scan_profile_summary(workspace.project)
            performance.scan_source_rows       = scan_source_profile(workspace.project)
            performance.scan_profile_refresh_at = time() + 0.5
        end
    elseif performance.scan_profile_project !== nothing
        performance.scan_profile_project = nothing
        empty!(performance.scan_kind_rows)
        empty!(performance.scan_source_rows)
    end

    scan_rows   = performance.scan_kind_rows
    source_rows = performance.scan_source_rows

    if isempty(scan_rows)
        ig.TextDisabled("No scan timing yet — cache hit or no scan run.")
        return nothing
    end

    source_items = sum(row.source_items for row in scan_rows)
    items_total  = sum(row.items for row in scan_rows)
    ig.Text(@sprintf("Last scan: %d source items, %d items", source_items, items_total))
    ig.TextDisabled("Total = elapsed per source item; process/stats = summed item work.")

    ig.Spacing()
    ig.TextUnformatted("Per-kind summary")
    _render_scan_kind_table(scan_rows, performance.scan_kind_grid)

    ig.Spacing()
    ig.TextUnformatted("Slow source items")
    _render_scan_source_table(source_rows, performance.scan_source_grid)
    return nothing
end

# ---------------------------------------------------------------------------
# Tab: Live Plots
# ---------------------------------------------------------------------------

"""
Build the redraw-timings figure once: one axis, four phase lines.

Each phase carries its own `(x, y)` Observable pair so a longer series never collides with the
others on a shared x.
"""
function _make_timings_figure(lp::LivePlotsState)::Figure
    fig = Figure(size=(680, 220))
    axis = Axis(fig[1, 1]; xlabel="redraw #", ylabel="ms", title="Plot redraw phases")
    lines!(axis, lp.load_x, lp.load_obs; color=:steelblue, label="load")
    lines!(axis, lp.setup_x, lp.setup_obs; color=:darkorange, label="setup")
    lines!(axis, lp.data_x, lp.data_obs; color=:forestgreen, label="draw data")
    lines!(axis, lp.total_x, lp.total_obs; color=:crimson, label="total")
    axislegend(axis; position=:lt, labelsize=10)
    return fig
end

"""
Build the six-panel live build dashboard once.

All series share elapsed seconds on x. Separate axes keep percentages, rates, times, and memory from
distorting each other.
"""
function _make_build_figure(lp::LivePlotsState)::Figure
    fig = Figure(size=(680, 560))

    progress_axis = Axis(fig[1, 1]; ylabel="%", title="Progress")
    lines!(progress_axis, lp.elapsed_obs, lp.scan_pct_obs; color=:royalblue, label="scan")
    lines!(progress_axis, lp.elapsed_obs, lp.analysis_pct_obs;
        color=:darkorange, label="analysis")
    axislegend(progress_axis; position=:lt, labelsize=10)

    throughput_axis = Axis(fig[1, 2]; ylabel="items / s", title="Analysis throughput")
    lines!(throughput_axis, lp.elapsed_obs, lp.throughput_obs; color=:purple)

    writers_axis = Axis(fig[2, 1]; ylabel="writers", title="Writers busy concurrently")
    lines!(writers_axis, lp.elapsed_obs, lp.writers_busy_obs; color=:black)

    item_write_axis = Axis(fig[2, 2]; ylabel="ms / item", title="Processed write cost")
    lines!(item_write_axis, lp.elapsed_obs, lp.per_item_write_obs; color=:teal)

    cumulative_axis = Axis(
        fig[3, 1]; xlabel="elapsed (s)", ylabel="s", title="Cumulative write time")
    lines!(cumulative_axis, lp.elapsed_obs, lp.interp_cum_obs; color=:orchid, label="interp")
    lines!(cumulative_axis, lp.elapsed_obs, lp.processed_cum_obs;
        color=:seagreen, label="processed")
    lines!(cumulative_axis, lp.elapsed_obs, lp.stats_cum_obs;
        color=:goldenrod, label="stats")
    axislegend(cumulative_axis; position=:lt, labelsize=10)

    rss_axis = Axis(fig[3, 2]; xlabel="elapsed (s)", ylabel="GiB", title="Peak RSS")
    lines!(rss_axis, lp.elapsed_obs, lp.rss_obs; color=:darkorange)
    return fig
end

"""Push `v` into a bounded ring buffer of `capacity` elements (mutates `buf` in place)."""
@inline function _ring_push!(buf::Vector{Float32}, v::Float32, capacity::Int)::Nothing
    push!(buf, v)
    length(buf) > capacity && popfirst!(buf)
    return nothing
end

"""
Refresh redraw-timing Observables only when a new redraw landed.

Each phase keeps its own x values so error-only total samples cannot cause dimension mismatches.
"""
function _update_live_timings!(lp::LivePlotsState, timings::Dict{Symbol,Vector{Float64}})::Nothing
    draw = get(timings, :plot_draw, Float64[])
    length(draw) == lp.timings_seen && return nothing
    lp.timings_seen = length(draw)
    for (key, x_values, y_values) in (
        (:plot_load, lp.load_x, lp.load_obs),
        (:plot_setup, lp.setup_x, lp.setup_obs),
        (:plot_data, lp.data_x, lp.data_obs),
        (:plot_draw, lp.total_x, lp.total_obs),
    )
        values = get(timings, key, Float64[])
        sample_count = length(values)
        retained = min(lp.capacity, sample_count)
        if retained == 0
            isempty(y_values[]) || (x_values[] = Float32[]; y_values[] = Float32[])
            continue
        end
        first_index = sample_count - retained + 1
        x_values[] = Float32.(first_index:sample_count)
        y_values[] = Float32.(@view values[first_index:sample_count])
    end
    return nothing
end

"""
Sample live build counters into bounded dashboard buffers at most once per 250 ms.
"""
const _BUILD_SAMPLE_INTERVAL_NS = UInt64(250_000_000)  # 250 ms

function _sample_build_progress!(lp::LivePlotsState, workspace)::Nothing
    now = time_ns()
    lp.t0_ns == 0 && (lp.t0_ns = now)
    now - lp.last_sample_ns < _BUILD_SAMPLE_INTERVAL_NS && return nothing
    lp.last_sample_ns = now
    elapsed = Float64(now - lp.t0_ns) / 1e9

    scan_pct = 0.0f0
    analysis_pct = 0.0f0
    completed = 0
    processed_write_ns = Int64(0)
    writer_busy_ns = Int64(0)
    interpreted_seconds = 0.0f0
    processed_seconds = 0.0f0
    stats_seconds = 0.0f0

    if workspace isa Workspace.Workspace
        progress = workspace.scan.progress
        total_si = progress.total_source_items
        scan_pct = total_si > 0 ?
            Float32(100 * progress.processed_source_items / total_si) : 0.0f0

        completed, total_processing = lock(workspace.processing.lock) do
            (workspace.processing.completed, workspace.processing.total)
        end
        analysis_pct = total_processing > 0 ?
            Float32(100 * completed / total_processing) : 0.0f0

        metrics = workspace.metrics
        processed_write_ns = metrics.processed_write_ns[]
        interpreted_seconds = Float32(metrics.interpreted_write_ns[] / 1e9)
        processed_seconds = Float32(metrics.processed_write_ns[] / 1e9)
        stats_seconds = Float32(metrics.stats_write_ns[] / 1e9)
        writer_busy_ns = workspace.cache.db.writer_busy_ns[]
    end

    elapsed_delta = elapsed - lp.last_elapsed_s
    completed_delta = completed - lp.last_completed
    busy_delta = writer_busy_ns - lp.last_writer_busy_ns
    write_delta = processed_write_ns - lp.last_processed_write_ns
    throughput = elapsed_delta > 0 && completed_delta >= 0 ?
        Float32(completed_delta / elapsed_delta) : 0.0f0
    writers_busy = elapsed_delta > 0 && busy_delta >= 0 ?
        Float32(busy_delta / 1e9 / elapsed_delta) : 0.0f0
    per_item_write = completed_delta > 0 && write_delta >= 0 ?
        Float32(write_delta / 1e6 / completed_delta) : 0.0f0
    peak_rss_gib = Float32(Sys.maxrss() / 1024^3)

    lp.last_elapsed_s = elapsed
    lp.last_completed = completed
    lp.last_processed_write_ns = processed_write_ns
    lp.last_writer_busy_ns = writer_busy_ns

    capacity = lp.capacity
    _ring_push!(lp.elapsed_buf, Float32(elapsed), capacity)
    _ring_push!(lp.scan_pct_buf, scan_pct, capacity)
    _ring_push!(lp.analysis_pct_buf, analysis_pct, capacity)
    _ring_push!(lp.throughput_buf, throughput, capacity)
    _ring_push!(lp.writers_busy_buf, writers_busy, capacity)
    _ring_push!(lp.per_item_write_buf, per_item_write, capacity)
    _ring_push!(lp.interp_cum_buf, interpreted_seconds, capacity)
    _ring_push!(lp.processed_cum_buf, processed_seconds, capacity)
    _ring_push!(lp.stats_cum_buf, stats_seconds, capacity)
    _ring_push!(lp.rss_buf, peak_rss_gib, capacity)

    lp.elapsed_obs[] = copy(lp.elapsed_buf)
    lp.scan_pct_obs[] = copy(lp.scan_pct_buf)
    lp.analysis_pct_obs[] = copy(lp.analysis_pct_buf)
    lp.throughput_obs[] = copy(lp.throughput_buf)
    lp.writers_busy_obs[] = copy(lp.writers_busy_buf)
    lp.per_item_write_obs[] = copy(lp.per_item_write_buf)
    lp.interp_cum_obs[] = copy(lp.interp_cum_buf)
    lp.processed_cum_obs[] = copy(lp.processed_cum_buf)
    lp.stats_cum_obs[] = copy(lp.stats_cum_buf)
    lp.rss_obs[] = copy(lp.rss_buf)
    return nothing
end

"""
Save a figure to PNG (or the user's chosen extension), then optionally dump the raw buffer series
to a companion CSV (same stem + `_data.csv`).
"""
function _export_live_figure(
    figure::Figure,
    default_name::String,
    series::Vector{Pair{String,Vector{Float32}}},
)::String
    path = save_file(default_name; filterlist="png,jpg,jpeg,svg,pdf")
    isempty(path) && return ""
    isempty(splitext(path)[2]) && (path *= ".png")
    try
        Makie.save(path, figure)
    catch err
        bt = catch_backtrace()
        summary = first(split(sprint(showerror, err), '\n'; limit=2))
        @error("Live-plot export failed\nOutput: $path", exception=(err, bt))
        return "Export failed: $summary. See the console for full details."
    end
    if !isempty(series)
        csv_path = splitext(path)[1] * "_data.csv"
        try
            open(csv_path, "w") do io
                println(io, join(first.(series), ","))
                rows = minimum(length(last(entry)) for entry in series)
                for i in 1:rows
                    println(io, join((string(last(s)[i]) for s in series), ","))
                end
            end
        catch err
            @warn "Could not write companion CSV" path=csv_path exception=err
        end
    end
    return ""
end

"""
Render redraw timings and aggregate build plots.

`MakieFigure` polls GLMakie for pending updates and renders only dirty frames.
"""
function _render_live_plots_tab!(state::BrowserState)::Nothing
    lp        = state.performance.live_plots
    workspace = state.workspace

    lp.timings_figure === nothing && (lp.timings_figure = _make_timings_figure(lp))
    _update_live_timings!(lp, state.performance.timings)
    ig.TextUnformatted("Plot redraw phase timings (ms per redraw)")
    ig.SameLine()
    if ig.Button("Export##live_timings_export")
        series = Pair{String,Vector{Float32}}[
            "redraw" => copy(lp.total_x[]),
            "load_ms" => copy(lp.load_obs[]),
            "setup_ms" => copy(lp.setup_obs[]),
            "data_ms" => copy(lp.data_obs[]),
            "total_ms" => copy(lp.total_obs[]),
        ]
        lp.timings_export_error = _export_live_figure(
            lp.timings_figure, "perf-timings.png", series)
    end
    isempty(lp.timings_export_error) || ig.TextWrapped(lp.timings_export_error)
    MakieFigure(
        "perf_live_timings", lp.timings_figure;
        auto_resize_x=true,
        auto_resize_y=false,
    )
    ig.Spacing()
    ig.Separator()
    ig.Spacing()

    lp.build_figure === nothing && (lp.build_figure = _make_build_figure(lp))
    _sample_build_progress!(lp, workspace)
    ig.TextUnformatted("Live build dashboard")
    ig.SameLine()
    if ig.Button("Export##live_build_export")
        series = Pair{String,Vector{Float32}}[
            "elapsed_s" => copy(lp.elapsed_obs[]),
            "scan_pct" => copy(lp.scan_pct_obs[]),
            "analysis_pct" => copy(lp.analysis_pct_obs[]),
            "throughput_per_s" => copy(lp.throughput_obs[]),
            "writers_busy" => copy(lp.writers_busy_obs[]),
            "per_item_write_ms" => copy(lp.per_item_write_obs[]),
            "interp_cum_s" => copy(lp.interp_cum_obs[]),
            "processed_cum_s" => copy(lp.processed_cum_obs[]),
            "stats_cum_s" => copy(lp.stats_cum_obs[]),
            "peak_rss_gib" => copy(lp.rss_obs[]),
        ]
        lp.build_export_error = _export_live_figure(
            lp.build_figure, "perf-build.png", series)
    end
    isempty(lp.build_export_error) || ig.TextWrapped(lp.build_export_error)

    MakieFigure("perf_live_build", lp.build_figure; auto_resize_x=true, auto_resize_y=false)

    return nothing
end

"""Run one internal profiler UI action and retain an actionable error."""
function _profile_action!(action::Function, profiler::Profiling.ProfileSession)::Nothing
    try
        action()
        profiler.error = ""
    catch error
        profiler.error = sprint(showerror, error)
        profiler.state = :error
        @error "Internal profiler action failed" exception=(error, catch_backtrace())
    end
    return nothing
end

"""Choose one exact structured-event dimension or show all values."""
function _profile_filter_combo!(
    label::String,
    current::Symbol,
    values::Vector{Symbol},
)::Symbol
    preview = current === :all ? "All" : String(current)
    ig.SetNextItemWidth(170.0f0)
    if ig.BeginCombo(label, preview)
        ig.Selectable("All", current === :all) && (current = :all)
        for value in values
            ig.Selectable(String(value), current === value) && (current = value)
        end
        ig.EndCombo()
    end
    return current
end

"""Render opt-in internal capture controls, summaries, events, and process counters."""
function _render_internal_profile_tab!(
    state::BrowserState,
    workspace::Workspace.Workspace,
)::Nothing
    profiler = workspace.profiler
    ig.Text("State: $(profiler.state)  Events: $(Profiling.event_count(profiler))  " *
            "Dropped: $(profiler.dropped_events[])")

    can_start = profiler.state in (:idle, :complete)
    can_start || ig.BeginDisabled()
    if ig.Button("Start capture")
        _profile_action!(profiler) do
            Workspace.start_internal_profile!(workspace)
        end
    end
    can_start || ig.EndDisabled()
    ig.SameLine()

    can_stop = profiler.state in (:preparing, :recording)
    can_stop || ig.BeginDisabled()
    if ig.Button("Stop capture")
        _profile_action!(profiler) do
            Workspace.stop_internal_profile!(workspace)
        end
    end
    can_stop || ig.EndDisabled()
    ig.SameLine()

    can_reset = profiler.state in (:complete, :error)
    can_reset || ig.BeginDisabled()
    if ig.Button("Reset")
        _profile_action!(profiler) do
            Workspace.reset_internal_profile!(workspace)
        end
    end
    can_reset || ig.EndDisabled()
    ig.SameLine()

    can_export = profiler.report isa Profiling.ProfileReport
    can_export || ig.BeginDisabled()
    if ig.Button("Export Perfetto")
        _profile_action!(profiler) do
            path = save_file("measurementbrowser-profile.json"; filterlist="json")
            isempty(path) || Workspace.export_internal_profile!(workspace, path)
        end
    end
    can_export || ig.EndDisabled()

    isempty(profiler.error) || ig.TextWrapped("Profiler error: $(profiler.error)")
    profiler.state === :preparing && ig.TextDisabled(
        "Canceling current work before a clean profiled rebuild.",
    )

    report = profiler.report
    if report isa Profiling.ProfileReport
        ig.Spacing()
        ig.Text(@sprintf(
            "Duration %.2f s  Counters %d  CPU samples %d",
            (report.stopped_ns - report.started_ns) / 1e9,
            length(report.counters),
            report.cpu === nothing ? 0 : report.cpu.total_samples,
        ))
        if _begin_perf_table("internal_summary", 10, 260.0f0)
            for (name, width) in (
                ("Operation", 160.0f0), ("n", 45.0f0), ("Total", 65.0f0),
                ("p50", 58.0f0), ("p90", 58.0f0), ("p99", 58.0f0),
                ("Max", 58.0f0), ("Wait", 65.0f0), ("Service", 65.0f0),
                ("Batch p50/90/max", 120.0f0),
            )
                ig.TableSetupColumn(name, ig.ImGuiTableColumnFlags_WidthFixed, width)
            end
            ig.TableHeadersRow()
            for row in report.summary
                ig.TableNextRow()
                for value in (
                    "$(row.category)/$(row.operation)", string(row.count),
                    @sprintf("%.1f", row.total_ms), @sprintf("%.2f", row.median_ms),
                    @sprintf("%.2f", row.p90_ms), @sprintf("%.2f", row.p99_ms),
                    @sprintf("%.2f", row.max_ms), @sprintf("%.1f", row.wait_ms),
                    @sprintf("%.1f", row.service_ms),
                    @sprintf("%.0f / %.0f / %.0f",
                        row.median_batch, row.p90_batch, row.max_batch),
                )
                    ig.TableNextColumn(); _table_text(value)
                end
            end
            ig.EndTable()
        end
    end

    counter = Profiling.latest_counter(profiler)
    if counter !== nothing
        ig.Text(@sprintf(
            "RSS %.1f MiB  GC %.1f ms / %d pauses (max %.1f ms)  Queue %d  Writer %.1f ms wait / %.1f ms service",
            counter.rss_bytes / 1024^2,
            counter.gc_total_ns / 1e6,
            counter.gc_pause_count,
            counter.gc_max_pause_ns / 1e6,
            counter.queue_depth,
            counter.writer_wait_ns / 1e6,
            counter.writer_busy_ns / 1e6,
        ))
    end

    available = report isa Profiling.ProfileReport ?
        report.events : Profiling.recent_events(profiler; limit=2_000)
    categories = sort!(unique([event.category for event in available]))
    operations = sort!(unique([event.operation for event in available]))
    state.performance.profile_category_filter = _profile_filter_combo!(
        "Category##profile_category",
        state.performance.profile_category_filter,
        categories,
    )
    ig.SameLine()
    state.performance.profile_operation_filter = _profile_filter_combo!(
        "Operation##profile_operation",
        state.performance.profile_operation_filter,
        operations,
    )
    category_filter = state.performance.profile_category_filter
    operation_filter = state.performance.profile_operation_filter
    shown = 0
    for index in length(available):-1:1
        event = available[index]
        ((category_filter === :all || event.category === category_filter) &&
         (operation_filter === :all || event.operation === operation_filter)) &&
            (shown += 1)
        shown >= 200 && break
    end
    if shown > 0 && _begin_perf_table("internal_events", 6, 220.0f0)
        for name in ("Category", "Operation", "ms", "Thread", "Status", "Batch")
            ig.TableSetupColumn(name)
        end
        ig.TableHeadersRow()
        rendered = 0
        for index in length(available):-1:1
            event = available[index]
            (category_filter === :all || event.category === category_filter) &&
                (operation_filter === :all || event.operation === operation_filter) ||
                continue
            ig.TableNextRow()
            for value in (
                String(event.category), String(event.operation),
                @sprintf("%.3f", event.duration_ns / 1e6),
                string(event.start_thread), String(event.status),
                string(event.attributes.batch_size),
            )
                ig.TableNextColumn(); _table_text(value)
            end
            rendered += 1
            rendered >= 200 && break
        end
        ig.EndTable()
    end

    return nothing
end

# ---------------------------------------------------------------------------
# Top-level entry point
# ---------------------------------------------------------------------------

"""
Render the Performance debug window with a tab bar.

Tabs:
- Live         — GLMakie line plots for rolling plot-redraw timings
- Plot Redraw  — load / setup / draw / total phase timings (table)
- Hotspots     — sampling profiler results (also the "Profile full rebuild" button)
- Operations   — all keyed `_append_perf_sample!` timings
- Scan         — per-kind and per-source-item scan timings
- Memory & GL  — process RSS, GC, I/O, OpenGL strings, frame counters
- Internal     — opt-in structured engine spans, counters, CPU samples, and export
"""
function render_perf_window(state::BrowserState)::Nothing
    state.show_performance_window || return nothing

    if ig.Begin("Performance###perf_window")
        workspace = state.workspace

        if ig.BeginTabBar("##perf_tabs")
            if ig.BeginTabItem("Live")
                ig.Spacing()
                _render_live_plots_tab!(state)
                ig.EndTabItem()
            end

            if ig.BeginTabItem("Plot Redraw")
                ig.Spacing()
                _render_plot_redraw_tab(state.performance)
                ig.EndTabItem()
            end

            if ig.BeginTabItem("Hotspots")
                ig.Spacing()
                _render_hotspots_tab(state.performance, workspace)
                ig.EndTabItem()
            end

            if ig.BeginTabItem("Operations")
                ig.Spacing()
                _render_operations_tab(state.performance)
                ig.EndTabItem()
            end

            if ig.BeginTabItem("Scan")
                ig.Spacing()
                _render_scan_tab(state)
                ig.EndTabItem()
            end

            if ig.BeginTabItem("Memory & GL")
                ig.Spacing()
                _render_memory_gl_tab(state)
                ig.EndTabItem()
            end

            if workspace isa Workspace.Workspace && workspace.profiler.enabled &&
               ig.BeginTabItem("Internal")
                ig.Spacing()
                _render_internal_profile_tab!(state, workspace)
                ig.EndTabItem()
            end

            ig.EndTabBar()
        end
    end
    ig.End()
    return nothing
end
