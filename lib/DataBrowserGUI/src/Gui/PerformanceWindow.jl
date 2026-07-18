using Printf
using Statistics: mean
import CImGui as ig
using NativeFileDialog: save_file

import DataBrowserProfiling as Profiling
using DataBrowserAPI: scan_profile_summary, scan_source_profile
import DataBrowserCore.Workspace

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
# Tab: Project — where project code (recipes) spends its time
# ---------------------------------------------------------------------------

"""Render mean/max/last/alloc for the plot-recipe phases (load, setup, draw, total)."""
function _render_plot_phase_table(performance::PerformanceState)::Nothing
    timings = performance.timings
    allocs  = performance.allocations

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

"""Copy the trailing `capacity` samples from `values` into `buf`."""
function _retain_tail!(buf::Vector{Float32}, values::Vector{Float64}, capacity::Int)::Nothing
    empty!(buf)
    retained = min(capacity, length(values))
    retained == 0 && return nothing
    append!(buf, Float32.(@view values[end - retained + 1:end]))
    return nothing
end

"""Render one bounded history as an ImGui sparkline with the latest value as overlay text."""
function _render_sparkline!(
    id::String,
    label::String,
    values::Vector{Float32},
    unit::String;
    height::Float32=60.0f0,
)::Nothing
    ig.TextUnformatted(label)
    if isempty(values)
        ig.TextDisabled(" (no samples)")
        return nothing
    end
    ig.SameLine()
    ig.TextDisabled(@sprintf("  %.2f %s", values[end], unit))
    width = max(120.0f0, ig.GetContentRegionAvail().x)
    # PlotLines takes positional arguments only: values_offset, overlay_text,
    # scale_min, scale_max, graph_size.
    ig.PlotLines(
        "##perf_spark_$id",
        values,
        length(values),
        0,
        C_NULL,
        0.0f0,
        ig.FLT_MAX,
        (width, height),
    )
    return nothing
end

"""Write aligned sparkline buffers to a CSV file chosen by the user."""
function _export_sparklines_csv(
    default_name::String,
    series::Vector{Pair{String,Vector{Float32}}},
)::String
    path = save_file(default_name; filterlist="csv")
    isempty(path) && return ""
    isempty(splitext(path)[2]) && (path *= ".csv")
    try
        open(path, "w") do io
            println(io, join(first.(series), ","))
            rows = minimum(length(last(entry)) for entry in series)
            for i in 1:rows
                println(io, join((string(last(s)[i]) for s in series), ","))
            end
        end
    catch err
        bt = catch_backtrace()
        summary = first(split(sprint(showerror, err), '\n'; limit=2))
        @error("Sparkline export failed\nOutput: $path", exception=(err, bt))
        return "Export failed: $summary. See the console for full details."
    end
    return ""
end

"""
Refresh redraw-timing ring buffers only when a new redraw landed.
"""
function _update_live_timings!(lp::LivePlotsState, timings::Dict{Symbol,Vector{Float64}})::Nothing
    draw = get(timings, :plot_draw, Float64[])
    length(draw) == lp.timings_seen && return nothing
    lp.timings_seen = length(draw)
    for (key, buf) in (
        (:plot_load, :load_buf),
        (:plot_setup, :setup_buf),
        (:plot_data, :data_buf),
        (:plot_draw, :total_buf),
    )
        _retain_tail!(getfield(lp, buf), get(timings, key, Float64[]), lp.capacity)
    end
    return nothing
end

"""
Render project-code timings: per-kind scan recipes, slow source items, and plot-recipe phases.

The per-kind and per-source-item tables reuse the existing `render_data_grid!` infrastructure
already wired up in Layout.jl; the grid states live in `PerformanceState`.
"""
function _render_project_tab!(state::BrowserState)::Nothing
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
    else
        source_items = sum(row.source_items for row in scan_rows)
        items_total  = sum(row.items for row in scan_rows)
        ig.Text(@sprintf("Last scan: %d source items, %d items", source_items, items_total))
        ig.TextDisabled("Total = elapsed per source item; process/analyze = summed item work.")

        ig.Spacing()
        ig.TextUnformatted("Per-kind summary")
        _render_scan_kind_table(scan_rows, performance.scan_kind_grid)

        ig.Spacing()
        ig.TextUnformatted("Slow source items")
        _render_scan_source_table(source_rows, performance.scan_source_grid)
    end

    ig.Spacing()
    ig.Separator()
    ig.TextUnformatted("Plot recipes (setup / draw callbacks)")
    _render_plot_phase_table(performance)

    lp = performance.live_plots
    _update_live_timings!(lp, performance.timings)
    if ig.Button("Export##live_timings_export")
        series = Pair{String,Vector{Float32}}[
            "load_ms" => copy(lp.load_buf),
            "setup_ms" => copy(lp.setup_buf),
            "data_ms" => copy(lp.data_buf),
            "total_ms" => copy(lp.total_buf),
        ]
        lp.timings_export_error = _export_sparklines_csv("perf-timings.csv", series)
    end
    isempty(lp.timings_export_error) || ig.TextWrapped(lp.timings_export_error)
    _render_sparkline!("load", "Load", lp.load_buf, "ms")
    _render_sparkline!("setup", "Setup", lp.setup_buf, "ms")
    _render_sparkline!("data", "Draw data", lp.data_buf, "ms")
    _render_sparkline!("total", "Total", lp.total_buf, "ms")
    return nothing
end

# ---------------------------------------------------------------------------
# Tab: Memory — cross-platform process and workspace memory
# ---------------------------------------------------------------------------

"""Render process RSS/GC totals and, when a workspace is open, its index footprint counts."""
function _render_memory_tab(state::BrowserState)::Nothing
    performance = state.performance

    rss = Profiling.process_rss_bytes()
    performance.memory_start_rss === nothing && (performance.memory_start_rss = rss)
    start_rss = something(performance.memory_start_rss)
    performance.memory_peak_rss = max(performance.memory_peak_rss, rss)
    gc = Base.gc_num()
    gc_live = Int64(Base.gc_live_bytes())

    ig.TextUnformatted("Process memory")
    if _begin_perf_table("memory", 2, 240.0f0)
        ig.TableSetupColumn("Field",  ig.ImGuiTableColumnFlags_WidthFixed,   130.0f0)
        ig.TableSetupColumn("Value",  ig.ImGuiTableColumnFlags_WidthStretch, 0.0f0)
        ig.TableHeadersRow()

        _mem_row(k, v) = begin
            ig.TableNextRow()
            ig.TableNextColumn(); _table_text(k)
            ig.TableNextColumn(); _table_text(v)
        end
        _mem_row("RSS",           _fmt_bytes(rss))
        _mem_row("RSS Δstart",    _fmt_bytes(rss - start_rss))
        _mem_row("RSS peak",      _fmt_bytes(performance.memory_peak_rss))
        _mem_row("RSS peak (OS)", _fmt_bytes(Int64(Sys.maxrss())))
        _mem_row("GC live",       _fmt_bytes(gc_live))
        _mem_row("Heap slack",    _fmt_bytes(max(Int64(0), rss - gc_live)))
        _mem_row("GC time",       @sprintf("%.1f s", gc.total_time / 1e9))
        _mem_row("GC pauses",     @sprintf(
            "%d (%d full, max %.1f ms)", gc.pause, gc.full_sweep, gc.max_pause / 1e6))
        _mem_row("Allocated",     _fmt_bytes(Int64(gc.total_allocd)))
        ig.EndTable()
    end

    workspace = state.workspace
    if workspace isa Workspace.Workspace
        index = workspace.index
        ig.Separator()
        ig.Text("Index: $(length(index.items)) items in " *
                "$(length(index.collections.records)) collections")
        ig.Text("Metadata held: $(length(index.item_metadata))   " *
                "Analysis errors: $(length(index.analysis_errors))")
    end
    return nothing
end

# ---------------------------------------------------------------------------
# Internal tab: Build — live engine pipeline
# ---------------------------------------------------------------------------

"""Push `v` into a bounded ring buffer of `capacity` elements (mutates `buf` in place)."""
@inline function _ring_push!(buf::Vector{Float32}, v::Float32, capacity::Int)::Nothing
    push!(buf, v)
    length(buf) > capacity && popfirst!(buf)
    return nothing
end

const _BUILD_SAMPLE_INTERVAL_NS = UInt64(250_000_000)  # 250 ms

"""
Sample live engine counters into the build sparkline buffers at most once per 250 ms.

Called whenever the Performance window is open (any tab), so the series covers the whole build
rather than the moments the Build tab happened to be visible.
"""
function _sample_build!(lp::LivePlotsState, workspace::Workspace.Workspace)::Nothing
    now = time_ns()
    lp.t0_ns == 0 && (lp.t0_ns = now)
    now - lp.last_sample_ns < _BUILD_SAMPLE_INTERVAL_NS && return nothing
    lp.last_sample_ns = now
    elapsed = Float64(now - lp.t0_ns) / 1e9

    completed, _total, active = Workspace.work_counts(workspace)
    staged = Workspace.cache_pending_counts(workspace.cache.db)

    elapsed_delta = elapsed - lp.last_elapsed_s
    completed_delta = completed - lp.last_completed
    throughput = elapsed_delta > 0 && completed_delta >= 0 ?
        Float32(completed_delta / elapsed_delta) : 0.0f0
    lp.last_elapsed_s = elapsed
    lp.last_completed = completed

    capacity = lp.capacity
    _ring_push!(lp.elapsed_buf, Float32(elapsed), capacity)
    _ring_push!(lp.throughput_buf, throughput, capacity)
    _ring_push!(lp.active_buf, Float32(active), capacity)
    _ring_push!(lp.pending_buf, Float32(staged.rows), capacity)
    _ring_push!(lp.rss_buf, Float32(Profiling.process_rss_bytes() / 1024^3), capacity)
    return nothing
end

"""Render live pipeline counters and the build sparklines."""
function _render_build_tab!(state::BrowserState, workspace::Workspace.Workspace)::Nothing
    lp = state.performance.live_plots

    completed, total, active = Workspace.work_counts(workspace)
    staged = Workspace.cache_pending_counts(workspace.cache.db)
    ig.Text(@sprintf("Analysis: %d / %d done   %d active", completed, total, active))
    ig.Text("Scan: $(workspace.scan.state)   Cache: $(workspace.cache_state)")
    ig.Text(@sprintf(
        "Cache write buffer: %d items / %d rows pending", staged.items, staged.rows))
    ig.Spacing()

    if ig.Button("Export##live_build_export")
        series = Pair{String,Vector{Float32}}[
            "elapsed_s" => copy(lp.elapsed_buf),
            "throughput_per_s" => copy(lp.throughput_buf),
            "active_work" => copy(lp.active_buf),
            "pending_rows" => copy(lp.pending_buf),
            "rss_gib" => copy(lp.rss_buf),
        ]
        lp.build_export_error = _export_sparklines_csv("perf-build.csv", series)
    end
    isempty(lp.build_export_error) || ig.TextWrapped(lp.build_export_error)
    _render_sparkline!("throughput", "Throughput", lp.throughput_buf, "items/s")
    _render_sparkline!("active", "Active work", lp.active_buf, "nodes")
    _render_sparkline!("pending", "Pending cache rows", lp.pending_buf, "rows")
    _render_sparkline!("rss", "RSS", lp.rss_buf, "GiB")
    return nothing
end

"""Render CPU hotspots from an explicit plot sampling session."""
function _render_hotspots_section(
    performance::PerformanceState,
    workspace,
)::Nothing
    profile = performance.plot_sampling_profile

    if profile === nothing
        ig.TextDisabled("No CPU profile captured yet.")
        return nothing
    end

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

    # Table: self% / total% / self ms / function / file:line
    if _begin_perf_table("hotspots", 5, 220.0f0)
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
# Internal tab: Frames — app-overhead debugging
# ---------------------------------------------------------------------------

"""Render frame rate, per-panel timings, and OpenGL strings (app overhead, not project code)."""
function _render_frames_tab(state::BrowserState)::Nothing
    performance = state.performance

    raw_io = ig.GetIO()
    fps = unsafe_load(raw_io.Framerate)
    if fps > 0
        ig.Text(@sprintf("FPS: %.1f   Frame: %.2f ms", fps, 1000 / fps))
    end
    ig.Text("Tree nodes rendered: $(performance.node_count)")
    ig.Text("Items visible / rendered: $(performance.item_rows_visible) / $(performance.item_rows_rendered)")
    ig.Spacing()

    timings = performance.timings
    allocs  = performance.allocations
    if isempty(timings)
        ig.TextDisabled("No operation timings yet.")
    else
        ig.TextDisabled(
            "Panel timings cover the Julia callback only; ImGui::Render() and OpenGL " *
            "draw run afterward and are not included.",
        )
        ig.Spacing()
        if _begin_perf_table("operations", 6, 260.0f0)
        ig.TableSetupScrollFreeze(0, 1)
        ig.TableSetupColumn("Operation", ig.ImGuiTableColumnFlags_WidthStretch, 0.0f0)
        ig.TableSetupColumn("n",         ig.ImGuiTableColumnFlags_WidthFixed,  40.0f0)
        ig.TableSetupColumn("Mean ms",   ig.ImGuiTableColumnFlags_WidthFixed,  72.0f0)
        ig.TableSetupColumn("Max ms",    ig.ImGuiTableColumnFlags_WidthFixed,  72.0f0)
        ig.TableSetupColumn("Last ms",   ig.ImGuiTableColumnFlags_WidthFixed,  72.0f0)
        ig.TableSetupColumn("AllocLast", ig.ImGuiTableColumnFlags_WidthFixed,  80.0f0)
        ig.TableHeadersRow()

        for key in sort!(collect(keys(timings)))
            samples = timings[key]
            isempty(samples) && continue
            bytes = get(allocs, key, Int[])
            last_alloc = isempty(bytes) ? 0 : bytes[end]

            ig.TableNextRow()
            ig.TableNextColumn(); _table_text(String(key))
            ig.TableNextColumn(); _table_text(string(length(samples)))
            ig.TableNextColumn(); _table_text(@sprintf("%.2f", mean(samples)))
            ig.TableNextColumn(); _table_text(@sprintf("%.2f", maximum(samples)))
            ig.TableNextColumn(); _table_text(@sprintf("%.2f", samples[end]))
            ig.TableNextColumn(); _table_text(_fmt_bytes(last_alloc))
        end
        ig.EndTable()
        end
    end
    if ig.Button("Clear timings")
        empty!(performance.timings)
        empty!(performance.allocations)
    end

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
# Top-level entry point
# ---------------------------------------------------------------------------

"""
Render the Performance window with normal project and workspace diagnostics:

Project diagnostics focus on where project code spends time:
- Project — per-kind recipe timings, slow source items, plot-recipe phases
- Memory  — cross-platform process RSS/GC and workspace index footprint

Workspace diagnostics show live pipeline counters, callback timings, and rendering overhead.
"""
function render_perf_window(state::BrowserState)::Nothing
    state.show_performance_window || return nothing

    if ig.Begin("Performance###perf_window")
        workspace = state.workspace
        workspace isa Workspace.Workspace && _sample_build!(state.performance.live_plots, workspace)

        if ig.BeginTabBar("##perf_tabs")
            if ig.BeginTabItem("Project")
                ig.Spacing()
                _render_project_tab!(state)
                ig.EndTabItem()
            end

            if ig.BeginTabItem("Memory")
                ig.Spacing()
                _render_memory_tab(state)
                ig.EndTabItem()
            end

            if workspace isa Workspace.Workspace && ig.BeginTabItem("Workspace")
                ig.Spacing()
                _render_build_tab!(state, workspace)
                ig.EndTabItem()
            end

            if ig.BeginTabItem("Callbacks")
                ig.Spacing()
                _render_frames_tab(state)
                ig.Separator()
                ig.TextUnformatted("Plot sampling")
                _render_hotspots_section(state.performance, workspace)
                ig.EndTabItem()
            end

            ig.EndTabBar()
        end
    end
    ig.End()
    return nothing
end
