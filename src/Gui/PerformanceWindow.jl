using Printf
using Statistics: mean
import CImGui as ig
using GLMakie: Axis, Figure, Observable, axislegend, lines!
import GLMakie.Makie as Makie
using NativeFileDialog: save_file

import ..Profiling
using ..Projects: scan_profile_summary, scan_source_profile
import ..Workspace
using .MakieImguiIntegration: MakieFigure

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

const _LIVE_PHASE_COLORS = (
    load  = :steelblue,
    setup = :darkorange,
    data  = :forestgreen,
    total = :crimson,
)

"""
Build the timings figure once, using Observables that are updated in place each frame.

Returns a `Figure` with one `Axis` showing the four plot-redraw phases as separate lines.
The Observables in `lp` are mutated directly by [`_update_live_timings!`](@ref).
"""
function _make_timings_figure(lp::LivePlotsState)::Figure
    fig = Figure(size=(600, 200))
    ax  = Axis(fig[1, 1]; xlabel="sample", ylabel="ms", title="Plot redraw phases")
    lines!(ax, lp.timings_x_obs, lp.load_obs;  color=_LIVE_PHASE_COLORS.load,  label="load")
    lines!(ax, lp.timings_x_obs, lp.setup_obs; color=_LIVE_PHASE_COLORS.setup, label="setup")
    lines!(ax, lp.timings_x_obs, lp.data_obs;  color=_LIVE_PHASE_COLORS.data,  label="draw data")
    lines!(ax, lp.timings_x_obs, lp.total_obs; color=_LIVE_PHASE_COLORS.total, label="total")
    axislegend(ax; position=:lt, labelsize=11)
    return fig
end

"""
Build the build-progress figure once.

Lines show scan %, analysis %, interp-write cumulative ms, processed-write cumulative ms, and RSS.
Because the scales differ wildly, scan/analysis are on the left axis (0–100 %) and
writer/RSS are on the right axis (ms or MiB).
"""
function _make_build_figure(lp::LivePlotsState)::Figure
    fig  = Figure(size=(600, 200))
    ax_l = Axis(fig[1, 1]; xlabel="sample", ylabel="%  /  ms", title="Build progress")
    ax_r = Axis(fig[1, 1]; ylabel="Peak RSS MiB",
                yaxisposition=:right, yticklabelcolor=:grey50,
                ytickcolor=:grey50, rightspinecolor=:grey50)
    # Link x axes so pan/zoom stay in sync
    lines!(ax_l, lp.build_x_obs, lp.scan_obs;             color=:royalblue,     label="scan %")
    lines!(ax_l, lp.build_x_obs, lp.analysis_obs;         color=:darkorange,    label="analysis %")
    lines!(ax_l, lp.build_x_obs, lp.interp_write_obs;     color=:orchid,        label="interp write ms")
    lines!(ax_l, lp.build_x_obs, lp.processed_write_obs;  color=:seagreen,      label="proc write ms")
    lines!(ax_r, lp.build_x_obs, lp.rss_obs;              color=:grey50,        label="RSS MiB")
    axislegend(ax_l; position=:lt, labelsize=11)
    return fig
end

"""Push `v` into a bounded ring buffer of `capacity` elements (mutates `buf` in place)."""
@inline function _ring_push!(buf::Vector{Float32}, v::Float32, capacity::Int)::Nothing
    push!(buf, v)
    length(buf) > capacity && popfirst!(buf)
    return nothing
end

"""
Copy the latest content of `state.performance.timings` into the timings figure Observables.

Only touches the four canonical plot-redraw phase keys; other keys are ignored.
"""
function _update_live_timings!(lp::LivePlotsState, timings::Dict{Symbol,Vector{Float64}})::Nothing
    phases = (
        (:plot_load,  lp.load_obs),
        (:plot_setup, lp.setup_obs),
        (:plot_data,  lp.data_obs),
        (:plot_draw,  lp.total_obs),
    )
    # All phase vectors share the same length as plot_load (or 0 when empty)
    ref_data = get(timings, :plot_load, Float64[])
    n = length(ref_data)
    if n == 0
        for (_, obs) in phases
            isempty(obs[]) || (obs[] = Float32[])
        end
        isempty(lp.timings_x_obs[]) || (lp.timings_x_obs[] = Float32[])
        return nothing
    end
    lp.timings_x_obs[] = Float32.(1:n)
    for (key, obs) in phases
        data = get(timings, key, Float64[])
        new_vals = isempty(data) ? zeros(Float32, n) : Float32.(data)
        obs[] = new_vals
    end
    return nothing
end

"""
Sample the current workspace's build counters into the ring buffers and update the build
figure Observables.

Called every frame from `_render_live_plots_tab!`. The ring buffers are the source of truth;
Observables are updated only when a new sample is appended or the buffers change.
"""
const _BUILD_SAMPLE_INTERVAL_NS = UInt64(250_000_000)  # 250 ms

function _sample_build_progress!(lp::LivePlotsState, workspace)::Nothing
    now = time_ns()
    now - lp.last_build_sample_ns < _BUILD_SAMPLE_INTERVAL_NS && return nothing
    lp.last_build_sample_ns = now

    scan_pct       = 0.0f0
    analysis_pct   = 0.0f0
    interp_ms      = 0.0f0
    processed_ms   = 0.0f0
    rss_mib        = Float32(Sys.maxrss()) / (1024.0f0^2)

    if workspace isa Workspace.Workspace
        progress = workspace.scan.progress
        total_si = progress.total_source_items
        scan_pct = total_si > 0 ?
            Float32(100 * progress.processed_source_items / total_si) : 0.0f0

        completed, total_proc = lock(workspace.processing.lock) do
            (workspace.processing.completed, workspace.processing.total)
        end
        analysis_pct = total_proc > 0 ?
            Float32(100 * completed / total_proc) : 0.0f0

        metrics = workspace.metrics
        interp_ms     = Float32(metrics.interpreted_write_ns[] / 1e6)
        processed_ms  = Float32(metrics.processed_write_ns[]  / 1e6)
    end

    cap = lp.capacity
    _ring_push!(lp.scan_pct_buf,          scan_pct,      cap)
    _ring_push!(lp.analysis_pct_buf,      analysis_pct,  cap)
    _ring_push!(lp.interp_write_ms_buf,   interp_ms,     cap)
    _ring_push!(lp.processed_write_ms_buf, processed_ms, cap)
    _ring_push!(lp.rss_mib_buf,           rss_mib,       cap)

    n = length(lp.scan_pct_buf)
    lp.build_x_obs[]          = Float32.(1:n)
    lp.scan_obs[]              = copy(lp.scan_pct_buf)
    lp.analysis_obs[]          = copy(lp.analysis_pct_buf)
    lp.interp_write_obs[]      = copy(lp.interp_write_ms_buf)
    lp.processed_write_obs[]   = copy(lp.processed_write_ms_buf)
    lp.rss_obs[]               = copy(lp.rss_mib_buf)
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
    # Write companion CSV
    if !isempty(series)
        csv_path = splitext(path)[1] * "_data.csv"
        try
            open(csv_path, "w") do io
                println(io, join(first.(series), ","))
                n = length(last(first(series)))
                for i in 1:n
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
Render the Live tab: two embedded Makie figures (timings + build progress) with Export buttons.

Figures are created once per `PerformanceState` lifetime and updated in place each frame via
Observables. `MakieFigure` polls GLMakie for pending updates and renders only dirty frames.
"""
function _render_live_plots_tab!(
    state::BrowserState;
    internal::Bool=false,
)::Nothing
    lp        = state.performance.live_plots
    workspace = state.workspace

    # Lazily create figures on first render
    if lp.timings_figure === nothing
        lp.timings_figure = _make_timings_figure(lp)
    end
    if lp.build_figure === nothing
        lp.build_figure = _make_build_figure(lp)
    end

    # Update timings Observables from the existing bounded vectors (no extra copy needed)
    _update_live_timings!(lp, state.performance.timings)

    # Sample build progress into ring buffers and push to Observables
    _sample_build_progress!(lp, workspace)

    if !internal
        ig.TextUnformatted("Plot redraw phase timings")
        ig.SameLine()
        if ig.Button("Export##live_timings_export")
            series = Pair{String,Vector{Float32}}[
                "load_ms"  => copy(lp.load_obs[]),
                "setup_ms" => copy(lp.setup_obs[]),
                "data_ms"  => copy(lp.data_obs[]),
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
        return nothing
    end

    # ---- Build progress figure ----
    ig.TextUnformatted("Build progress (scan %, analysis %, writer ms, peak RSS MiB)")
    ig.SameLine()
    if ig.Button("Export##live_build_export")
        series = Pair{String,Vector{Float32}}[
            "scan_pct"        => copy(lp.scan_obs[]),
            "analysis_pct"    => copy(lp.analysis_obs[]),
            "interp_write_ms" => copy(lp.interp_write_obs[]),
            "proc_write_ms"   => copy(lp.processed_write_obs[]),
            "peak_rss_mib"    => copy(lp.rss_obs[]),
        ]
        lp.build_export_error = _export_live_figure(
            lp.build_figure, "perf-build.png", series)
    end
    isempty(lp.build_export_error) || ig.TextWrapped(lp.build_export_error)

    MakieFigure("perf_live_build", lp.build_figure; auto_resize_x=true, auto_resize_y=false)

    return nothing
end

"""Run one internal profiler UI action and retain an actionable error."""
function _profile_action!(profiler::Profiling.ProfileSession, action::Function)::Nothing
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
    filtered = Profiling.ProfileEvent[
        event for event in available
        if (category_filter === :all || event.category === category_filter) &&
           (operation_filter === :all || event.operation === operation_filter)
    ]
    recent = length(filtered) <= 200 ? filtered : filtered[end-199:end]
    if !isempty(recent) && _begin_perf_table("internal_events", 6, 220.0f0)
        for name in ("Category", "Operation", "ms", "Thread", "Status", "Batch")
            ig.TableSetupColumn(name)
        end
        ig.TableHeadersRow()
        for event in Iterators.reverse(recent)
            ig.TableNextRow()
            for value in (
                String(event.category), String(event.operation),
                @sprintf("%.3f", event.duration_ns / 1e6),
                string(event.start_thread), String(event.status),
                string(event.attributes.batch_size),
            )
                ig.TableNextColumn(); _table_text(value)
            end
        end
        ig.EndTable()
    end

    ig.Separator()
    _render_live_plots_tab!(state; internal=true)
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
