using Printf: @sprintf
import CImGui as ig
import DataBrowserCore.Workspace
import TimerOutputs

# Performance window: the MAIN_TIMER section tree (Timings), live throughput sparklines
# (Throughput), and workspace pipeline counters (Workspace). For deep multi-task engine timing
# use `@timed_dbg` + DataBrowserProfiling instead (docs/profiling.md).

@inline function _table_text(s::AbstractString)::Nothing
    ig.TextUnformatted(s)
    return nothing
end

"""Append `v` to an ordered history ring, dropping the oldest sample past `PERF_HISTORY_CAP`."""
@inline function _push_capped!(buf::Vector{Float32}, v::Real)::Nothing
    push!(buf, Float32(v))
    length(buf) > PERF_HISTORY_CAP && popfirst!(buf)
    return nothing
end

"""
Record one throughput sample into `h`, throttled to a 0.25s minimum interval. The first admitted
call seeds the `last_*` baselines only, so rates never spike on open. Returns whether it was admitted.
"""
function _record_throughput!(
    h::PerfHistory,
    now::Float64,
    completed::Int,
    active::Int,
    pending_rows::Int,
    discovered::Int,
    cached::Int,
)::Bool
    dt = now - h.last_sample_t
    dt >= 0.25 || return false
    if h.last_sample_t != 0.0
        _push_capped!(h.items_per_s, max(completed - h.last_completed, 0) / dt)
        _push_capped!(h.active, active)
        _push_capped!(h.pending_rows, pending_rows)
        _push_capped!(h.scan_per_s, max(discovered - h.last_discovered, 0) / dt)
        _push_capped!(h.cache_per_s, max(cached - h.last_cached, 0) / dt)
    end
    h.last_sample_t = now
    h.last_completed = completed
    h.last_discovered = discovered
    h.last_cached = cached
    return true
end

"""Sample the active workspace's throughput counters into `state.performance.history`."""
function _sample_throughput!(state::BrowserState)::Nothing
    workspace = state.workspace
    workspace isa Workspace.Workspace || return nothing
    completed, _total, active = Workspace.work_counts(workspace)
    _record_throughput!(
        state.performance.history,
        time(),
        completed,
        active,
        Workspace.cache_pending_counts(workspace.cache.db).rows,
        workspace.scan.discovered[],
        Workspace.cache_stage_summary(workspace.cache.db).cached_sources,
    )
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

"""Render frame stats, the MAIN_TIMER section tree, and the OpenGL strings."""
function _render_timings_tab(state::BrowserState)::Nothing
    performance = state.performance

    raw_io = ig.GetIO()
    fps = unsafe_load(raw_io.Framerate)
    if fps > 0
        ig.Text(@sprintf("FPS: %.1f   Frame: %.2f ms", fps, 1000 / fps))
    end
    ig.Text("Tree nodes rendered: $(performance.node_count)")
    ig.Text("Items visible / rendered: $(performance.item_rows_visible) / $(performance.item_rows_rendered)")
    ig.Spacing()

    if ig.Button("Reset")
        state.performance.reset_main_timer = true   # applied at frame top, outside any @timed section
    end
    ig.SameLine()
    if ig.Button("Copy")
        # Copy the unicode tree — the paste target isn't limited to ImGui's ASCII-only font.
        ig.SetClipboardText(sprint(show, MAIN_TIMER))
    end
    ig.SameLine()
    ig.TextDisabled("Cumulative since the first frame; ncalls counts frames.")
    ig.Spacing()

    if isempty(MAIN_TIMER.inner_timers)
        ig.TextDisabled("No sections recorded yet.")
    else
        # Reuse TimerOutputs' own tree rendering. The default ImGui font is fixed-width but ASCII-only
        # (no box-drawing glyphs, no micro sign), so use ascii rules and render microseconds as "us".
        text = sprint(io -> TimerOutputs.print_timer(io, MAIN_TIMER; linechars=:ascii))
        _table_text(replace(rstrip(text), 'μ' => 'u'))
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
# Throughput tab — item-throughput sparklines
# ---------------------------------------------------------------------------

"""Draw a labeled sparkline of `values` (Dear ImGui `PlotLines`, no ImPlot) with a current-value readout."""
function _sparkline(label::String, values::Vector{Float32}, unit::String)::Nothing
    current = isempty(values) ? 0.0f0 : values[end]
    ig.Text(@sprintf("%s: %.1f%s", label, current, unit))
    if !isempty(values)
        lo = min(0.0f0, minimum(values))
        hi = max(1.0f-6, maximum(values))
        ig.PlotLines("##$label", values, length(values), 0, C_NULL, lo, hi, (0.0f0, 40.0f0))
    end
    return nothing
end

"""Render live item-throughput sparklines. Samples on entry; needs an open workspace."""
function _render_throughput_tab!(state::BrowserState)::Nothing
    workspace = state.workspace
    if !(workspace isa Workspace.Workspace)
        ig.TextDisabled("Open a workspace to see throughput.")
        return nothing
    end
    _sample_throughput!(state)
    h = state.performance.history

    _sparkline("Items analyzed", h.items_per_s, " /s")
    ig.Separator()
    ig.TextDisabled("Backlog")
    _sparkline("Active tasks", h.active, "")
    _sparkline("Pending cache rows", h.pending_rows, "")
    ig.Separator()
    _sparkline("Scan discovery", h.scan_per_s, " /s")
    _sparkline("Cache flush", h.cache_per_s, " /s")
    return nothing
end

"""Render live pipeline counters and index footprint for the active workspace."""
function _render_workspace_tab!(workspace::Workspace.Workspace)::Nothing
    completed, total, active = Workspace.work_counts(workspace)
    staged = Workspace.cache_pending_counts(workspace.cache.db)
    ig.Text(@sprintf("Analysis: %d / %d done   %d active", completed, total, active))
    ig.Text("Scan: $(workspace.scan.state)   Cache: $(workspace.cache_state)")
    ig.Text(@sprintf("Cache write buffer: %d items / %d rows pending", staged.items, staged.rows))
    ig.Spacing()

    index = workspace.index
    ig.Text("Index: $(length(index.items)) items, $(length(index.collections.records)) collections")
    ig.Text("Metadata held: $(length(index.item_metadata))")
    return nothing
end

"""
Render the Performance window: the main-task section tree (Timings), live item throughput
(Throughput), and, when a workspace is open, its pipeline counters (Workspace). Toggled from the
Debug menu. All figures come from data the render loop already collects.
"""
function render_perf_window(state::BrowserState)::Nothing
    state.show_performance_window || return nothing

    if ig.Begin("Performance###perf_window", C_NULL, ig.ImGuiWindowFlags_HorizontalScrollbar)
        workspace = state.workspace
        if ig.BeginTabBar("##perf_tabs")
            if ig.BeginTabItem("Timings")
                ig.Spacing()
                _render_timings_tab(state)
                ig.EndTabItem()
            end
            if ig.BeginTabItem("Throughput")
                ig.Spacing()
                _render_throughput_tab!(state)
                ig.EndTabItem()
            end
            if workspace isa Workspace.Workspace && ig.BeginTabItem("Workspace")
                ig.Spacing()
                _render_workspace_tab!(workspace)
                ig.EndTabItem()
            end
            ig.EndTabBar()
        end
    end
    ig.End()
    return nothing
end
