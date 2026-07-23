using Printf: @sprintf
using Statistics: mean
import CImGui as ig
import DataBrowserCore.Workspace

# ---------------------------------------------------------------------------
# Performance window
#
# A lightweight, always-on diagnostic driven entirely by data the render loop
# already collects: the per-panel `_time!` timings/allocations in
# `state.performance`, and live workspace pipeline counters. It has no dependency
# on the (dev-only) instrumentation profiler — for deep engine timing use
# `@timed_dbg` + DataBrowserProfiling instead (see docs/profiling.md).
# ---------------------------------------------------------------------------

@inline function _table_text(s::AbstractString)::Nothing
    ig.TextUnformatted(s)
    return nothing
end

# ---------------------------------------------------------------------------
# Throughput sampling
#
# Cheap, always-readable engine counters sampled at ~4 Hz into fixed-cap ring
# buffers. No async timing and no new instrumentation — just differences of
# counters the pipeline already maintains.
# ---------------------------------------------------------------------------

"""Append `v` to an ordered history ring, dropping the oldest sample past `PERF_HISTORY_CAP`."""
@inline function _push_capped!(buf::Vector{Float32}, v::Real)::Nothing
    push!(buf, Float32(v))
    length(buf) > PERF_HISTORY_CAP && popfirst!(buf)
    return nothing
end

"""
Record one throughput sample into `h`, throttled to a 0.25s minimum interval.

Pure over its inputs (no workspace/ImGui) so it is unit-testable. The first admitted call only
seeds the `last_*` baselines (no data point, so rates never spike on open). Returns `true` when the
sample was admitted (baselines advanced), `false` when throttled.
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

"""Render frame rate, per-panel callback timings, and OpenGL strings (app overhead)."""
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
Render the Performance window: app render-loop overhead (Frames) and, when a
workspace is open, its live pipeline counters (Workspace). Toggled from the Debug
menu. All figures come from data the render loop already collects.
"""
function render_perf_window(state::BrowserState)::Nothing
    state.show_performance_window || return nothing

    if ig.Begin("Performance###perf_window")
        workspace = state.workspace
        if ig.BeginTabBar("##perf_tabs")
            if ig.BeginTabItem("Frames")
                ig.Spacing()
                _render_frames_tab(state)
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
