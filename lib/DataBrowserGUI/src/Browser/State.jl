using DataBrowserAnnotations
import CImGui as ig

import DataBrowserCore.Workspace
using DataBrowserAPI: Project
using DataBrowserSources
using DataBrowserCore: InspectorTable

"""Saved tree controls from `databrowser.toml`."""
Base.@kwdef struct PersistedTreeView
    expanded::Vector{String} = String[]
    selected::Vector{String} = String[]
    filter::String = ""
end

"""Saved item controls from `databrowser.toml`."""
Base.@kwdef struct PersistedItemsView
    selected::Vector{String} = String[]
    filter::String = ""
end

"""
The browser state written inside one project directory.

Only stable ids and text are persisted. Runtime objects, loaded data, figures, and background work
remain in memory.
"""
Base.@kwdef struct PersistedProjectView
    project::String = ""
    tree::PersistedTreeView = PersistedTreeView()
    items::PersistedItemsView = PersistedItemsView()
    extensions::Dict{String,Dict{String,Any}} = Dict{String,Dict{String,Any}}()
end

const TABLE_INSPECTOR_PATH_BUFFER_SIZE = 1024

"""
Transient selection and scroll state for a DataGrid widget.

Owned by the consumer of the grid (e.g. TableInspectorState); passed by reference each frame.
Column widths are persisted across restarts via the imgui.ini file, keyed by the grid `id`.
Row mode uses `selected_rows`. Cell mode stores the opposite corners of one rectangular selection in
`cell_anchor` and `cell_focus`. Use `scroll_to_row` to request a programmatic scroll on the next frame.
"""
Base.@kwdef mutable struct DataGridState
    selected_rows::Vector{Int}           = Int[]
    scroll_to_row::Union{Nothing,Int}    = nothing
    cell_anchor::Union{Nothing,Tuple{Int,Int}} = nothing
    cell_focus::Union{Nothing,Tuple{Int,Int}} = nothing
    dragging_cells::Bool                 = false
    focused::Bool                        = false
end

"""State for the generic table-inspection window."""
Base.@kwdef mutable struct TableInspectorState
    visible::Bool = false
    # Item-data view (primary mode)
    inspector_table::Union{Nothing,InspectorTable} = nothing
    inspector_warnings::Vector{String} = String[]
    inspector_key::Union{Nothing,Tuple} = nothing  # (item_ids..., show_provenance)
    grid::DataGridState = DataGridState()
    show_provenance_column::Bool = false
    # current_kind drives the per-kind DataGrid table id so imgui.ini keys column widths per kind
    current_kind::Union{Nothing,Symbol} = nothing
    # Raw file-preview mode (secondary): file → DataGrid
    preview::Union{Nothing,DataBrowserSources.TabularFileSource} = nothing
    file_grid::DataGridState = DataGridState()
    live::Bool = true
    path_buffer::Vector{UInt8} = fill(UInt8(0), TABLE_INSPECTOR_PATH_BUFFER_SIZE)
    error::String = ""
end

"""Newest-last cap for the throughput history ring buffers (~1 minute at the 0.25s sample rate)."""
const PERF_HISTORY_CAP = 240

"""
Bounded history of item-throughput samples for the Performance window's Throughput tab.

Each vector is an ordered ring buffer (oldest first, newest last, capped at `PERF_HISTORY_CAP`).
Rates are per-second deltas between samples; `active`/`pending_rows` are instantaneous levels.
The `last_*` fields hold the previous sample's raw counters so the next rate can be differenced.
"""
Base.@kwdef mutable struct PerfHistory
    items_per_s::Vector{Float32}  = Float32[]
    active::Vector{Float32}       = Float32[]   # background tasks in flight  ┐ backlog
    pending_rows::Vector{Float32} = Float32[]   # cache write queue depth     ┘
    scan_per_s::Vector{Float32}   = Float32[]
    cache_per_s::Vector{Float32}  = Float32[]
    last_sample_t::Float64 = 0.0
    last_completed::Int = 0
    last_discovered::Int = 0
    last_cached::Int = 0
end

"""Frame counters and item-throughput history gathered during the render loop."""
Base.@kwdef mutable struct PerformanceState
    frame::Int = 0
    """Monotonic `time()` when the first non-blank frame was submitted, or `NaN` until then."""
    first_frame_at::Float64 = NaN
    gl_info::Dict{Symbol,String} = Dict{Symbol,String}()
    node_count::Int = 0
    item_rows_visible::Int = 0
    item_rows_rendered::Int = 0
    history::PerfHistory = PerfHistory()
    # Removed in the Performance-window rewrite; unwritten since the @timed migration.
    timings::Dict{Symbol,Vector{Float64}} = Dict{Symbol,Vector{Float64}}()
    allocations::Dict{Symbol,Vector{Int}} = Dict{Symbol,Vector{Int}}()
end

"""
Runtime state owned by the browser.

The workspace owns items, selection, cache state, data, and background work. This type owns
only controls, windows, local persistence, annotations, and rendering state.
"""
Base.@kwdef mutable struct BrowserState
    workspace::Union{Nothing,Workspace.Workspace} = nothing
    project_locked::Bool = false
    table_inspector::TableInspectorState = TableInspectorState()
    performance::PerformanceState = PerformanceState()
    project_preference::String = "auto"
    saved_project_view::PersistedProjectView = PersistedProjectView()
    expanded_collection_ids::Vector{String} = String[]
    tree_filter::String = ""
    item_filter::String = ""
    tree_filter_widget::Union{Nothing,Ptr{ig.lib.ImGuiTextFilter}} = nothing
    item_filter_widget::Union{Nothing,Ptr{ig.lib.ImGuiTextFilter}} = nothing
    reset_project_filters::Bool = false
    scroll_to_collection_id::Union{Nothing,String} = nothing
    scroll_to_item_id::Union{Nothing,String} = nothing
    show_bad::Bool = true
    tag_state::Union{Nothing,DataBrowserAnnotations.Tags.TagState} = nothing
    tag_state_error::String = ""
    show_project_window::Bool = false
    show_performance_window::Bool = false
    show_imgui_metrics::Bool = false
    show_imgui_debug_log::Bool = false
    show_imgui_id_stack::Bool = false
    show_imgui_style_editor::Bool = false
    show_imgui_user_guide::Bool = false
    show_imgui_about::Bool = false
    show_imgui_demo::Bool = false
    show_implot_metrics::Bool = false
    show_implot_style_editor::Bool = false
    show_implot_user_guide::Bool = false
    show_implot_demo::Bool = false
    implot_context::Ptr{ig.lib.ImPlotContext} = C_NULL
    collection_metadata_modal::Bool = true
    modal_root_path::String = ""
    cache_rebuild_modal::Bool = false
    cache_rebuild_path::String = ""
    cache_rebuild_project::Union{Nothing,Project} = nothing
    cache_rebuild_error::String = ""
    shutdown_complete::Bool = false
    """Set by `close_browser!` so the render loop exits without a GLFW close click."""
    exit_requested::Bool = false
    extensions::Vector{GuiExtension} = GuiExtension[]
end

"""Non-blocking `open_browser` handle: the render `task` plus its `BrowserState`."""
struct BrowserSession
    task::Task
    state::BrowserState
end

"""Call `shutdown!` on every extension instance."""
function _shutdown_extensions!(state::BrowserState)::Nothing
    for ext in state.extensions
        shutdown!(ext, state)
    end
    return nothing
end

"""Call `reset!` on every extension instance."""
function _reset_extensions!(state::BrowserState)::Nothing
    for ext in state.extensions
        reset!(ext, state)
    end
    return nothing
end

"""True when every extension reports ready (warmup gating)."""
function _extensions_ready(state::BrowserState)::Bool
    return all(is_ready(ext, state) for ext in state.extensions)
end
