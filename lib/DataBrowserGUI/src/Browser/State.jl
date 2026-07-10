using DataBrowserAnnotations
import CImGui as ig

import DataBrowserCore.Workspace
import DataBrowserProfiling as Profiling
using DataBrowserAPI: KindProfileRow, Project, SourceProfileRow
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

"""
Ring-buffer state for the Performance window sparklines.

Timing buffers hold rolling plot-redraw phase durations (load / setup / draw-data / total).
Build buffers are sampled while the window is open; throughput is a finite difference of the
completed-work counter between ticks.
"""
Base.@kwdef mutable struct LivePlotsState
    capacity::Int = 600

    load_buf::Vector{Float32} = Float32[]
    setup_buf::Vector{Float32} = Float32[]
    data_buf::Vector{Float32} = Float32[]
    total_buf::Vector{Float32} = Float32[]
    timings_seen::Int = -1
    timings_export_error::String = ""

    elapsed_buf::Vector{Float32} = Float32[]
    throughput_buf::Vector{Float32} = Float32[]
    active_buf::Vector{Float32} = Float32[]
    pending_buf::Vector{Float32} = Float32[]
    rss_buf::Vector{Float32} = Float32[]
    build_export_error::String = ""

    t0_ns::UInt64 = UInt64(0)
    last_sample_ns::UInt64 = UInt64(0)
    last_elapsed_s::Float64 = 0.0
    last_completed::Int = 0
end

"""Counters and samples shown by the performance window."""
Base.@kwdef mutable struct PerformanceState
    frame::Int = 0
    gl_info::Dict{Symbol,String} = Dict{Symbol,String}()
    node_count::Int = 0
    item_rows_visible::Int = 0
    item_rows_rendered::Int = 0
    memory_start_rss::Union{Nothing,Int64} = nothing
    memory_peak_rss::Int64 = 0
    timings::Dict{Symbol,Vector{Float64}} = Dict{Symbol,Vector{Float64}}()
    allocations::Dict{Symbol,Vector{Int}} = Dict{Symbol,Vector{Int}}()
    scan_profile_project::Union{Nothing,Project} = nothing
    scan_profile_refresh_at::Float64 = 0.0
    scan_kind_rows::Vector{KindProfileRow} = KindProfileRow[]
    scan_source_rows::Vector{SourceProfileRow} = SourceProfileRow[]
    scan_kind_grid::DataGridState = DataGridState()
    scan_source_grid::DataGridState = DataGridState()
    plot_sampling_profile::Union{Nothing,Profiling.SamplingProfile} = nothing
    profile_category_filter::Symbol = :all
    profile_operation_filter::Symbol = :all
    live_plots::LivePlotsState = LivePlotsState()
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
    expanded_collection_paths::Vector{String} = String[]
    tree_filter::String = ""
    item_filter::String = ""
    tree_filter_widget::Union{Nothing,Ptr{ig.lib.ImGuiTextFilter}} = nothing
    item_filter_widget::Union{Nothing,Ptr{ig.lib.ImGuiTextFilter}} = nothing
    reset_project_filters::Bool = false
    scroll_to_collection_path::Union{Nothing,String} = nothing
    scroll_to_item_id::Union{Nothing,String} = nothing
    show_bad::Bool = true
    tag_state::Union{Nothing,DataBrowserAnnotations.Tags.TagState} = nothing
    tag_state_error::String = ""
    show_project_window::Bool = false
    show_performance_window::Bool = false
    collection_metadata_modal::Bool = true
    modal_root_path::String = ""
    cache_rebuild_modal::Bool = false
    cache_rebuild_path::String = ""
    cache_rebuild_project::Union{Nothing,Project} = nothing
    cache_rebuild_error::String = ""
    shutdown_complete::Bool = false
    extensions::Vector{GuiExtension} = GuiExtension[]
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
