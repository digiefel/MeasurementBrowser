using Annotations
using GLMakie: Figure, Observable
import CImGui as ig

import ..Workspace
using ..Projects: KindProfileRow, Project, SourceProfileRow
using ..TableInspector: TablePreview, InspectorTable
using ..Visualization: PlotKind

"""
Runtime state for one plot window.

The same type represents the main plot and detached plots so drawing, export, persistence, and
error handling follow one path.
"""
Base.@kwdef mutable struct PlotViewState
    id::String
    title::String
    live::Bool
    item_ids::Vector{String} = String[]
    kind::Union{Nothing,Symbol} = nothing
    plot_kind::Union{Nothing,Type{<:PlotKind}} = nothing
    figure::Union{Nothing,Figure} = nothing
    last_key::Union{Nothing,Tuple} = nothing
    error::String = ""
    export_error::String = ""
end

"""
All plot-related browser state.

Plot choices are stored as Julia types while the app runs. Persistence converts them to names only
when writing `measurementbrowser.toml`.
"""
Base.@kwdef mutable struct PlotState
    main::PlotViewState = PlotViewState(id="main", title="Plot Area", live=true)
    windows::Vector{PlotViewState} = PlotViewState[]
    kind_by_item::Dict{Symbol,DataType} = Dict{Symbol,DataType}()
    next_window_id::Int = 0
    runtime_warmed::Bool = false
    warmup_figure::Union{Nothing,Figure} = nothing
end

"""Saved tree controls from `measurementbrowser.toml`."""
Base.@kwdef struct PersistedTreeView
    expanded::Vector{String} = String[]
    selected::Vector{String} = String[]
    filter::String = ""
end

"""Saved item controls from `measurementbrowser.toml`."""
Base.@kwdef struct PersistedItemsView
    selected::Vector{String} = String[]
    filter::String = ""
end

"""Saved plot-window state from `measurementbrowser.toml`."""
Base.@kwdef struct PersistedPlotView
    id::String = ""
    title::String = ""
    plot_kind::String = ""
    live::Bool = false
    items::Vector{String} = String[]
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
    plot_kinds::Dict{String,String} = Dict{String,String}()
    main_plot::PersistedPlotView =
        PersistedPlotView(id="main", title="Plot Area", live=true)
    plot_windows::Vector{PersistedPlotView} = PersistedPlotView[]
end

"""One recently opened source root and its user preferences."""
Base.@kwdef struct RecentProject
    path::String
    project_preference::String = "auto"
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
    plot_selected_only::Bool = false
    show_provenance_column::Bool = false
    # current_kind drives the per-kind DataGrid table id so imgui.ini keys column widths per kind
    current_kind::Union{Nothing,Symbol} = nothing
    # Quick-plot state (driven from inspector_table)
    x_column::Int = 1
    y_column::Int = 2
    figure::Union{Nothing,Figure} = nothing
    plot_key::Union{Nothing,Tuple} = nothing
    plot_error::String = ""
    # Raw file-preview mode (secondary): file → DataGrid
    preview::Union{Nothing,TablePreview} = nothing
    file_grid::DataGridState = DataGridState()
    live::Bool = true
    path_buffer::Vector{UInt8} = fill(UInt8(0), TABLE_INSPECTOR_PATH_BUFFER_SIZE)
    error::String = ""
end

"""
Ring-buffer state and Makie figure handles for the live-plots tab in the Performance window.

Two figures are created once and updated in place via Observables:
- `timings_figure`: rolling plot-redraw phase durations (load / setup / draw-data / total).
- `build_figure`: a six-panel dashboard of the live build — progress, throughput, writer concurrency,
  per-item write cost, cumulative write time, and RSS — the same quantities `bench/visualize_build.jl`
  reconstructs from the CSVs, but computed live from the sampled counters.

Each line carries its own `_x`/`_obs` pair so series of different length never collide (the plot-error
path records only `:plot_draw`, so the phase vectors do diverge). `_buf` fields are bounded ring
buffers sampled once per poll tick. Derived rates (throughput, writers-busy, per-item write) are
finite differences of the cumulative counters between ticks, so `last_*` hold the previous snapshot.
"""
Base.@kwdef mutable struct LivePlotsState
    capacity::Int = 600

    # ---- timings figure: one independent (x, y) ring buffer per redraw phase ----
    timings_figure::Union{Nothing,Figure} = nothing
    load_x::Observable{Vector{Float32}}  = Observable(Float32[])
    load_obs::Observable{Vector{Float32}} = Observable(Float32[])
    setup_x::Observable{Vector{Float32}} = Observable(Float32[])
    setup_obs::Observable{Vector{Float32}} = Observable(Float32[])
    data_x::Observable{Vector{Float32}}  = Observable(Float32[])
    data_obs::Observable{Vector{Float32}} = Observable(Float32[])
    total_x::Observable{Vector{Float32}} = Observable(Float32[])
    total_obs::Observable{Vector{Float32}} = Observable(Float32[])
    timings_seen::Int = -1
    timings_export_error::String = ""

    # ---- build dashboard: shared elapsed-seconds x, one ring buffer per series ----
    build_figure::Union{Nothing,Figure} = nothing
    elapsed_buf::Vector{Float32}        = Float32[]
    scan_pct_buf::Vector{Float32}       = Float32[]
    analysis_pct_buf::Vector{Float32}   = Float32[]
    throughput_buf::Vector{Float32}     = Float32[]   # analysed items / s
    writers_busy_buf::Vector{Float32}   = Float32[]   # avg writers busy concurrently
    per_item_write_buf::Vector{Float32} = Float32[]   # processed write ms / item
    interp_cum_buf::Vector{Float32}     = Float32[]   # cumulative interp-write s
    processed_cum_buf::Vector{Float32}  = Float32[]   # cumulative processed-write s
    stats_cum_buf::Vector{Float32}      = Float32[]   # cumulative stats-write s
    rss_buf::Vector{Float32}            = Float32[]   # RSS GiB

    elapsed_obs::Observable{Vector{Float32}}      = Observable(Float32[])
    scan_pct_obs::Observable{Vector{Float32}}     = Observable(Float32[])
    analysis_pct_obs::Observable{Vector{Float32}} = Observable(Float32[])
    throughput_obs::Observable{Vector{Float32}}   = Observable(Float32[])
    writers_busy_obs::Observable{Vector{Float32}} = Observable(Float32[])
    per_item_write_obs::Observable{Vector{Float32}} = Observable(Float32[])
    interp_cum_obs::Observable{Vector{Float32}}   = Observable(Float32[])
    processed_cum_obs::Observable{Vector{Float32}} = Observable(Float32[])
    stats_cum_obs::Observable{Vector{Float32}}    = Observable(Float32[])
    rss_obs::Observable{Vector{Float32}}          = Observable(Float32[])
    build_export_error::String = ""

    # ---- sampling state (finite-difference bases + throttle) ----
    t0_ns::UInt64 = UInt64(0)
    last_sample_ns::UInt64 = UInt64(0)
    last_elapsed_s::Float64 = 0.0
    last_completed::Int = 0
    last_processed_write_ns::Int64 = 0
    last_writer_busy_ns::Int64 = 0
end

"""Counters and samples shown by the performance window."""
Base.@kwdef mutable struct PerformanceState
    frame::Int = 0
    gl_info::Dict{Symbol,String} = Dict{Symbol,String}()
    node_count::Int = 0
    item_rows_visible::Int = 0
    item_rows_rendered::Int = 0
    memory_start_rss_kb::Union{Nothing,Int} = nothing
    memory_peak_rss_kb::Union{Nothing,Int} = nothing
    memory_start_read_bytes::Union{Nothing,Int} = nothing
    timings::Dict{Symbol,Vector{Float64}} = Dict{Symbol,Vector{Float64}}()
    allocations::Dict{Symbol,Vector{Int}} = Dict{Symbol,Vector{Int}}()
    scan_profile_project::Union{Nothing,Project} = nothing
    scan_profile_refresh_at::Float64 = 0.0
    scan_kind_rows::Vector{KindProfileRow} = KindProfileRow[]
    scan_source_rows::Vector{SourceProfileRow} = SourceProfileRow[]
    scan_kind_grid::DataGridState = DataGridState()
    scan_source_grid::DataGridState = DataGridState()
    sampling_grid::DataGridState = DataGridState()
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
    plots::PlotState = PlotState()
    table_inspector::TableInspectorState = TableInspectorState()
    performance::PerformanceState = PerformanceState()
    project_preference::String = "auto"
    recent_projects::Vector{RecentProject} = RecentProject[]
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
    tag_state::Union{Nothing,Annotations.Tags.TagState} = nothing
    tag_state_error::String = ""
    show_project_window::Bool = false
    show_performance_window::Bool = false
    # When set, the next plot redraw is wrapped in the sampling profiler + cache timers and the full
    # breakdown is printed to the console. Set it via the toolbar button or `request_plot_profile!`.
    profile_next_plot::Bool = false
    collection_parameters_modal::Bool = true
    modal_root_path::String = ""
    shutdown_complete::Bool = false
end
