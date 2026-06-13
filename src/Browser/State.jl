using Annotations
using GLMakie: Figure
import CImGui as ig

import ..Workspace
using ..MeasurementIndex: MeasurementInfo
using ..TableInspector: TablePreview
using ..Visualization: PlotKind
using ..MeasurementBrowser:
    NamedMeasurementGroup,
    _FigureScriptFactIndex

"""
Runtime state for one plot window.

The same type represents the main plot and detached plots so drawing, export, persistence, and
error handling follow one path.
"""
Base.@kwdef mutable struct PlotViewState
    id::String
    title::String
    live::Bool
    measurement_ids::Vector{String} = String[]
    measurement_kind::Union{Nothing,Symbol} = nothing
    plot_kind::Union{Nothing,Type{<:PlotKind}} = nothing
    figure::Union{Nothing,Figure} = nothing
    last_key::Union{Nothing,Tuple} = nothing
    error::String = ""
    export_error::String = ""
    debug::Bool = false
end

"""
All plot-related browser state.

Plot choices are stored as Julia types while the app runs. Persistence converts them to names only
when writing `measurementbrowser.toml`.
"""
Base.@kwdef mutable struct PlotState
    main::PlotViewState = PlotViewState(id="main", title="Plot Area", live=true)
    windows::Vector{PlotViewState} = PlotViewState[]
    kind_by_measurement::Dict{Symbol,DataType} = Dict{Symbol,DataType}()
    debug::Bool = false
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

"""Saved measurement controls from `measurementbrowser.toml`."""
Base.@kwdef struct PersistedMeasurementsView
    selected::Vector{String} = String[]
    filter::String = ""
end

"""Saved plot-window state from `measurementbrowser.toml`."""
Base.@kwdef struct PersistedPlotView
    id::String = ""
    title::String = ""
    plot_kind::String = ""
    live::Bool = false
    measurements::Vector{String} = String[]
end

"""
The browser state written inside one project directory.

Only stable ids and text are persisted. Runtime objects, loaded data, figures, and background work
remain in memory.
"""
Base.@kwdef struct PersistedProjectView
    project::String = ""
    tree::PersistedTreeView = PersistedTreeView()
    measurements::PersistedMeasurementsView = PersistedMeasurementsView()
    plot_kinds::Dict{String,String} = Dict{String,String}()
    main_plot::PersistedPlotView =
        PersistedPlotView(id="main", title="Plot Area", live=true)
    plot_windows::Vector{PersistedPlotView} = PersistedPlotView[]
end

"""One recently opened source root and its user preferences."""
Base.@kwdef struct RecentProject
    path::String
    project_preference::String = "auto"
    figure_script_output_dir::String = ""
end

"""
Temporary state for the figure-script interface.

This state remains isolated because the Workflow phase will replace figure scripts without changing
the rest of the browser.
"""
Base.@kwdef mutable struct FigureScriptState
    visible::Bool = false
    root_path::String = ""
    output_dir_buffer::Vector{UInt8} = fill(UInt8(0), FIGURE_OUTPUT_DIR_BUFFER_SIZE)
    script_name_buffer::Vector{UInt8} = fill(UInt8(0), FIGURE_SCRIPT_NAME_BUFFER_SIZE)
    group_name_buffer::Vector{UInt8} = fill(UInt8(0), FIGURE_GROUP_NAME_BUFFER_SIZE)
    groups::Vector{NamedMeasurementGroup} = NamedMeasurementGroup[]
    group_matches::Dict{String,Vector{MeasurementInfo}} =
        Dict{String,Vector{MeasurementInfo}}()
    group_matches_valid::Bool = false
    fact_index::Union{Nothing,_FigureScriptFactIndex} = nothing
    fact_index_valid::Bool = false
    selected_group::Int = 0
    error::String = ""
    status::String = ""
    overwrite_confirm::String = ""
    job_state::Symbol = :idle
    job_sequence::Int = 0
    active_job_id::Int = 0
    job_events::Union{Nothing,Channel{NamedTuple}} = nothing
    job_profiles::Vector{NamedTuple} = NamedTuple[]
end

const TABLE_INSPECTOR_PATH_BUFFER_SIZE = 1024

"""State for the generic table-inspection window."""
Base.@kwdef mutable struct TableInspectorState
    visible::Bool = false
    live::Bool = true
    path_buffer::Vector{UInt8} = fill(UInt8(0), TABLE_INSPECTOR_PATH_BUFFER_SIZE)
    preview::Union{Nothing,TablePreview} = nothing
    error::String = ""
    x_column::Int = 1
    y_column::Int = 2
    figure::Union{Nothing,Figure} = nothing
    plot_key::Union{Nothing,Tuple} = nothing
    plot_error::String = ""
end

"""Counters and samples shown by the performance window."""
Base.@kwdef mutable struct PerformanceState
    frame::Int = 0
    gl_info::Dict{Symbol,String} = Dict{Symbol,String}()
    node_count::Int = 0
    measurement_rows_visible::Int = 0
    measurement_rows_rendered::Int = 0
    memory_start_rss_kb::Union{Nothing,Int} = nothing
    memory_peak_rss_kb::Union{Nothing,Int} = nothing
    memory_start_read_bytes::Union{Nothing,Int} = nothing
    timings::Dict{Symbol,Vector{Float64}} = Dict{Symbol,Vector{Float64}}()
    allocations::Dict{Symbol,Vector{Int}} = Dict{Symbol,Vector{Int}}()
end

"""
Runtime state owned by the browser.

The workspace owns measurements, selection, cache state, data, and background work. This type owns
only controls, windows, local persistence, annotations, and rendering state.
"""
Base.@kwdef mutable struct BrowserState
    workspace::Union{Nothing,Workspace.Workspace} = nothing
    project_locked::Bool = false
    plots::PlotState = PlotState()
    figure_scripts::FigureScriptState = FigureScriptState()
    table_inspector::TableInspectorState = TableInspectorState()
    performance::PerformanceState = PerformanceState()
    project_preference::String = "auto"
    recent_projects::Vector{RecentProject} = RecentProject[]
    saved_project_view::PersistedProjectView = PersistedProjectView()
    expanded_device_paths::Vector{String} = String[]
    tree_filter::String = ""
    measurement_filter::String = ""
    tree_filter_widget::Union{Nothing,Ptr{ig.lib.ImGuiTextFilter}} = nothing
    measurement_filter_widget::Union{Nothing,Ptr{ig.lib.ImGuiTextFilter}} = nothing
    reset_project_filters::Bool = false
    scroll_to_device_path::Union{Nothing,String} = nothing
    scroll_to_measurement_id::Union{Nothing,String} = nothing
    show_bad::Bool = true
    tag_state::Union{Nothing,Annotations.Tags.TagState} = nothing
    tag_state_error::String = ""
    show_project_window::Bool = false
    show_performance_window::Bool = false
    device_info_modal::Bool = true
    modal_root_path::String = ""
    shutdown_complete::Bool = false
end
const FIGURE_SCRIPT_NAME_BUFFER_SIZE = 192
const FIGURE_GROUP_NAME_BUFFER_SIZE = 192
const FIGURE_OUTPUT_DIR_BUFFER_SIZE = 512
