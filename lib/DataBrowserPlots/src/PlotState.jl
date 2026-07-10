using GLMakie: Figure

using DataBrowserAPI: PlotKind

"""Runtime state for one plot window."""
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

"""All plot-related browser state."""
Base.@kwdef mutable struct PlotState
    main::PlotViewState = PlotViewState(id="main", title="Plot Area", live=true)
    windows::Vector{PlotViewState} = PlotViewState[]
    kind_by_item::Dict{Symbol,DataType} = Dict{Symbol,DataType}()
    next_window_id::Int = 0
    runtime_warmed::Bool = false
    warmup_figure::Union{Nothing,Figure} = nothing
end

"""Saved plot-window state under `extensions.PlotsExtension` in `databrowser.toml`."""
Base.@kwdef struct PersistedPlotView
    id::String = ""
    title::String = ""
    plot_kind::String = ""
    live::Bool = false
    items::Vector{String} = String[]
end

"""State for the workspace-selection table plot window."""
Base.@kwdef mutable struct TablePlotState
    visible::Bool = false
    x_column::Int = 1
    y_column::Int = 2
    figure::Union{Nothing,Figure} = nothing
    plot_key::Union{Nothing,Tuple} = nothing
    plot_error::String = ""
end
