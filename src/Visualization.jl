module Visualization

using GLMakie: Figure
using InteractiveUtils: subtypes

import ..MeasurementIndex: MeasurementInfo
import ..Projects: project_name
import ..Workspace

"""
Base type for visualizers selectable by the browser and Julia API.

`PlotKind` is an internal identity used by the engine to dispatch plotting. Projects do not define
their own subtypes; they call `register_plot!`, and the engine represents each registered plot as a
`RegistryPlot{kind}` (below).
"""
abstract type PlotKind end

"""
Internal plot identity for a plot registered with `register_plot!`, parameterized by the measurement
kind it draws. One registered plot per measurement kind, so the kind alone identifies the plot.
"""
struct RegistryPlot{Kind} <: PlotKind end

"""The measurement-kind symbol a `RegistryPlot` draws."""
plot_kind_symbol(::Type{RegistryPlot{Kind}}) where {Kind} = Kind

"""Stable name used to persist a plot choice in `measurementbrowser.toml`."""
plot_kind_name(plot_kind::Type{<:PlotKind})::String = String(nameof(plot_kind))
plot_kind_name(::Type{RegistryPlot{Kind}}) where {Kind} = String(Kind)

"""Resolve a persisted plot name back to its internal identity (a registered plot kind)."""
plot_kind_from_name(name::AbstractString)::Type{<:PlotKind} = RegistryPlot{Symbol(name)}

"""Create the figure layout required by a visualizer."""
function setup_plot(
    workspace::Workspace.Workspace,
    plot_kind::Type{<:PlotKind},
    measurements::Vector{MeasurementInfo},
)::Figure
    error(
        "No setup_plot implementation for $(project_name(workspace.project)) " *
        "plot kind '$plot_kind'",
    )
end

"""Draw selected measurements into an existing figure."""
function plot_data!(
    workspace::Workspace.Workspace,
    plot_kind::Type{<:PlotKind},
    measurements::Vector{MeasurementInfo},
    figure::Figure,
)::Nothing
    error(
        "No plot_data! implementation for $(project_name(workspace.project)) " *
        "plot kind '$plot_kind'",
    )
end

# --- Registered-plot bridge ------------------------------------------------
# A RegistryPlot{kind} dispatches the engine's plot calls to the project's registered setup/draw
# callbacks. `project.plots` is the recipe table filled by register_plot!.

"""
Build the figure for a registered plot by running its `setup` callback.

The package resolves processed data (through its cache) before invoking the recipe, so `setup` can
size the figure layout to the data without calling `process_measurement_data` itself.
"""
function setup_plot(
    workspace::Workspace.Workspace,
    ::Type{RegistryPlot{Kind}},
    measurements::Vector{MeasurementInfo},
)::Figure where {Kind}
    recipe = get(workspace.project.plots, Kind, nothing)
    recipe === nothing && error("No plot registered for measurement kind :$Kind")
    processed = Workspace.process_measurement_data(workspace, measurements)
    return recipe.setup(workspace, measurements, processed)::Figure
end

"""
Draw a registered plot into its figure by running its `draw` callback.

The package resolves processed data (through its cache) before invoking the recipe, so `draw`
callbacks receive the processed `DataFrame`s directly and never call `process_measurement_data`
themselves.
"""
function plot_data!(
    workspace::Workspace.Workspace,
    ::Type{RegistryPlot{Kind}},
    measurements::Vector{MeasurementInfo},
    figure::Figure,
)::Nothing where {Kind}
    recipe = get(workspace.project.plots, Kind, nothing)
    recipe === nothing && error("No plot registered for measurement kind :$Kind")
    processed = Workspace.process_measurement_data(workspace, measurements)
    recipe.draw(workspace, measurements, processed, figure)
    return nothing
end

"""The plot kind registered for a measurement kind, or `nothing` if none is registered."""
registered_plot_kind(project, measurement_kind::Symbol)::Union{Nothing,Type{<:PlotKind}} =
    haskey(project.plots, measurement_kind) ? RegistryPlot{measurement_kind} : nothing

"""Human label for a registered plot kind, taken from its recipe."""
function plot_kind_label(project, ::Type{RegistryPlot{Kind}})::String where {Kind}
    recipe = get(project.plots, Kind, nothing)
    return recipe === nothing ? string(Kind) : recipe.label
end

"""
Return every loaded concrete visualizer type.

Visualizers are discovered from Julia's type hierarchy, so external project packages require no
registration with MeasurementBrowser.
"""
function plot_kinds()::Vector{Type{<:PlotKind}}
    kinds = Type{<:PlotKind}[]
    pending = Type[PlotKind]
    while !isempty(pending)
        for child in subtypes(pop!(pending))
            isabstracttype(child) ? push!(pending, child) : push!(kinds, child)
        end
    end
    sort!(kinds; by=kind -> String(nameof(kind)))
    return kinds
end

"""Draw source data directly while bypassing package caches and processed data."""
function debug_plot(
    ::Workspace.Workspace,
    measurements::Vector{MeasurementInfo},
    loaded;
    kwargs...,
)
    error("Debug plots are not implemented for this project")
end

end
