module Visualization

using GLMakie: Figure
using InteractiveUtils: subtypes

import ..ItemIndex: ItemRecord
import ..Projects: project_name
import ..Workspace

"""
Base type for visualizers selectable by the browser and Julia API.

`PlotKind` is an internal identity used by the engine to dispatch plotting. Projects do not define
their own subtypes; they call `register_plot!`, and the engine represents each registered plot as a
`RegistryPlot{kind,label}` (below).
"""
abstract type PlotKind end

"""
Internal plot identity for a plot registered with `register_plot!`, parameterized by the measurement
kind it applies to and the label that distinguishes plots for that kind.
"""
struct RegistryPlot{Kind,Label} <: PlotKind end

"""The measurement-kind symbol a `RegistryPlot` draws."""
plot_kind_symbol(::Type{RegistryPlot{Kind,Label}}) where {Kind,Label} = Kind

"""Stable name used to persist a plot choice in `measurementbrowser.toml`."""
plot_kind_name(plot_kind::Type{<:PlotKind})::String = String(nameof(plot_kind))
plot_kind_name(::Type{RegistryPlot{Kind,Label}}) where {Kind,Label} =
    string(Kind, "::", Label)

"""Resolve a persisted plot name back to its internal identity (a registered plot kind)."""
function plot_kind_from_name(name::AbstractString)::Union{Nothing,Type{<:PlotKind}}
    parts = split(name, "::"; limit=2)
    length(parts) == 2 || return nothing
    kind, label = parts
    return RegistryPlot{Symbol(kind),Symbol(label)}
end

"""Create the figure layout required by a visualizer."""
function setup_plot(
    workspace::Workspace.Workspace,
    plot_kind::Type{<:PlotKind},
    measurements::Vector{ItemRecord},
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
    measurements::Vector{ItemRecord},
    figure::Figure,
)::Nothing
    error(
        "No plot_data! implementation for $(project_name(workspace.project)) " *
        "plot kind '$plot_kind'",
    )
end

# --- Registered-plot bridge ------------------------------------------------
# A RegistryPlot{kind,label} dispatches the engine's plot calls to the project's registered setup/draw
# callbacks. `project.plots` is the recipe table filled by register_plot!.

"""
Build the figure for a registered plot by running its `setup` callback.

The package resolves processed data (through its cache) before invoking the recipe, so `setup` can
size the figure layout to the data without calling `process_item_data` itself.
"""
function setup_plot(
    workspace::Workspace.Workspace,
    ::Type{RegistryPlot{Kind,Label}},
    measurements::Vector{ItemRecord},
)::Figure where {Kind,Label}
    recipe = workspace.project.plots[Kind][String(Label)]
    processed = Workspace.process_item_data(workspace, measurements)
    return recipe.setup(workspace, measurements, processed)::Figure
end

"""
Draw a registered plot into its figure by running its `draw` callback.

The package resolves processed data (through its cache) before invoking the recipe, so `draw`
callbacks receive the processed `DataFrame`s directly and never call `process_item_data`
themselves.
"""
function plot_data!(
    workspace::Workspace.Workspace,
    ::Type{RegistryPlot{Kind,Label}},
    measurements::Vector{ItemRecord},
    figure::Figure,
)::Nothing where {Kind,Label}
    recipe = workspace.project.plots[Kind][String(Label)]
    processed = Workspace.process_item_data(workspace, measurements)
    recipe.draw(workspace, measurements, processed, figure)
    return nothing
end

"""Every plot registered for a measurement kind, sorted by label."""
function registered_plot_kinds(project, kind::Symbol)::Vector{Type{<:PlotKind}}
    recipes = get(project.plots, kind, nothing)
    recipes === nothing && return Type{<:PlotKind}[]
    return Type{<:PlotKind}[
        RegistryPlot{kind,Symbol(label)}
        for label in sort!(collect(keys(recipes)))
    ]
end

"""Human label for a registered plot kind, taken from its recipe."""
function plot_kind_label(project, ::Type{RegistryPlot{Kind,Label}})::String where {Kind,Label}
    return project.plots[Kind][String(Label)].label
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
    measurements::Vector{ItemRecord},
    loaded;
    kwargs...,
)
    error("Debug plots are not implemented for this project")
end

end
