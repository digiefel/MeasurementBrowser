module Visualization

using GLMakie: Figure
using InteractiveUtils: subtypes

import ..MeasurementIndex: MeasurementInfo
import ..Projects: project_name
import ..Workspace

"""
Base type for visualizers selectable by the browser and Julia API.
"""
abstract type PlotKind end

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
