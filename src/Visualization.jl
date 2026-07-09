module Visualization

using GLMakie: Figure

import DataBrowserAPI
import DataBrowserAPI: AbstractDataItem, PlotKind, RegisteredPlot, plot_data!, setup_plot
import ..ItemIndex: ItemRecord
import DataBrowserProfiling as Profiling
import ..Workspace

"""Materialize records and create the figure layout required by a visualizer."""
function setup_plot(
    workspace::Workspace.Workspace,
    plot_kind::Type{<:PlotKind},
    records::Vector{ItemRecord},
)::Figure
    return setup_plot(
        workspace, plot_kind, Workspace.materialize_items(workspace, records))
end

"""Materialize records and draw their items into an existing figure."""
function plot_data!(
    workspace::Workspace.Workspace,
    plot_kind::Type{<:PlotKind},
    records::Vector{ItemRecord},
    figure::Figure,
)::Nothing
    plot_data!(
        workspace, plot_kind, Workspace.materialize_items(workspace, records), figure)
    return nothing
end

"""
Build the figure for a registered plot by running its `setup` callback.

The package materializes the loaded `items` (each with `item.data`) before invoking the recipe, so
`setup` sizes the figure layout to the data without resolving it itself.
"""
function setup_plot(
    workspace::Workspace.Workspace,
    ::Type{RegisteredPlot{Kind,Label}},
    items::Vector{<:AbstractDataItem},
)::Figure where {Kind,Label}
    recipe = workspace.project.plots[Kind][String(Label)]
    return Profiling.@profile_span workspace.profiler :visualization :setup_plot Profiling.ProfileAttributes(
        kind=Kind,
        items=Int64(length(items)),
    ) begin
        recipe.setup(workspace, items)::Figure
    end
end

"""
Draw a registered plot into its figure by running its `draw` callback.

The package materializes the loaded `items` (each with `item.data`) before invoking the recipe, so
`draw` reads `item.data` directly and never resolves data itself.
"""
function plot_data!(
    workspace::Workspace.Workspace,
    ::Type{RegisteredPlot{Kind,Label}},
    items::Vector{<:AbstractDataItem},
    figure::Figure,
)::Nothing where {Kind,Label}
    recipe = workspace.project.plots[Kind][String(Label)]
    Profiling.@profile_span workspace.profiler :visualization :draw_plot Profiling.ProfileAttributes(
        kind=Kind,
        items=Int64(length(items)),
    ) begin
        recipe.draw(workspace, items, figure)
    end
    return nothing
end

end
