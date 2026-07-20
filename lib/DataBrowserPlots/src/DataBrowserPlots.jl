"""GLMakie plot rendering as a GUI extension."""
module DataBrowserPlots

using GLMakie: Figure

import DataBrowserAPI
import DataBrowserAPI:
    AbstractDataItem,
    Project
import DataBrowserAPI.ItemIndex: ItemRecord
import DataBrowserCore.Workspace
using DataBrowserAPI: @time_dbg
import DataBrowserGUI
const Browser = DataBrowserGUI.Browser

include("plots.jl")

include("MakieIntegration.jl")
using .MakieImguiIntegration

include("PlotState.jl")
include("PlotPanel.jl")
include("TablePlotPanel.jl")
include("PlotsExtension.jl")

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

function setup_plot(
    workspace::Workspace.Workspace,
    ::Type{RegisteredPlot{Kind,Label}},
    items::Vector{<:AbstractDataItem},
)::Figure where {Kind,Label}
    recipe = _plot_recipes(workspace.project)[Kind][String(Label)]
    return @time_dbg "setup_plot" recipe.setup(workspace, items)::Figure
end

function plot_data!(
    workspace::Workspace.Workspace,
    ::Type{RegisteredPlot{Kind,Label}},
    items::Vector{<:AbstractDataItem},
    figure::Figure,
)::Nothing where {Kind,Label}
    recipe = _plot_recipes(workspace.project)[Kind][String(Label)]
    @time_dbg "draw_plot" recipe.draw(workspace, items, figure)
    return nothing
end

function __init__()
    Browser.register_gui_extension!(PlotsExtension)
end

export PlotsExtension,
    plots_extension,
    PlotState,
    PlotViewState,
    PersistedPlotView,
    PlotKind,
    PlotRecipe,
    RegisteredPlot,
    register_plot!,
    setup_plot,
    plot_data!,
    registered_plot_kinds,
    plot_kind_label,
    plot_kind_name,
    plot_kind_from_name,
    plot_kinds

end
