"""GLMakie plot rendering as a GUI extension."""
module DataBrowserPlots

using GLMakie: Figure

import DataBrowserAPI
import DataBrowserAPI:
    AbstractDataItem,
    PlotKind,
    RegisteredPlot,
    plot_data!,
    setup_plot
import DataBrowserAPI.ItemIndex: ItemRecord
import DataBrowserCore.Workspace
import DataBrowserProfiling as Profiling
import DataBrowserGUI
const Browser = DataBrowserGUI.Browser

include("MakieIntegration.jl")
using .MakieImguiIntegration

include("PlotState.jl")
include("PlotPanel.jl")
include("TablePlotPanel.jl")
include("PlotsExtension.jl")

"""Materialize records and create the figure layout required by a visualizer."""
function DataBrowserAPI.setup_plot(
    workspace::Workspace.Workspace,
    plot_kind::Type{<:PlotKind},
    records::Vector{ItemRecord},
)::Figure
    return setup_plot(
        workspace, plot_kind, Workspace.materialize_items(workspace, records))
end

"""Materialize records and draw their items into an existing figure."""
function DataBrowserAPI.plot_data!(
    workspace::Workspace.Workspace,
    plot_kind::Type{<:PlotKind},
    records::Vector{ItemRecord},
    figure::Figure,
)::Nothing
    plot_data!(
        workspace, plot_kind, Workspace.materialize_items(workspace, records), figure)
    return nothing
end

function DataBrowserAPI.setup_plot(
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

function DataBrowserAPI.plot_data!(
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

function __init__()
    Browser.register_gui_extension!(PlotsExtension)
end

export PlotsExtension, plots_extension, request_plot_profile!
export PlotState, PlotViewState, PersistedPlotView

end
