module Visualization

using GLMakie: Figure
using InteractiveUtils: subtypes

import ..ItemIndex: ItemRecord
import ..Projects: AbstractDataItem, source_label
import ..Workspace

"""
Base type for visualizers selectable by the browser and Julia API.

`PlotKind` is an internal identity used by the engine to dispatch plotting. Projects do not define
their own subtypes; they call `register_plot!`, and the engine represents each registered plot as a
`RegistryPlot{kind,label}` (below).
"""
abstract type PlotKind end

"""
Internal plot identity for a plot registered with `register_plot!`, parameterized by the item kind
it applies to and the label that distinguishes plots for that kind.
"""
struct RegistryPlot{Kind,Label} <: PlotKind end

"""The item-kind symbol a `RegistryPlot` draws."""
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
    records::Vector{ItemRecord},
)::Figure
    error(
        "No setup_plot implementation for $(source_label(workspace.source)) " *
        "plot kind '$plot_kind'",
    )
end

"""Draw selected items into an existing figure."""
function plot_data!(
    workspace::Workspace.Workspace,
    plot_kind::Type{<:PlotKind},
    records::Vector{ItemRecord},
    figure::Figure,
)::Nothing
    error(
        "No plot_data! implementation for $(source_label(workspace.source)) " *
        "plot kind '$plot_kind'",
    )
end

# --- Registered-plot bridge ------------------------------------------------
# A RegistryPlot{kind,label} dispatches the engine's plot calls to the project's registered setup/draw
# callbacks. `project.plots` is the recipe table filled by register_plot!.

"""
Materialize the loaded, data-bearing items for a selection of records (engine-owned).

Each item carries its payload as `item.data`, so callbacks read `item.data` instead of a parallel
array. Recipe-API items resolve through the processed-data cache; type-API items rebuild from source.
"""
_data_items(workspace::Workspace.Workspace, records::Vector{ItemRecord}) =
    Workspace.materialize_items(workspace, records)

"""
Build the figure for a registered plot by running its `setup` callback.

The package materializes the loaded `items` (each with `item.data`) before invoking the recipe, so
`setup` sizes the figure layout to the data without resolving it itself.
"""
function setup_plot(
    workspace::Workspace.Workspace,
    ::Type{RegistryPlot{Kind,Label}},
    records::Vector{ItemRecord},
)::Figure where {Kind,Label}
    recipe = workspace.source.project.plots[Kind][String(Label)]
    return recipe.setup(workspace, _data_items(workspace, records))::Figure
end

"""
Draw a registered plot into its figure by running its `draw` callback.

The package materializes the loaded `items` (each with `item.data`) before invoking the recipe, so
`draw` reads `item.data` directly and never resolves data itself.
"""
function plot_data!(
    workspace::Workspace.Workspace,
    ::Type{RegistryPlot{Kind,Label}},
    records::Vector{ItemRecord},
    figure::Figure,
)::Nothing where {Kind,Label}
    recipe = workspace.source.project.plots[Kind][String(Label)]
    recipe.draw(workspace, _data_items(workspace, records), figure)
    return nothing
end

"""Every plot registered for an item kind, sorted by label."""
function registered_plot_kinds(source, kind::Symbol)::Vector{Type{<:PlotKind}}
    recipes = get(source.project.plots, kind, nothing)
    recipes === nothing && return Type{<:PlotKind}[]
    return Type{<:PlotKind}[
        RegistryPlot{kind,Symbol(label)}
        for label in sort!(collect(keys(recipes)))
    ]
end

"""Human label for a registered plot kind, taken from its recipe."""
function plot_kind_label(source, ::Type{RegistryPlot{Kind,Label}})::String where {Kind,Label}
    return source.project.plots[Kind][String(Label)].label
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

"""
Draw raw-data items directly while bypassing processed data.

`items` carry the unprocessed payload as `item.data`, so a debug callback can tune the analysis live.
"""
function debug_plot(
    ::Workspace.Workspace,
    items::Vector{<:AbstractDataItem};
    kwargs...,
)
    error("Debug plots are not implemented for this project")
end

end
