using InteractiveUtils: subtypes

using DataBrowserAPI: AbstractDataItem, Project

"""
One registered plot recipe for an item kind.

`setup(workspace, items)` returns the figure handle; `draw(workspace, items, figure)` fills it.
`items` are the loaded, data-bearing items for the selection (`Vector{<:AbstractDataItem}`);
`process(item)`, if registered, already returned the item shape these callbacks receive. Neither
callback resolves data itself.
"""
struct PlotRecipe
    kind::Symbol
    label::String
    setup::Function
    draw::Function
end

const PROJECT_PLOT_RECIPES = WeakKeyDict{Project, Dict{Symbol, Dict{String, PlotRecipe}}}()

"""Plot recipes registered for one project."""
function _plot_recipes(project::Project)::Dict{Symbol, Dict{String, PlotRecipe}}
    return get!(PROJECT_PLOT_RECIPES, project) do
        Dict{Symbol, Dict{String, PlotRecipe}}()
    end
end

"""
Base type for visualizers selectable by the browser and Julia API.

`PlotKind` is an internal identity used by the engine to dispatch plotting. Projects do not define
their own subtypes; they call `register_plot!`, and the engine represents each registered plot as a
`RegisteredPlot{kind,label}`.
"""
abstract type PlotKind end

"""
Internal plot identity for a plot registered with `register_plot!`, parameterized by the item kind
it applies to and the label that distinguishes plots for that kind.
"""
struct RegisteredPlot{Kind,Label} <: PlotKind end

"""The item-kind symbol a `RegisteredPlot` draws."""
plot_kind_symbol(::Type{RegisteredPlot{Kind,Label}}) where {Kind,Label} = Kind

"""Stable name used to persist a plot choice in workspace view state."""
plot_kind_name(plot_kind::Type{<:PlotKind})::String = String(nameof(plot_kind))
plot_kind_name(::Type{RegisteredPlot{Kind,Label}}) where {Kind,Label} =
    string(Kind, "::", Label)

"""Resolve a persisted plot name back to its internal identity (a registered plot kind)."""
function plot_kind_from_name(name::AbstractString)::Union{Nothing,Type{<:PlotKind}}
    parts = split(name, "::"; limit=2)
    length(parts) == 2 || return nothing
    kind, label = parts
    return RegisteredPlot{Symbol(kind),Symbol(label)}
end

"""Materialize records and create the figure layout required by a visualizer."""
function setup_plot end

"""Materialize records and draw their items into an existing figure."""
function plot_data! end

"""
Register (or replace) one plot recipe for an item kind.

`setup(workspace, items)` builds and returns the figure; `draw(workspace, items, figure)` fills it
in. `items` are the loaded, data-bearing items for the selection (`Vector{<:AbstractDataItem}`);
each has already passed through `process(item)`, if registered, so neither callback resolves item
data itself. A kind may have multiple plots; re-registering the same `label` for the same kind
replaces that plot, which keeps REPL iteration stable.
"""
function register_plot!(
    project::Project,
    kind::Symbol;
    label::AbstractString,
    setup::Function,
    draw::Function,
)::Project
    label_string = String(label)
    recipes = get!(_plot_recipes(project), kind) do
        Dict{String,PlotRecipe}()
    end
    recipes[label_string] = PlotRecipe(kind, label_string, setup, draw)
    return project
end

"""Every plot registered for an item kind, sorted by label."""
function registered_plot_kinds(project, kind::Symbol)::Vector{Type{<:PlotKind}}
    recipes = get(_plot_recipes(project), kind, nothing)
    recipes === nothing && return Type{<:PlotKind}[]
    return Type{<:PlotKind}[
        RegisteredPlot{kind,Symbol(label)}
        for label in sort!(collect(keys(recipes)))
    ]
end

"""Human label for a registered plot kind, taken from its recipe."""
function plot_kind_label(project, ::Type{RegisteredPlot{Kind,Label}})::String where {Kind,Label}
    return _plot_recipes(project)[Kind][String(Label)].label
end

"""
Return every loaded concrete visualizer type.

Visualizers are discovered from Julia's type hierarchy, so external project packages require no
registration with the browser engine.
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
