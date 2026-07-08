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

"""Every plot registered for an item kind, sorted by label."""
function registered_plot_kinds(project, kind::Symbol)::Vector{Type{<:PlotKind}}
    recipes = get(project.plots, kind, nothing)
    recipes === nothing && return Type{<:PlotKind}[]
    return Type{<:PlotKind}[
        RegisteredPlot{kind,Symbol(label)}
        for label in sort!(collect(keys(recipes)))
    ]
end

"""Human label for a registered plot kind, taken from its recipe."""
function plot_kind_label(project, ::Type{RegisteredPlot{Kind,Label}})::String where {Kind,Label}
    return project.plots[Kind][String(Label)].label
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
