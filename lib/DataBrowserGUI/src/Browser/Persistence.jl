using TOML

using DataBrowserAPI: project_name
using DataBrowserCore.ItemIndex: ItemRecord
import DataBrowserCore.Workspace
using DataBrowserAPI:
    PlotKind,
    plot_kind_from_name,
    plot_kind_name

"""Convert a user path into the stable absolute form used by persistence."""
function _normalize_project_path(path::AbstractString)::String
    return normpath(abspath(expanduser(String(path))))
end

# ---------------------------------------------------------------------------
# Project-local browser state
# ---------------------------------------------------------------------------

const PROJECT_VIEW_FILENAME = "databrowser.toml"

function _project_view_file_path(root_path::AbstractString)::String
    return joinpath(_normalize_project_path(root_path), PROJECT_VIEW_FILENAME)
end

"""
Decode TOML values according to the fields of the requested browser-state type.

Missing fields use the type defaults; incorrectly typed fields are errors.
"""
function _project_view_from_toml(::Type{String}, value::Any)::String
    value isa AbstractString || error("Expected string, got $(typeof(value))")
    return String(value)
end

function _project_view_from_toml(::Type{Bool}, value::Any)::Bool
    value isa Bool || error("Expected boolean, got $(typeof(value))")
    return value
end

function _project_view_from_toml(::Type{Vector{T}}, value::Any)::Vector{T} where {T}
    value isa AbstractVector || error("Expected array, got $(typeof(value))")
    return [_project_view_from_toml(T, item) for item in value]
end

function _project_view_from_toml(::Type{Dict{String,String}}, value::Any)::Dict{String,String}
    value isa AbstractDict || error("Expected table, got $(typeof(value))")
    return Dict(String(key) => _project_view_from_toml(String, item) for (key, item) in value)
end

function _project_view_from_toml(
    ::Type{Dict{String,Dict{String,Any}}},
    value::Any,
)::Dict{String,Dict{String,Any}}
    value isa AbstractDict || error("Expected table, got $(typeof(value))")
    return Dict(
        String(key) => Dict{String,Any}(
            String(inner_key) => inner_value
            for (inner_key, inner_value) in item
        )
        for (key, item) in value
    )
end

function _project_view_from_toml(::Type{T}, data::Any)::T where {T}
    data isa AbstractDict || error("Expected table for $T, got $(typeof(data))")
    defaults = T()
    return T(; (
        name => begin
            key = String(name)
            haskey(data, key) ?
                _project_view_from_toml(fieldtype(T, name), data[key]) :
                getfield(defaults, name)
        end
        for name in fieldnames(T)
    )...)
end

_project_view_to_toml(value::AbstractString)::String = String(value)
_project_view_to_toml(value::Bool)::Bool = value
_project_view_to_toml(value::AbstractVector)::Vector = [_project_view_to_toml(item) for item in value]
_project_view_to_toml(value::AbstractDict)::Dict{String,Any} =
    Dict{String,Any}(String(key) => _project_view_to_toml(item) for (key, item) in value)

"""Convert persisted browser-state structs into values accepted by `TOML.print`."""
function _project_view_to_toml(value::Any)::Dict{String,Any}
    return Dict{String,Any}(
        String(name) => _project_view_to_toml(getfield(value, name))
        for name in fieldnames(typeof(value))
    )
end

"""Read project-local browser state when it exists."""
function _load_project_view(root_path::AbstractString)::PersistedProjectView
    path = _project_view_file_path(root_path)
    isfile(path) || return PersistedProjectView()
    return _project_view_from_toml(
        PersistedProjectView,
        TOML.parsefile(path),
    )
end

"""Write project-local browser state beside the source data."""
function _save_project_view(root_path::AbstractString, view::PersistedProjectView)::Nothing
    open(_project_view_file_path(root_path), "w") do io
        TOML.print(io, _project_view_to_toml(view))
    end
    return nothing
end

"""Resolve stable item ids against the current workspace index."""
function _items_for_ids(
    state::BrowserState,
    ids::Vector{String},
)::Vector{ItemRecord}
    workspace = state.workspace
    workspace isa Workspace.Workspace || return ItemRecord[]
    index = workspace.index.items
    return [index[id] for id in ids if haskey(index, id)]
end

"""Convert one runtime plot window into its persisted form."""
function _persisted_plot_view(view::PlotViewState)::PersistedPlotView
    return PersistedPlotView(
        id=view.id,
        title=view.title,
        plot_kind=view.plot_kind === nothing ? "" : plot_kind_name(view.plot_kind),
        live=view.live,
        items=copy(view.item_ids),
    )
end

"""Capture the browser controls that belong in the source-root state file."""
function _project_view_from_browser(
    state::BrowserState,
)::Union{Nothing,PersistedProjectView}
    workspace = state.workspace
    workspace isa Workspace.Workspace || return nothing
    plots = state.plots

    return PersistedProjectView(
        project=project_name(workspace.project),
        tree=PersistedTreeView(
            expanded=copy(state.expanded_collection_paths),
            selected=copy(workspace.selection.collection_paths),
            filter=state.tree_filter,
        ),
        items=PersistedItemsView(
            selected=copy(workspace.selection.item_ids),
            filter=state.item_filter,
        ),
        plot_kinds=Dict(
            String(kind) => plot_kind_name(plot_kind)
            for (kind, plot_kind) in plots.kind_by_item
        ),
        main_plot=_persisted_plot_view(plots.main),
        plot_windows=[_persisted_plot_view(view) for view in plots.windows],
        extensions=_persisted_extensions(state),
    )
end

"""Merge extension save views with any persisted ids not loaded this session."""
function _persisted_extensions(state::BrowserState)::Dict{String,Dict{String,Any}}
    extensions = Dict{String,Dict{String,Any}}()
    for (id, view) in state.saved_project_view.extensions
        extensions[id] = copy(view)
    end
    for ext in state.extensions
        extensions[extension_id(ext)] = save_view(ext, state)
    end
    return extensions
end

"""Apply loaded project-local state to browser controls and workspace selection."""
function _apply_project_view!(
    state::BrowserState,
    view::PersistedProjectView,
)::Nothing
    workspace = state.workspace::Workspace.Workspace
    plots = state.plots
    state.expanded_collection_paths = copy(view.tree.expanded)
    workspace.selection.collection_paths = copy(view.tree.selected)
    workspace.selection.item_ids = copy(view.items.selected)
    state.tree_filter = view.tree.filter
    state.item_filter = view.items.filter
    empty!(plots.kind_by_item)
    for (kind, plot_kind_name) in view.plot_kinds
        isempty(plot_kind_name) && continue
        plot_kind = plot_kind_from_name(plot_kind_name)
        plot_kind === nothing && continue
        plots.kind_by_item[Symbol(kind)] = plot_kind
    end
    plots.main = PlotViewState(
        id="main",
        title="Plot Area",
        live=view.main_plot.live,
        item_ids=copy(view.main_plot.items),
        plot_kind=isempty(view.main_plot.plot_kind) ?
            nothing :
            plot_kind_from_name(view.main_plot.plot_kind),
    )
    plots.windows = [
        PlotViewState(
            id=plot_view.id,
            title=isempty(plot_view.title) ? "Plot" : plot_view.title,
            live=plot_view.live,
            item_ids=copy(plot_view.items),
            plot_kind=isempty(plot_view.plot_kind) ?
                nothing :
                plot_kind_from_name(plot_view.plot_kind),
        )
        for plot_view in view.plot_windows
    ]
    counters = [
        parse(Int, only(match_result.captures))
        for plot_view in view.plot_windows
        for match_result in (match(r"^plot_(\d+)$", plot_view.id),)
        if match_result !== nothing
    ]
    plots.next_window_id = isempty(counters) ? 0 : maximum(counters)

    for ext in state.extensions
        id = extension_id(ext)
        haskey(view.extensions, id) && load_view!(ext, state, view.extensions[id])
    end

    _reset_project_filter_widgets!(state)
    return nothing
end

"""Write project-local state only when its serialized value changed."""
function _save_project_view_if_changed!(state::BrowserState)::Nothing
    view = _project_view_from_browser(state)
    view === nothing && return nothing
    _project_view_to_toml(state.saved_project_view) ==
        _project_view_to_toml(view) && return nothing
    workspace = state.workspace::Workspace.Workspace
    hasproperty(workspace.source, :root_path) || return nothing
    _save_project_view(workspace.source.root_path, view)
    state.saved_project_view = view
    return nothing
end
