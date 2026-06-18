using TOML

using ..Projects:
    PROJECTS,
    project_name,
    source_label
using ..ItemIndex: ItemRecord, item_record_key
import ..Workspace
using ..Visualization:
    PlotKind,
    plot_kind_from_name,
    plot_kind_name

# ---------------------------------------------------------------------------
# App preferences
# ---------------------------------------------------------------------------

const _MAX_RECENT_PROJECTS = 12

function _prefs_path()::String
    return joinpath(homedir(), ".config", "MeasurementBrowser", "prefs.toml")
end

"""Read app-wide preferences, returning an empty table before the first save."""
function _load_prefs()::Dict{String,Any}
    path = _prefs_path()
    isfile(path) || return Dict{String,Any}()
    return TOML.parsefile(path)
end

"""Write app-wide preferences outside any source project."""
function _save_prefs(data::Dict{String,Any})::Nothing
    path = _prefs_path()
    mkpath(dirname(path))
    open(path, "w") do io
        TOML.print(io, data)
    end
    return nothing
end

"""Convert a user path into the stable absolute form used by persistence."""
function _normalize_project_path(path::AbstractString)::String
    return normpath(abspath(expanduser(String(path))))
end

"""Accept only the automatic choice or a currently registered project name."""
function _sanitize_project_preference(pref::AbstractString)::String
    pref == "auto" && return "auto"
    for project in PROJECTS
        project_name(project) == pref && return String(pref)
    end
    return "auto"
end

"""Decode recent-project entries from app preferences."""
function _parse_recent_projects(
    prefs::Dict{String,Any},
)::Vector{RecentProject}
    recents = RecentProject[]
    raw = get(prefs, "recent_projects", Any[])
    raw isa Vector || return recents

    for entry in raw
        entry isa AbstractDict || continue
        path = get(entry, "path", nothing)
        path isa AbstractString || continue
        path = strip(path)
        isempty(path) && continue
        pref = get(entry, "project_preference", "auto")
        pref = pref isa AbstractString ? pref : "auto"
        push!(
            recents,
            RecentProject(
                path=_normalize_project_path(path),
                project_preference=_sanitize_project_preference(pref),
            ),
        )
    end

    return recents
end

"""Move an opened project to the front of the bounded recent-project list."""
function _update_recent_projects!(
    recents::Vector{RecentProject},
    path::AbstractString,
    pref::AbstractString,
)::Nothing
    norm_path = _normalize_project_path(path)
    filter!(entry -> entry.path != norm_path, recents)
    pushfirst!(
        recents,
        RecentProject(
            path=norm_path,
            project_preference=String(pref),
        ),
    )
    length(recents) > _MAX_RECENT_PROJECTS && resize!(recents, _MAX_RECENT_PROJECTS)
    return nothing
end

"""Save the current project choice and optional source-root preferences."""
function _persist_preferences!(
    state::BrowserState;
    path::Union{Nothing,String}=nothing,
)::Nothing
    prefs = _load_prefs()
    pref = _sanitize_project_preference(state.project_preference)
    state.project_preference = pref
    prefs["project"] = pref

    recents = _parse_recent_projects(prefs)
    if path !== nothing && !isempty(path)
        _update_recent_projects!(recents, path, pref)
        prefs["recent_projects"] = [
            Dict{String,String}(
                "path" => recent.path,
                "project_preference" => recent.project_preference,
            )
            for recent in recents
        ]
    end

    _save_prefs(prefs)
    state.recent_projects = recents
    return nothing
end

"""Find the recent-project entry for one normalized source root."""
function _recent_project_entry_for_path(
    state::BrowserState,
    path::String,
)::Union{Nothing,RecentProject}
    for entry in state.recent_projects
        entry.path == path && return entry
    end
    return nothing
end

"""Return the project choice saved for a source root."""
function _project_preference_for_path(state::BrowserState, path::String)::String
    entry = _recent_project_entry_for_path(state, path)
    if entry !== nothing
        return _sanitize_project_preference(entry.project_preference)
    end
    return _sanitize_project_preference(state.project_preference)
end

"""Persist preferences associated with the currently open source root."""
function _persist_current_project_preferences!(state::BrowserState)::Nothing
    workspace = state.workspace
    workspace isa Workspace.Workspace || return nothing
    hasproperty(workspace.source, :root_path) || return nothing
    _persist_preferences!(state; path=workspace.source.root_path)
    return nothing
end

# ---------------------------------------------------------------------------
# Project-local browser state
# ---------------------------------------------------------------------------

const PROJECT_VIEW_FILENAME = "measurementbrowser.toml"

function _project_view_file_path(root_path::AbstractString)::String
    return joinpath(_normalize_project_path(root_path), PROJECT_VIEW_FILENAME)
end

"""
Decode TOML values according to the fields of the requested browser-state type.

Missing or incorrectly typed fields are errors because no compatibility format is supported.
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

function _project_view_from_toml(::Type{T}, data::Any)::T where {T}
    data isa AbstractDict || error("Expected table for $T, got $(typeof(data))")
    return T(; (
        name => _project_view_from_toml(fieldtype(T, name), data[String(name)])
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
    return _project_view_from_toml(PersistedProjectView, TOML.parsefile(path))
end

"""Write project-local browser state beside the source data."""
function _save_project_view(root_path::AbstractString, view::PersistedProjectView)::Nothing
    open(_project_view_file_path(root_path), "w") do io
        TOML.print(io, _project_view_to_toml(view))
    end
    return nothing
end

"""Resolve stable item keys against the current workspace index."""
function _items_for_keys(
    state::BrowserState,
    keys::Vector{String},
)::Vector{ItemRecord}
    workspace = state.workspace
    workspace isa Workspace.Workspace || return ItemRecord[]
    index = workspace.index.items
    return [index[key] for key in keys if haskey(index, key)]
end

"""Convert one runtime plot window into its persisted form."""
function _persisted_plot_view(view::PlotViewState)::PersistedPlotView
    return PersistedPlotView(
        id=view.id,
        title=view.title,
        plot_kind=view.plot_kind === nothing ? "" : plot_kind_name(view.plot_kind),
        live=view.live,
        items=copy(view.item_keys),
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
        project=source_label(workspace.source),
        tree=PersistedTreeView(
            expanded=copy(state.expanded_collection_paths),
            selected=copy(workspace.selection.collection_paths),
            filter=state.tree_filter,
        ),
        items=PersistedItemsView(
            selected=copy(workspace.selection.item_keys),
            filter=state.item_filter,
        ),
        plot_kinds=Dict(
            String(kind) => plot_kind_name(plot_kind)
            for (kind, plot_kind) in plots.kind_by_item
        ),
        main_plot=_persisted_plot_view(plots.main),
        plot_windows=[_persisted_plot_view(view) for view in plots.windows],
    )
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
    workspace.selection.item_keys = copy(view.items.selected)
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
        item_keys=copy(view.main_plot.items),
        plot_kind=isempty(view.main_plot.plot_kind) ?
            nothing :
            plot_kind_from_name(view.main_plot.plot_kind),
    )
    plots.windows = [
        PlotViewState(
            id=plot_view.id,
            title=isempty(plot_view.title) ? "Plot" : plot_view.title,
            live=plot_view.live,
            item_keys=copy(plot_view.items),
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
