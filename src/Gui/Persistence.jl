# ---------------------------------------------------------------------------
# App preferences
# ---------------------------------------------------------------------------

const _MAX_RECENT_PROJECTS = 12

function _prefs_path()::String
    return joinpath(homedir(), ".config", "MeasurementBrowser", "prefs.toml")
end

function _load_prefs()::Dict{String,Any}
    path = _prefs_path()
    isfile(path) || return Dict{String,Any}()
    return TOML.parsefile(path)
end

function _save_prefs(data::Dict)::Nothing
    path = _prefs_path()
    mkpath(dirname(path))
    open(path, "w") do io
        TOML.print(io, data)
    end
    return nothing
end

function _normalize_project_path(path::AbstractString)::String
    return normpath(abspath(expanduser(String(path))))
end

function _sanitize_project_preference(pref::AbstractString)::String
    pref == "auto" && return "auto"
    for project in PROJECTS
        project_name(project) == pref && return String(pref)
    end
    return "auto"
end

function _sanitize_figure_script_output_dir(value::Any)::String
    value isa AbstractString || return ""
    return String(strip(String(value)))
end

function _parse_recent_projects(prefs::Dict{String,Any})::Vector{Dict{String,String}}
    recents = Dict{String,String}[]
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
        figure_script_output_dir = _sanitize_figure_script_output_dir(get(entry, "figure_script_output_dir", ""))
        push!(recents, Dict{String,String}(
            "path" => _normalize_project_path(path),
            "project_preference" => _sanitize_project_preference(pref),
            "figure_script_output_dir" => figure_script_output_dir,
        ))
    end

    return recents
end

function _update_recent_projects(
    recents::Vector{Dict{String,String}},
    path::AbstractString,
    pref::AbstractString,
    figure_script_output_dir::AbstractString,
)::Vector{Dict{String,String}}
    norm_path = _normalize_project_path(path)
    filter!(entry -> get(entry, "path", "") != norm_path, recents)
    pushfirst!(recents, Dict{String,String}(
        "path" => norm_path,
        "project_preference" => String(pref),
        "figure_script_output_dir" => _sanitize_figure_script_output_dir(figure_script_output_dir),
    ))
    length(recents) > _MAX_RECENT_PROJECTS && resize!(recents, _MAX_RECENT_PROJECTS)
    return recents
end

function _current_figure_script_output_dir(ui_state)::String
    haskey(ui_state, :figure_script_output_dir_buffer) || return ""
    return _sanitize_figure_script_output_dir(_buffer_string(ui_state[:figure_script_output_dir_buffer]))
end

function _persist_preferences!(
    ui_state;
    path::Union{Nothing,String}=nothing,
)::Nothing
    prefs = _load_prefs()
    pref = _sanitize_project_preference(string(get(ui_state, :project_preference, "auto")))
    ui_state[:project_preference] = pref
    prefs["project"] = pref

    recents = _parse_recent_projects(prefs)
    if path !== nothing && !isempty(path)
        _update_recent_projects(
            recents,
            path,
            pref,
            _current_figure_script_output_dir(ui_state),
        )
        prefs["recent_projects"] = recents
    end

    _save_prefs(prefs)
    ui_state[:recent_projects] = recents
    return nothing
end

function _recent_project_entry_for_path(ui_state, path::String)::Union{Nothing,Dict{String,String}}
    recents = get(ui_state, :recent_projects, Dict{String,String}[])
    for entry in recents
        get(entry, "path", "") == path && return entry
    end
    return nothing
end

function _project_preference_for_path(ui_state, path::String)::String
    entry = _recent_project_entry_for_path(ui_state, path)
    if entry !== nothing
        pref = get(entry, "project_preference", "auto")
        return _sanitize_project_preference(pref)
    end
    pref = string(get(ui_state, :project_preference, "auto"))
    return _sanitize_project_preference(pref)
end

function _figure_script_output_dir_for_path(ui_state, path::String)::String
    entry = _recent_project_entry_for_path(ui_state, path)
    entry === nothing && return ""
    return _sanitize_figure_script_output_dir(get(entry, "figure_script_output_dir", ""))
end

function _persist_current_project_preferences!(ui_state)::Nothing
    workspace = get(ui_state, :workspace, nothing)
    workspace isa Workspace.Workspace || return nothing
    _persist_preferences!(ui_state; path=workspace.root_path)
    return nothing
end

# ---------------------------------------------------------------------------
# Project-local browser state
# ---------------------------------------------------------------------------

const PROJECT_VIEW_FILENAME = "measurementbrowser.toml"

"""
    PersistedTreeView(; expanded, selected, filter) -> PersistedTreeView

Saved tree state from `measurementbrowser.toml`.
"""
Base.@kwdef struct PersistedTreeView
    expanded::Vector{String} = String[]
    selected::Vector{String} = String[]
    filter::String = ""
end

"""
    PersistedMeasurementsView(; selected, filter) -> PersistedMeasurementsView

Saved measurement-list state from `measurementbrowser.toml`.
"""
Base.@kwdef struct PersistedMeasurementsView
    selected::Vector{String} = String[]
    filter::String = ""
end

"""
    PersistedPlotView(; id, title, plot_kind, live, measurements) -> PersistedPlotView

Saved plot-window state from `measurementbrowser.toml`.
"""
Base.@kwdef struct PersistedPlotView
    id::String = ""
    title::String = ""
    plot_kind::String = ""
    live::Bool = false
    measurements::Vector{String} = String[]
end

"""
    PersistedProjectView(; kwargs...) -> PersistedProjectView

Saved browser state for one project. It contains only stable ids and strings, never cache data,
runtime jobs, figures, `MeasurementInfo`, `HierarchyNode`, or ImGui objects.
"""
Base.@kwdef struct PersistedProjectView
    project::String = ""
    tree::PersistedTreeView = PersistedTreeView()
    measurements::PersistedMeasurementsView = PersistedMeasurementsView()
    plot_kinds::Dict{String,String} = Dict{String,String}()
    main_plot::PersistedPlotView = PersistedPlotView(id="main", title="Plot Area", live=true)
    plot_windows::Vector{PersistedPlotView} = PersistedPlotView[]
end

function _project_view_file_path(root_path::AbstractString)::String
    return joinpath(_normalize_project_path(root_path), PROJECT_VIEW_FILENAME)
end

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

function _project_view_to_toml(value)::Dict{String,Any}
    return Dict{String,Any}(
        String(name) => _project_view_to_toml(getfield(value, name))
        for name in fieldnames(typeof(value))
    )
end

function _same_project_view(a::PersistedProjectView, b::PersistedProjectView)::Bool
    return _project_view_to_toml(a) == _project_view_to_toml(b)
end

function _load_project_view(root_path::AbstractString)::PersistedProjectView
    path = _project_view_file_path(root_path)
    isfile(path) || return PersistedProjectView()
    return _project_view_from_toml(PersistedProjectView, TOML.parsefile(path))
end

function _save_project_view(root_path::AbstractString, view::PersistedProjectView)::Nothing
    open(_project_view_file_path(root_path), "w") do io
        TOML.print(io, _project_view_to_toml(view))
    end
    return nothing
end

function _load_project_view!(ui_state, root_path::AbstractString, project::AbstractProject)::Nothing
    view = _load_project_view(root_path)
    if !isempty(view.project) && view.project != project_name(project)
        view = PersistedProjectView(project=project_name(project))
    end
    ui_state[:project_view_loaded] = view
    ui_state[:project_view_saved] = view
    return nothing
end

_measurement_ids(measurements)::Vector{String} =
    [measurement.unique_id for measurement in measurements if measurement isa MeasurementInfo]

function _measurements_for_ids(ui_state, ids::Vector{String})::Vector{MeasurementInfo}
    workspace = get(ui_state, :workspace, nothing)
    workspace isa Workspace.Workspace || return MeasurementInfo[]
    index = workspace.index.measurements
    return [index[id] for id in ids if haskey(index, id)]
end

function _plot_kind_name_for_state(state, key::Symbol)::String
    kind = get(state, key, nothing)
    kind isa Type && kind <: PlotKind || return ""
    return String(nameof(kind))
end

function _plot_measurement_ids(entry)::Vector{String}
    get(entry, :live, false) === true && return _measurement_ids(get(entry, :measurements, MeasurementInfo[]))
    ids = get(entry, :measurement_ids, nothing)
    ids isa Vector{String} && return copy(ids)
    return _measurement_ids(get(entry, :measurements, MeasurementInfo[]))
end

function _refresh_plot_measurement_refs!(ui_state)::Nothing
    ids = get(ui_state, :main_plot_measurement_ids, nothing)
    ids isa Vector{String} && (ui_state[:main_plot_measurements] = _measurements_for_ids(ui_state, ids))

    open_plots = get(ui_state, :open_plot_windows, nothing)
    open_plots isa AbstractVector || return nothing
    for entry in open_plots
        entry isa AbstractDict || continue
        ids = get(entry, :measurement_ids, nothing)
        ids isa Vector{String} && (entry[:measurements] = _measurements_for_ids(ui_state, ids))
    end
    return nothing
end

function _persisted_plot_view(entry)::PersistedPlotView
    return PersistedPlotView(
        id=String(get(entry, :target_id, "")),
        title=String(get(entry, :title, "Plot")),
        plot_kind=_plot_kind_name_for_state(entry, :plot_kind),
        live=get(entry, :live, false) === true,
        measurements=_plot_measurement_ids(entry),
    )
end

function _project_view_from_ui_state(ui_state)::Union{Nothing,PersistedProjectView}
    workspace = get(ui_state, :workspace, nothing)
    workspace isa Workspace.Workspace || return nothing

    open_plots = get(ui_state, :open_plot_windows, Dict{Symbol,Any}[])
    main_plot_state = Dict{Symbol,Any}(
        :live => get(ui_state, :main_plot_live, true) === true,
        :measurement_ids => get(ui_state, :main_plot_measurement_ids, nothing),
        :measurements => get(ui_state, :main_plot_measurements, MeasurementInfo[]),
    )

    return PersistedProjectView(
        project=project_name(workspace.project),
        tree=PersistedTreeView(
            expanded=copy(get(ui_state, :expanded_device_paths, String[])),
            selected=copy(workspace.selection.device_paths),
            filter=String(get(ui_state, :tree_filter, "")),
        ),
        measurements=PersistedMeasurementsView(
            selected=copy(workspace.selection.measurement_ids),
            filter=String(get(ui_state, :measurement_filter, "")),
        ),
        plot_kinds=copy(get(ui_state, :plot_kind_by_measurement_kind, Dict{String,String}())),
        main_plot=PersistedPlotView(
            id="main",
            title="Plot Area",
            plot_kind=_plot_kind_name_for_state(ui_state, :main_plot_kind),
            live=main_plot_state[:live],
            measurements=_plot_measurement_ids(main_plot_state),
        ),
        plot_windows=[_persisted_plot_view(entry) for entry in open_plots if entry isa AbstractDict],
    )
end

function _set_plot_kind_from_name!(state, key::Symbol, name::AbstractString)::Nothing
    kind = _plot_kind_from_name(name)
    kind === nothing ? delete!(state, key) : (state[key] = kind)
    return nothing
end

function _reset_project_filter_widgets!(ui_state)::Nothing
    for key in (:_imgui_text_filter_tree, :_imgui_text_filter_meas)
        filter = get(ui_state, key, nothing)
        filter === nothing || ig.Destroy(filter)
        delete!(ui_state, key)
    end
    return nothing
end

function _restore_plot_windows!(ui_state, views::Vector{PersistedPlotView})::Nothing
    ui_state[:open_plot_windows] = [begin
        entry = Dict{Symbol,Any}(
            :target_id => view.id,
            :title => isempty(view.title) ? "Plot" : view.title,
            :live => view.live,
            :measurement_ids => copy(view.measurements),
            :measurements => _measurements_for_ids(ui_state, view.measurements),
        )
        _set_plot_kind_from_name!(entry, :plot_kind, view.plot_kind)
        entry
    end for view in views]
    counters = [parse(Int, only(m.captures)) for view in views for m in (match(r"^plot_(\d+)$", view.id),) if m !== nothing]
    ui_state[:_plot_window_counter] = max(get(ui_state, :_plot_window_counter, 0), isempty(counters) ? 0 : maximum(counters))
    return nothing
end

function _apply_project_view!(ui_state, view::PersistedProjectView)::Nothing
    workspace = ui_state[:workspace]::Workspace.Workspace
    ui_state[:expanded_device_paths] = copy(view.tree.expanded)
    workspace.selection.device_paths = copy(view.tree.selected)
    workspace.selection.measurement_ids = copy(view.measurements.selected)
    ui_state[:tree_filter] = view.tree.filter
    ui_state[:measurement_filter] = view.measurements.filter
    ui_state[:plot_kind_by_measurement_kind] = copy(view.plot_kinds)
    ui_state[:main_plot_live] = view.main_plot.live
    ui_state[:main_plot_measurement_ids] = copy(view.main_plot.measurements)
    ui_state[:main_plot_measurements] = _measurements_for_ids(ui_state, view.main_plot.measurements)
    _set_plot_kind_from_name!(ui_state, :main_plot_kind, view.main_plot.plot_kind)
    delete!(ui_state, :main_plot_measurement_kind)
    delete!(ui_state, :_last_plot_key)
    _restore_plot_windows!(ui_state, view.plot_windows)
    _reset_project_filter_widgets!(ui_state)
    return nothing
end

function _save_project_view_if_changed!(ui_state)::Nothing
    view = _project_view_from_ui_state(ui_state)
    view === nothing && return nothing
    saved = get(ui_state, :project_view_saved, nothing)
    saved isa PersistedProjectView && _same_project_view(saved, view) && return nothing
    workspace = ui_state[:workspace]::Workspace.Workspace
    _save_project_view(workspace.root_path, view)
    ui_state[:project_view_saved] = view
    return nothing
end

function _imgui_text_filter_text(filter)::String
    filter === nothing && return ""
    bytes = UInt8[]
    for char in unsafe_load(filter).InputBuf
        char == 0 && break
        push!(bytes, UInt8(mod(Int(char), 256)))
    end
    return String(bytes)
end

function _sync_imgui_text_filter!(ui_state, key::Symbol, filter)::Nothing
    ui_state[key] = _imgui_text_filter_text(filter)
    return nothing
end

function _imgui_text_filter_for_state!(ui_state, filter_key::Symbol, text_key::Symbol)
    return get!(ui_state, filter_key) do
        ig.ImGuiTextFilter_ImGuiTextFilter(String(get(ui_state, text_key, "")))
    end
end
