const PLOT_HELP_TEXT = "Live follows the browser selection.\nDetach opens an independent plot window.\nExport saves the current figure.\nScroll zooms, right-drag pans, Ctrl-click resets limits."

"""
Find a plot type previously written to project state.

An empty name means no plot was selected. Any other unknown name is an invalid project state.
"""
function _plot_kind_from_name(name::AbstractString)::Union{Nothing,Type{<:PlotKind}}
    isempty(name) && return nothing
    for plot_kind in plot_kinds()
        String(nameof(plot_kind)) == name && return plot_kind
    end
    error("Unknown plot kind '$name' in measurementbrowser.toml")
end

"""Remember the plot chosen for each kind of selected measurement."""
function _remember_plot_kind!(
    plots::PlotState,
    measurements::Vector{MeasurementInfo},
    plot_kind::Type{<:PlotKind},
)::Nothing
    for measurement in measurements
        plots.kind_by_measurement[measurement.measurement_kind] = plot_kind
    end
    return nothing
end

"""Reserve the next unique id for a detached plot window."""
function _next_plot_window_id!(plots::PlotState)::String
    plots.next_window_id += 1
    return "plot_$(plots.next_window_id)"
end

"""Create a detached plot for one measurement."""
function _measurement_plot_window(
    plots::PlotState,
    measurement::MeasurementInfo,
)::PlotViewState
    return PlotViewState(
        id=_next_plot_window_id!(plots),
        title=measurement.clean_title,
        live=false,
        measurement_ids=[measurement.unique_id],
        measurements=[measurement],
        plot_kind=get(plots.kind_by_measurement, measurement.measurement_kind, nothing),
    )
end

"""Discard rendered figures so every open plot is redrawn."""
function _clear_plot_views!(plots::PlotState)::Nothing
    for view in (plots.main, plots.windows...)
        view.figure = nothing
        view.last_key = nothing
        view.error = ""
        view.export_error = ""
    end
    return nothing
end

"""
Draw one plot and keep a short UI error when drawing fails.

The console receives the complete exception together with the project, plot type, and source files.
"""
function _draw_plot_view!(
    state::BrowserState,
    view::PlotViewState,
    plot_key::Tuple,
    measurements::Vector{MeasurementInfo},
    plot_kind::Type{<:PlotKind},
)::Nothing
    workspace = state.workspace::Workspace.Workspace
    project = workspace.project
    plots = state.plots
    started_ns = time_ns()
    draw_alloc = 0
    try
        figure = nothing
        draw_alloc = @allocated begin
            figure = if plots.debug
                source_data = DataFrame[]
                sizehint!(source_data, length(measurements))
                for measurement in measurements
                    source_file = index_source_file(measurement.filepath)
                    push!(source_data, load_source_data(project, source_file; measurement))
                end
                debug_plot(workspace, measurements, source_data; plot_kind)
            else
                result = setup_plot(workspace, plot_kind, measurements)
                plot_data!(workspace, plot_kind, measurements, result)
                result
            end
        end
        figure === nothing && error("Plot renderer returned no figure.")
        _append_perf_sample!(
            state,
            :plot_draw,
            (time_ns() - started_ns) / 1e6,
            draw_alloc,
        )
        view.figure = figure
        view.last_key = plot_key
        view.debug = plots.debug
        view.error = ""
    catch err
        bt = catch_backtrace()
        _append_perf_sample!(
            state,
            :plot_draw,
            (time_ns() - started_ns) / 1e6,
            draw_alloc,
        )
        view.figure = nothing
        view.last_key = plot_key
        view.debug = plots.debug
        summary = first(split(sprint(showerror, err), '\n'; limit=2))
        view.error = "Plot failed: $summary. See the console for full details."
        measurement_context = join(
            [
                "$(measurement.clean_title) ($(measurement.measurement_kind))\n" *
                "  $(measurement.filepath)"
                for measurement in measurements
            ],
            "\n",
        )
        @error(
            "Plot drawing failed\n" *
            "Project: $(project_name(project))\n" *
            "Plot: $(nameof(plot_kind))\n" *
            "Measurements: $(length(measurements))\n$measurement_context",
            exception=(err, bt),
        )
    end
    return nothing
end

"""Render the controls shared by the main and detached plot windows."""
function _render_plot_toolbar!(
    state::BrowserState,
    view::PlotViewState,
    measurements::Vector{MeasurementInfo},
    selected_measurements::Vector{MeasurementInfo},
)::Nothing
    plots = state.plots
    current = view.plot_kind
    if ig.BeginCombo(
        "##plot_kind_$(view.id)",
        current === nothing ? "Choose plot kind" : String(nameof(current)),
    )
        for candidate in plot_kinds()
            if ig.Selectable(String(nameof(candidate)), current === candidate)
                view.plot_kind = candidate
                view.last_key = nothing
                view.error = ""
                _remember_plot_kind!(plots, measurements, candidate)
            end
        end
        ig.EndCombo()
    end

    live = view.live
    ig.SameLine()
    if @c ig.Checkbox("Live##live_$(view.id)", &live)
        view.live = live
        view.last_key = nothing
        if !live
            view.measurements = copy(selected_measurements)
            view.measurement_ids =
                [measurement.unique_id for measurement in selected_measurements]
        end
    end

    if view === plots.main
        can_detach = current !== nothing && !isempty(measurements)
        ig.SameLine()
        !can_detach && ig.BeginDisabled()
        if ig.Button("Detach") && can_detach
            push!(
                plots.windows,
                PlotViewState(
                    id=_next_plot_window_id!(plots),
                    title=length(measurements) == 1 ?
                          only(measurements).clean_title :
                          "$(length(measurements)) measurements",
                    live=false,
                    measurement_ids=[
                        measurement.unique_id for measurement in measurements
                    ],
                    measurements=copy(measurements),
                    plot_kind=current,
                ),
            )
        end
        !can_detach && ig.EndDisabled()
    end

    can_export = view.figure !== nothing
    ig.SameLine()
    !can_export && ig.BeginDisabled()
    if ig.Button("Export##export_$(view.id)") && can_export
        name = current === nothing ? "plot" : lowercase(String(nameof(current)))
        default_name = length(measurements) == 1 ?
            "$(splitext(basename(only(measurements).filepath))[1])-$name.png" :
            "$(length(measurements))-measurements-$name.png"
        path = save_file(default_name; filterlist="png,jpg,jpeg,svg,pdf")
        if !isempty(path)
            isempty(splitext(path)[2]) && (path *= ".png")
            try
                Makie.save(path, view.figure)
                view.export_error = ""
            catch err
                bt = catch_backtrace()
                summary = first(split(sprint(showerror, err), '\n'; limit=2))
                view.export_error =
                    "Export failed: $summary. See the console for full details."
                @error("Plot export failed\nOutput: $path", exception=(err, bt))
            end
        end
    end
    !can_export && ig.EndDisabled()

    ig.SameLine()
    if ig.Button("?##plot_help_$(view.id)")
        ig.OpenPopup("plot_help_popup_$(view.id)")
    end
    if ig.BeginItemTooltip()
        ig.TextUnformatted("Help")
        ig.EndTooltip()
    end
    if ig.BeginPopup("plot_help_popup_$(view.id)")
        ig.PushTextWrapPos(ig.GetFontSize() * 38.0)
        ig.TextUnformatted(PLOT_HELP_TEXT)
        ig.PopTextWrapPos()
        ig.EndPopup()
    end
    isempty(view.export_error) || ig.TextWrapped(view.export_error)
    return nothing
end

"""Render a figure or the reason no figure is available."""
function _render_plot_body!(
    state::BrowserState,
    view::PlotViewState,
    status::Symbol,
)::Nothing
    plots = state.plots
    if plots.debug
        ig.TextColored((0.2, 0.8, 0.2, 1.0), "Debug Plot Mode")
        ig.SameLine()
        _helpmarker("Debug mode bypasses cache and reads source files directly.")
    end

    if view.figure !== nothing
        makie_id = view === plots.main ?
            "measurement_plot" :
            "measurement_plot_$(replace(view.id, '/' => '_'))"
        _time!(state, :makie_fig) do
            MakieFigure(
                makie_id,
                view.figure;
                auto_resize_x=true,
                auto_resize_y=true,
            )
        end
    elseif status == :needs_kind
        ig.TextDisabled("No plot kind selected")
    elseif status == :empty
        ig.TextDisabled("No measurement selected")
    else
        ig.TextColored((1.0, 0.4, 0.4, 1.0), "Plot generation failed")
        isempty(view.error) || ig.TextWrapped(view.error)
    end
    return nothing
end

"""
Update and render one plot window.

Live views copy the current browser selection. Static views retain their own measurement ids and
references.
"""
function _render_plot_view!(
    state::BrowserState,
    view::PlotViewState,
    selected_measurements::Vector{MeasurementInfo},
)::Nothing
    plots = state.plots

    if view.live
        view.measurements = copy(selected_measurements)
        view.measurement_ids = [
            measurement.unique_id for measurement in selected_measurements
        ]
    end
    measurements = view.measurements

    if view === plots.main && view.live
        measurement_kind =
            !isempty(measurements) &&
            all(
                measurement ->
                    measurement.measurement_kind == first(measurements).measurement_kind,
                measurements,
            ) ?
            first(measurements).measurement_kind :
            nothing
        if view.measurement_kind !== measurement_kind
            view.measurement_kind = measurement_kind
            if measurement_kind !== nothing
                view.plot_kind =
                    get(plots.kind_by_measurement, measurement_kind, nothing)
            end
            view.last_key = nothing
        end
    end

    status = isempty(measurements) ?
        :empty :
        view.plot_kind === nothing ? :needs_kind : :ready
    if status == :ready
        workspace = state.workspace::Workspace.Workspace
        plot_key = (
            project_name(workspace.project),
            view.id,
            nameof(view.plot_kind),
            [measurement.unique_id for measurement in measurements],
            plots.debug,
        )
        view.last_key == plot_key ||
            _draw_plot_view!(
                state,
                view,
                plot_key,
                measurements,
                view.plot_kind,
            )
    elseif view.live
        view.figure = nothing
        view.last_key = nothing
    end

    _render_plot_toolbar!(state, view, measurements, selected_measurements)
    ig.Separator()
    _render_plot_body!(state, view, status)
    return nothing
end

"""Warm GLMakie once in a hidden window before the first visible plot."""
function _ensure_plot_runtime_warmed!(state::BrowserState)::Nothing
    plots = state.plots
    plots.runtime_warmed && return nothing
    flags = ig.ImGuiWindowFlags_NoDecoration | ig.ImGuiWindowFlags_NoInputs |
            ig.ImGuiWindowFlags_NoBackground | ig.ImGuiWindowFlags_NoSavedSettings |
            ig.ImGuiWindowFlags_NoNav | ig.ImGuiWindowFlags_NoFocusOnAppearing
    ig.SetNextWindowPos((-10_000.0, -10_000.0), ig.ImGuiCond_Always)
    ig.SetNextWindowSize((8.0, 8.0), ig.ImGuiCond_Always)
    if ig.Begin("###plot_runtime_warmup", C_NULL, flags)
        if plots.warmup_figure === nothing
            figure = Figure(size=(64, 48))
            axis = Axis(figure[1, 1], xlabel="x", ylabel="y", title="warm")
            lines!(axis, [0.0, 1.0], [0.0, 1.0], color=:blue)
            plots.warmup_figure = figure
        end
        MakieFigure(
            "_plot_runtime_warmup",
            plots.warmup_figure;
            auto_resize_x=false,
            auto_resize_y=false,
        )
        plots.runtime_warmed = true
    end
    ig.End()
    return nothing
end

"""Render the main plot window."""
function render_plot_window(state::BrowserState)::Nothing
    _, selected, _ = _project_visible_selection(state)
    main = state.plots.main
    if ig.Begin(main.title)
        _render_plot_view!(state, main, selected)
    end
    ig.End()
    return nothing
end

"""Render every detached plot window and remove windows closed by the user."""
function render_additional_plot_windows(state::BrowserState)::Nothing
    plots = state.plots
    _, selected, _ = _project_visible_selection(state)
    kept = PlotViewState[]

    for view in plots.windows
        open = Ref(true)
        window_id = replace(view.id, '/' => '_')
        if ig.Begin("Plot: $(view.title)###plot_window_$window_id", open)
            view.debug == plots.debug || (view.last_key = nothing)
            _render_plot_view!(state, view, selected)
        end
        ig.End()
        open[] && push!(kept, view)
    end
    plots.windows = kept
    return nothing
end
