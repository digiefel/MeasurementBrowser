const PLOT_HELP_TEXT = "Live follows the browser selection.\nDetach opens an independent plot window.\nExport saves the current figure.\nScroll zooms, right-drag pans, Ctrl-click resets limits."

"""Return the plot choice stored in `state[key]`, or `nothing` when no choice is set."""
function _plot_kind_from_state(state::Dict{Symbol,Any}, key::Symbol)::Union{Nothing,Type{<:PlotKind}}
    haskey(state, key) || return nothing
    plot_kind = state[key]
    plot_kind isa Type && plot_kind <: PlotKind && return plot_kind
    error("Invalid plot kind stored in $key: $(repr(plot_kind))")
end

"""Find the plot type saved as text, for example `RuO2PUNDPlot`."""
function _plot_kind_from_name(name::AbstractString)::Union{Nothing,Type{<:PlotKind}}
    for plot_kind in plot_kinds()
        String(nameof(plot_kind)) == String(name) && return plot_kind
    end
    return nothing
end

"""Return the current plot request key used to avoid redrawing unchanged views."""
function _plot_key(
    project::AbstractProject,
    view_id::AbstractString,
    plot_kind::Type{<:PlotKind},
    measurements::Vector{MeasurementInfo},
    debug::Bool,
)::Tuple
    return (
        project_name(project),
        String(view_id),
        nameof(plot_kind),
        [measurement.unique_id for measurement in measurements],
        debug,
    )
end

"""Save `measurement kind => plot kind` for the selected measurements."""
function _remember_plot_kind!(ui_state::Dict{Symbol,Any}, measurements::Vector{MeasurementInfo}, plot_kind::Type{<:PlotKind})::Nothing
    prefs = get!(ui_state, :plot_kind_by_measurement_kind) do
        Dict{String,String}()
    end
    for measurement in measurements
        prefs[String(measurement.measurement_kind)] = String(nameof(plot_kind))
    end
    return nothing
end

"""Return a new detached plot window id."""
function _next_plot_window_id!(ui_state::Dict{Symbol,Any})::String
    counter = get(ui_state, :_plot_window_counter, 0) + 1
    ui_state[:_plot_window_counter] = counter
    return "plot_$counter"
end

"""Build the detached-window entry for a measurement double-click."""
function _measurement_plot_window_entry(ui_state::Dict{Symbol,Any}, measurement::MeasurementInfo)::Dict{Symbol,Any}
    prefs = get(ui_state, :plot_kind_by_measurement_kind, Dict{String,String}())
    remembered = prefs isa AbstractDict ? get(prefs, String(measurement.measurement_kind), "") : ""
    return Dict{Symbol,Any}(
        :target_id => _next_plot_window_id!(ui_state),
        :filepath => measurement.filepath,
        :measurement_ids => [measurement.unique_id],
        :measurements => [measurement],
        :title => measurement.clean_title,
        :plot_kind => remembered isa AbstractString ? _plot_kind_from_name(remembered) : nothing,
        :live => false,
        :params => _measurement_parameters(measurement),
    )
end

"""Return the plot error text shown in the GUI and console."""
function _plot_error_text(
    project::AbstractProject,
    plot_kind::Type{<:PlotKind},
    measurements::Vector{MeasurementInfo},
    err,
    bt;
    phase::AbstractString,
)::String
    lines = String[
        "Plot failed while $phase: $(sprint(showerror, err))",
        "Project: $(project_name(project))",
        "Plot: $(nameof(plot_kind))",
        "Measurements: $(length(measurements))",
    ]
    for measurement in first(measurements, min(3, length(measurements)))
        push!(lines, "Measurement: $(measurement.clean_title) ($(measurement.measurement_kind))")
        push!(lines, "File: $(measurement.filepath)")
    end
    length(measurements) > 3 && push!(lines, "More files: $(length(measurements) - 3)")
    if bt !== nothing
        stack_lines = String[]
        for frame in stacktrace(bt)
            file = String(frame.file)
            if !(occursin("GlfwOpenGLBackend", file) || occursin("CImGui", file) || occursin("GLFW", file) || endswith(file, "client.jl"))
                push!(stack_lines, "  $(frame.func) at $(basename(file)):$(frame.line)")
                length(stack_lines) >= 10 && break
            end
        end
        isempty(stack_lines) || push!(lines, "Stack:\n" * join(stack_lines, "\n"))
    end
    return join(lines, "\n")
end

"""Forget generated plot figures so the next render builds them again."""
function _clear_plot_views!(ui_state::Dict{Symbol,Any})::Nothing
    for key in (:plot_figure, :_last_plot_key, :plot_error, :plot_export_error)
        delete!(ui_state, key)
    end
    open_plots = get(ui_state, :open_plot_windows, nothing)
    if open_plots isa AbstractVector
        for entry in open_plots
            entry isa Dict{Symbol,Any} || continue
            for key in (:figure, :_last_plot_key, :plot_error, :plot_export_error)
                delete!(entry, key)
            end
        end
    end
    return nothing
end

"""Build and store one plot directly when the view needs a new figure."""
function _draw_plot_view!(
    ui_state::Dict{Symbol,Any},
    view::Dict{Symbol,Any},
    figure_key::Symbol,
    plot_key::Tuple,
    measurements::Vector{MeasurementInfo},
    plot_kind::Type{<:PlotKind},
)::Nothing
    project = ui_state[:project]
    debug = get(ui_state, :debug_plot_mode, false)
    started_ns = time_ns()
    draw_alloc = 0
    try
        figure = nothing
        draw_alloc = @allocated begin
            figure = if debug
                source_data = DataFrame[]
                sizehint!(source_data, length(measurements))
                for measurement in measurements
                    push!(source_data, load_source_data(project, index_source_file(measurement.filepath); measurement))
                end
                debug_plot(project, measurements, source_data; plot_kind)
            else
                fig = setup_plot(project, plot_kind, measurements)
                plot_data!(project, plot_kind, measurements, fig)
                fig
            end
        end
        figure === nothing && error("Plot renderer returned no figure.")
        _append_perf_sample!(ui_state, :plot_draw, (time_ns() - started_ns) / 1e6, draw_alloc)
        view[figure_key] = figure
        view[:_last_plot_key] = plot_key
        view[:debug] = debug
        delete!(view, :plot_error)
    catch err
        _append_perf_sample!(ui_state, :plot_draw, (time_ns() - started_ns) / 1e6, draw_alloc)
        message = _plot_error_text(project, plot_kind, measurements, err, catch_backtrace(); phase="drawing")
        delete!(view, figure_key)
        view[:_last_plot_key] = plot_key
        view[:plot_error] = message
        view[:debug] = debug
        @error "Plot drawing failed\n$message"
    end
    return nothing
end

"""Render one plot view inside an open ImGui window."""
function _render_plot_view!(
    ui_state::Dict{Symbol,Any},
    view::Dict{Symbol,Any},
    measurements::Vector{MeasurementInfo},
    view_id::AbstractString,
    figure_key::Symbol,
    makie_id::AbstractString;
    live_key::Symbol,
    live_default::Bool,
    plot_kind_key::Symbol,
    selected_measurements::Vector{MeasurementInfo},
)::Nothing
    plot_kind = _plot_kind_from_state(view, plot_kind_key)
    status = isempty(measurements) ? :empty : plot_kind === nothing ? :needs_kind : :ready
    plot_key = status == :ready ?
        _plot_key(ui_state[:project], view_id, plot_kind, measurements, get(ui_state, :debug_plot_mode, false)) :
        nothing

    if status == :ready && get(view, :_last_plot_key, nothing) != plot_key
        _draw_plot_view!(ui_state, view, figure_key, plot_key, measurements, plot_kind)
    elseif status != :ready && get(view, live_key, live_default) === true
        delete!(view, figure_key)
        delete!(view, :_last_plot_key)
    end

    _render_plot_toolbar!(
        ui_state,
        view,
        measurements,
        plot_kind,
        view_id;
        live_key,
        live_default,
        plot_kind_key,
        selected_measurements,
    )
    ig.Separator()
    _render_plot_body!(ui_state, view, status, figure_key, makie_id)
    return nothing
end

"""Show a compact error line with the full error in a tooltip."""
function _render_plot_error(message)::Nothing
    isempty(message) && return nothing
    summary = first(split(String(message), '\n'))
    length(summary) > 140 && (summary = first(summary, 137) * "...")
    ig.TextDisabled(summary)
    if ig.BeginItemTooltip()
        ig.PushTextWrapPos(ig.GetFontSize() * 45.0)
        ig.TextUnformatted(String(message))
        ig.PopTextWrapPos()
        ig.EndTooltip()
    end
    return nothing
end

"""Render shared plot toolbar controls."""
function _render_plot_toolbar!(ui_state::Dict{Symbol,Any}, state::Dict{Symbol,Any}, measurements::Vector{MeasurementInfo}, plot_kind::Union{Nothing,Type{<:PlotKind}}, id::AbstractString; live_key::Symbol, live_default::Bool, plot_kind_key::Symbol, selected_measurements::Vector{MeasurementInfo}=MeasurementInfo[])::Nothing
    current = _plot_kind_from_state(state, plot_kind_key)
    if ig.BeginCombo("##plot_kind_$id", current === nothing ? "Choose plot kind" : String(nameof(current)))
        for candidate in plot_kinds()
            if ig.Selectable(String(nameof(candidate)), current === candidate)
                state[plot_kind_key] = candidate
                delete!(state, :_last_plot_key)
                delete!(state, :plot_error)
                _remember_plot_kind!(ui_state, measurements, candidate)
            end
        end
        ig.EndCombo()
    end

    live = get(state, live_key, live_default) === true
    ig.SameLine()
    if @c ig.Checkbox("Live##live_$id", &live)
        state[live_key] = live
        delete!(state, :_last_plot_key)
        if !live
            frozen = copy(isempty(selected_measurements) ? measurements : selected_measurements)
            state[id == "main" ? :main_plot_measurements : :measurements] = frozen
            state[:measurement_ids] = [measurement.unique_id for measurement in frozen]
        end
    end

    if id == "main"
        can_detach = plot_kind !== nothing && !isempty(measurements)
        ig.SameLine()
        !can_detach && ig.BeginDisabled()
        if ig.Button("Detach") && can_detach
            open_plots = get!(ui_state, :open_plot_windows) do
                Vector{Dict{Symbol,Any}}()
            end
            push!(open_plots, Dict{Symbol,Any}(
                :title => length(measurements) == 1 ? only(measurements).clean_title : "$(length(measurements)) measurements",
                :target_id => _next_plot_window_id!(ui_state),
                :measurement_ids => [measurement.unique_id for measurement in measurements],
                :measurements => copy(measurements),
                :plot_kind => plot_kind,
                :live => false,
            ))
        end
        !can_detach && ig.EndDisabled()
    end

    figure = get(state, id == "main" ? :plot_figure : :figure, nothing)
    can_export = figure !== nothing
    ig.SameLine()
    !can_export && ig.BeginDisabled()
    if ig.Button("Export##export_$id") && can_export
        name = plot_kind === nothing ? "plot" : lowercase(String(nameof(plot_kind)))
        default_name = length(measurements) == 1 ?
            "$(splitext(basename(only(measurements).filepath))[1])-$name.png" :
            "$(length(measurements))-measurements-$name.png"
        path = save_file(default_name; filterlist="png,jpg,jpeg,svg,pdf")
        if !isempty(path)
            isempty(splitext(path)[2]) && (path *= ".png")
            try
                Makie.save(path, figure)
                delete!(state, :plot_export_error)
            catch err
                message = "Export failed for $path: $(sprint(showerror, err))"
                state[:plot_export_error] = message
                @error "Plot export failed\n$message"
            end
        end
    end
    !can_export && ig.EndDisabled()

    ig.SameLine()
    if ig.Button("?##plot_help_$id")
        ig.OpenPopup("plot_help_popup_$id")
    end
    if ig.BeginItemTooltip()
        ig.TextUnformatted("Help")
        ig.EndTooltip()
    end
    if ig.BeginPopup("plot_help_popup_$id")
        ig.PushTextWrapPos(ig.GetFontSize() * 38.0)
        ig.TextUnformatted(PLOT_HELP_TEXT)
        ig.PopTextWrapPos()
        ig.EndPopup()
    end
    haskey(state, :plot_export_error) && _render_plot_error(state[:plot_export_error])
    return nothing
end

"""Render a figure, empty text, or error text for one plot target."""
function _render_plot_body!(ui_state::Dict{Symbol,Any}, state::Dict{Symbol,Any}, status::Symbol, figure_key::Symbol, makie_id::String)::Nothing
    if get(ui_state, :debug_plot_mode, false)
        ig.TextColored((0.2, 0.8, 0.2, 1.0), "Debug Plot Mode")
        ig.SameLine()
        _helpmarker("Debug mode bypasses cache and reads source files directly.")
    end

    figure = get(state, figure_key, nothing)
    if figure !== nothing
        _time!(ui_state, :makie_fig) do
            MakieFigure(makie_id, figure; auto_resize_x=true, auto_resize_y=true, tooltip=false)
        end
    elseif status == :needs_kind
        ig.TextDisabled("No plot kind selected")
    elseif status == :empty
        ig.TextDisabled("No measurement selected")
    else
        ig.TextColored((1.0, 0.4, 0.4, 1.0), "Plot generation failed")
        _render_plot_error(get(state, :plot_error, ""))
    end
    return nothing
end

"""Warm GLMakie once in a hidden tiny window to reduce the first visible plot stall."""
function _ensure_plot_runtime_warmed!(ui_state::Dict{Symbol,Any})::Nothing
    get(ui_state, :plot_runtime_warmed, false) && return nothing
    flags = ig.ImGuiWindowFlags_NoDecoration | ig.ImGuiWindowFlags_NoInputs |
            ig.ImGuiWindowFlags_NoBackground | ig.ImGuiWindowFlags_NoSavedSettings |
            ig.ImGuiWindowFlags_NoNav | ig.ImGuiWindowFlags_NoFocusOnAppearing
    ig.SetNextWindowPos((-10_000.0, -10_000.0), ig.ImGuiCond_Always)
    ig.SetNextWindowSize((8.0, 8.0), ig.ImGuiCond_Always)
    if ig.Begin("###plot_runtime_warmup", C_NULL, flags)
        fig = get(ui_state, :plot_runtime_warmup_figure, nothing)
        if fig === nothing
            fig = Figure(size=(64, 48))
            ax = Axis(fig[1, 1], xlabel="x", ylabel="y", title="warm")
            lines!(ax, [0.0, 1.0], [0.0, 1.0], color=:blue)
            ui_state[:plot_runtime_warmup_figure] = fig
        end
        MakieFigure("_plot_runtime_warmup", fig; auto_resize_x=false, auto_resize_y=false, tooltip=false, stats=false)
        ui_state[:plot_runtime_warmed] = true
    end
    ig.End()
    return nothing
end

"""Render the main plot window."""
function render_plot_window(ui_state::Dict{Symbol,Any})::Nothing
    selected = get(ui_state, :selected_measurements, MeasurementInfo[])
    measurements = get(ui_state, :main_plot_live, true) === true ? selected : get(ui_state, :main_plot_measurements, selected)
    measurements isa Vector{MeasurementInfo} || (measurements = MeasurementInfo[])
    get(ui_state, :main_plot_live, true) === true && (ui_state[:main_plot_measurements] = copy(measurements))
    if get(ui_state, :main_plot_live, true) === true
        if isempty(measurements) || any(m -> m.measurement_kind !== first(measurements).measurement_kind, measurements)
            delete!(ui_state, :main_plot_measurement_kind)
        else
            measurement_kind = first(measurements).measurement_kind
            if get(ui_state, :main_plot_measurement_kind, nothing) !== measurement_kind
                ui_state[:main_plot_measurement_kind] = measurement_kind
                prefs = get(ui_state, :plot_kind_by_measurement_kind, Dict{String,String}())
                remembered = prefs isa AbstractDict ? get(prefs, String(measurement_kind), "") : ""
                plot_kind = remembered isa AbstractString ? _plot_kind_from_name(remembered) : nothing
                plot_kind === nothing ? delete!(ui_state, :main_plot_kind) : (ui_state[:main_plot_kind] = plot_kind)
                delete!(ui_state, :_last_plot_key)
            end
        end
    end

    if ig.Begin("Plot Area")
        _render_plot_view!(
            ui_state,
            ui_state,
            measurements,
            "main",
            :plot_figure,
            "measurement_plot";
            live_key=:main_plot_live,
            live_default=true,
            plot_kind_key=:main_plot_kind,
            selected_measurements=selected,
        )
    end
    ig.End()
    return nothing
end

"""Render detached plot windows."""
function render_additional_plot_windows(ui_state::Dict{Symbol,Any})::Nothing
    open_plots = get(ui_state, :open_plot_windows, nothing)
    open_plots === nothing && return nothing
    selected = get(ui_state, :selected_measurements, MeasurementInfo[])
    kept = Vector{Dict{Symbol,Any}}()

    for entry in open_plots
        target_id = string(get!(entry, :target_id) do
            _next_plot_window_id!(ui_state)
        end)
        entry[:target_id] = target_id
        open_ref = Ref(true)
        window_id = replace(target_id, '/' => '_')

        if ig.Begin("Plot: $(get(entry, :title, "Plot"))###plot_window_$window_id", open_ref)
            measurements = get(entry, :live, false) === true ? selected : get(entry, :measurements, MeasurementInfo[])
            measurements isa Vector{MeasurementInfo} || (measurements = MeasurementInfo[])
            get(entry, :live, false) === true && (entry[:measurements] = copy(measurements))
            if get(entry, :debug, false) != get(ui_state, :debug_plot_mode, false)
                delete!(entry, :_last_plot_key)
            end
            _render_plot_view!(
                ui_state,
                entry,
                measurements,
                target_id,
                :figure,
                "measurement_plot_$window_id";
                live_key=:live,
                live_default=false,
                plot_kind_key=:plot_kind,
                selected_measurements=selected,
            )
        end
        ig.End()
        open_ref[] && push!(kept, entry)
    end
    ui_state[:open_plot_windows] = kept
    return nothing
end
