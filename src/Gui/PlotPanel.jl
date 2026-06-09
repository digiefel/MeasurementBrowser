function _plot_kind_name(kind::Type{<:PlotKind})::String
    return String(nameof(kind))
end

function _plot_kind_from_state!(state::Dict{Symbol,Any}, key::Symbol)::Union{Nothing,Type{<:PlotKind}}
    kinds = plot_kinds()
    isempty(kinds) && return nothing
    current = get(state, key, nothing)
    if current isa Type && current <: PlotKind && current in kinds
        return current
    end
    delete!(state, key)
    return nothing
end

function _next_plot_window_id!(ui_state::Dict{Symbol,Any})::String
    counter = get(ui_state, :_plot_window_counter, 0) + 1
    ui_state[:_plot_window_counter] = counter
    return "plot_$counter"
end

function _plot_view_key(
    plot_kind::Type{<:PlotKind},
    measurements::Vector{MeasurementInfo},
)::Tuple
    return (
        nameof(plot_kind),
        sort([measurement.unique_id for measurement in measurements]),
        _plot_parameters_key([measurement.parameters for measurement in measurements]),
    )
end

function _measurement_plot_window_entry(
    ui_state::Dict{Symbol,Any},
    measurement::MeasurementInfo,
)::Dict{Symbol,Any}
    plot_kind = _plot_kind_from_state!(ui_state, :main_plot_kind)
    return Dict{Symbol,Any}(
        :target_id => _next_plot_window_id!(ui_state),
        :filepath => measurement.filepath,
        :measurements => [measurement],
        :title => measurement.clean_title,
        :plot_kind => plot_kind,
        :live => false,
        :params => _measurement_parameters(measurement),
    )
end

function _plot_running(ui_state)
    return get(ui_state, :plot_state, :idle) in (:loading, :analyzing, :drawing, :canceling)
end

function _clear_plot_jobs!(ui_state)
    ui_state[:plot_state] = :idle
    ui_state[:plot_error] = ""
    ui_state[:plot_events] = nothing
    ui_state[:plot_cancel_token] = nothing
    ui_state[:active_plot_job] = nothing
    ui_state[:pending_plot_job] = nothing
    ui_state[:active_plot_id] = get(ui_state, :active_plot_id, 0) + 1
end

function _request_plot_cancel!(ui_state)
    _plot_running(ui_state) || return
    token = get(ui_state, :plot_cancel_token, nothing)
    token !== nothing && Base.Threads.atomic_xchg!(token, true)
    ui_state[:plot_state] = :canceling
end

function _plot_job(
    ui_state,
    proj,
    measurements::Vector{MeasurementInfo},
    plot_kind::Type{<:PlotKind},
    target::Symbol;
    target_id::String,
    plot_key=nothing,
)
    isempty(measurements) && error("Plot job requires at least one measurement")
    debug = get(ui_state, :debug_plot_mode, false)
    device_params = [merge(measurement.device_info.parameters, measurement.parameters) for measurement in measurements]
    return PlotJob(
        _plot_job_key(proj, target, target_id, measurements, plot_kind, device_params, debug),
        plot_key,
        proj,
        target,
        target_id,
        measurements,
        plot_kind,
        device_params,
        debug,
    )
end

function _launch_plot_job!(ui_state, job::PlotJob)
    plot_id = get(ui_state, :plot_seq, 0) + 1
    ui_state[:plot_seq] = plot_id
    ui_state[:active_plot_id] = plot_id
    ui_state[:active_plot_job] = job
    ui_state[:plot_state] = :loading
    ui_state[:plot_error] = ""
    ui_state[:pending_plot_job] = nothing

    events = Channel{NamedTuple}(8)
    ui_state[:plot_events] = events
    cancel_token = Base.Threads.Atomic{Bool}(false)
    ui_state[:plot_cancel_token] = cancel_token

    Base.Threads.@spawn begin
        try
            data = _run_plot_job(job, () -> cancel_token[])
            cancel_token[] && throw(JobCancelled())
            put!(events, (kind=:data, plot_id=plot_id, data=data))
        catch err
            if err isa JobCancelled
                put!(events, (kind=:canceled, plot_id=plot_id))
            else
                put!(events, (kind=:error, plot_id=plot_id, error=err, bt=catch_backtrace()))
            end
        finally
            close(events)
        end
    end
end

function _queue_plot_job!(ui_state, job::PlotJob)
    if _plot_running(ui_state)
        active_job = get(ui_state, :active_plot_job, nothing)
        pending_job = get(ui_state, :pending_plot_job, nothing)
        if active_job isa PlotJob && active_job.job_key == job.job_key
            return
        end
        if pending_job isa PlotJob && pending_job.job_key == job.job_key
            return
        end
        ui_state[:pending_plot_job] = job
        _request_plot_cancel!(ui_state)
        return
    end
    _launch_plot_job!(ui_state, job)
end

function _queue_plot_job!(ui_state, proj, measurements::Vector{MeasurementInfo}, plot_kind::Type{<:PlotKind};
                          target::Symbol, target_id::String, plot_key=nothing)
    job = _plot_job(ui_state, proj, measurements, plot_kind, target; target_id, plot_key)
    _queue_plot_job!(ui_state, job)
    return job
end

function _finalize_plot_job!(ui_state)
    ui_state[:plot_events] = nothing
    ui_state[:plot_cancel_token] = nothing
    ui_state[:active_plot_job] = nothing
    pending_job = get(ui_state, :pending_plot_job, nothing)
    ui_state[:pending_plot_job] = nothing
    if pending_job isa PlotJob
        _launch_plot_job!(ui_state, pending_job)
        return
    end
end

function _finish_plot_job!(ui_state, job::PlotJob; fig=nothing, error::AbstractString="")
    message = isempty(error) && fig === nothing ? "Plot renderer returned no figure." : String(error)
    failed = fig === nothing || !isempty(message)

    if job.target == :main
        fig === nothing ? delete!(ui_state, :plot_figure) : (ui_state[:plot_figure] = fig)
        job.plot_key === nothing ? delete!(ui_state, :_last_plot_key) : (ui_state[:_last_plot_key] = job.plot_key)
        ui_state[:plot_error] = failed ? message : ""
        return failed
    end

    open_plots = get(ui_state, :open_plot_windows, nothing)
    open_plots === nothing && return failed

    for entry in open_plots
        entry_target_id = get(entry, :target_id, get(entry, :filepath, ""))
        entry_target_id == job.target_id || continue
        fig === nothing ? delete!(entry, :figure) : (entry[:figure] = fig)
        job.plot_key === nothing ? delete!(entry, :_last_plot_key) : (entry[:_last_plot_key] = job.plot_key)
        failed ? (entry[:plot_error] = message) : delete!(entry, :plot_error)
        entry[:debug] = job.debug
        return failed
    end
    return failed
end

function _plot_measurement_context(measurements::Vector{MeasurementInfo})::String
    lines = String["Measurements: $(length(measurements))"]
    shown = min(length(measurements), 3)
    for measurement in first(measurements, shown)
        push!(lines, "Measurement: $(measurement.clean_title) ($(measurement.measurement_kind))")
        push!(lines, "File: $(measurement.filepath)")
    end
    length(measurements) > shown && push!(lines, "More files: $(length(measurements) - shown)")
    return join(lines, "\n")
end

function _plot_stack_noise(file::AbstractString)::Bool
    return occursin("GlfwOpenGLBackend", file) ||
        occursin("CImGui", file) ||
        occursin("GLFW", file) ||
        endswith(file, "client.jl")
end

function _plot_backtrace_summary(bt)::String
    bt === nothing && return ""
    lines = String[]
    for frame in stacktrace(bt)
        file = String(frame.file)
        _plot_stack_noise(file) && continue
        push!(lines, "  $(frame.func) at $(basename(file)):$(frame.line)")
        length(lines) >= 10 && break
    end
    isempty(lines) && return ""
    return "Stack:\n" * join(lines, "\n")
end

function _format_plot_error(job::PlotJob, err, bt; phase::AbstractString)::String
    parts = String[
        "Plot failed while $phase: $(sprint(showerror, err))",
        "Project: $(project_name(job.project))",
        "Plot: $(nameof(job.plot_kind))",
        _plot_measurement_context(job.measurements),
    ]
    stack = _plot_backtrace_summary(bt)
    !isempty(stack) && push!(parts, stack)
    return join(parts, "\n")
end

function _poll_plot_events!(ui_state)
    events = get(ui_state, :plot_events, nothing)
    events === nothing && return

    while isready(events)
        msg = take!(events)
        msg.plot_id == get(ui_state, :active_plot_id, 0) || continue

        if msg.kind == :data
            job = get(ui_state, :active_plot_job, nothing)
            job isa PlotJob || continue
            ui_state[:plot_state] = :drawing
            try
                fig = _draw_plot_job(job, msg.data)
                ui_state[:plot_state] = _finish_plot_job!(ui_state, job; fig) ? :error : :done
            catch err
                bt = catch_backtrace()
                message = _format_plot_error(job, err, bt; phase="drawing")
                _finish_plot_job!(ui_state, job; error=message)
                ui_state[:plot_state] = :error
                @error "Plot drawing failed" details = message
            end
            _finalize_plot_job!(ui_state)
        elseif msg.kind == :canceled
            ui_state[:plot_state] = :canceled
            _finalize_plot_job!(ui_state)
        elseif msg.kind == :error
            job = get(ui_state, :active_plot_job, nothing)
            job isa PlotJob || continue
            message = _format_plot_error(job, msg.error, msg.bt; phase="loading data")
            _finish_plot_job!(ui_state, job; error=message)
            ui_state[:plot_state] = :error
            @error "Plot job failed" details = message
            _finalize_plot_job!(ui_state)
        end
    end
end

function _plot_target_loading(ui_state, target::Symbol; target_id::String="")
    _plot_running(ui_state) || return false
    active_job = get(ui_state, :active_plot_job, nothing)
    pending_job = get(ui_state, :pending_plot_job, nothing)
    for job in (active_job, pending_job)
        job isa PlotJob || continue
        job.target == target || continue
        job.target_id == target_id && return true
    end
    return false
end

function _render_plot_indicator!(ui_state)
    _plot_running(ui_state) || return

    state = get(ui_state, :plot_state, :idle)
    if state == :loading
        ig.TextDisabled("Plot: loading data...")
    elseif state == :analyzing
        ig.TextDisabled("Plot: analyzing data...")
    elseif state == :drawing
        ig.TextDisabled("Plot: drawing figure...")
    elseif state == :canceling
        ig.TextDisabled("Plot: canceling...")
    end
end

function _render_plot_error_detail!(message)
    isempty(message) && return
    lines = split(String(message), '\n')
    summary = isempty(lines) ? "" : first(lines)
    if length(summary) > 140
        summary = first(summary, 137) * "..."
    end
    !isempty(summary) && ig.TextDisabled(summary)
    if ig.BeginItemTooltip()
        ig.PushTextWrapPos(ig.GetFontSize() * 45.0)
        ig.TextUnformatted(String(message))
        ig.PopTextWrapPos()
        ig.EndTooltip()
    end
end

function _render_plot_kind_combo!(
    state::Dict{Symbol,Any},
    key::Symbol,
    label::AbstractString,
)::Bool
    current = _plot_kind_from_state!(state, key)
    preview = current === nothing ? "Choose plot kind" : _plot_kind_name(current)
    changed = false

    ig.SetNextItemWidth(240)
    if ig.BeginCombo(String(label), preview)
        for plot_kind in plot_kinds()
            if ig.Selectable(_plot_kind_name(plot_kind), current === plot_kind)
                state[key] = plot_kind
                changed = true
            end
        end
        ig.EndCombo()
    end
    return changed
end

function _render_live_checkbox!(
    state::Dict{Symbol,Any},
    key::Symbol,
    default::Bool,
    label::AbstractString,
)::Bool
    live = get(state, key, default) === true
    changed = @c ig.Checkbox(String(label), &live)
    state[key] = live
    return changed
end

function _plot_help_text()::String
    return """
    Plot controls:
    Live keeps this plot attached to the current browser selection.
    Detach opens the current view as an independent plot window.
    Export saves the current figure.

    Figure controls:
    Scroll to zoom.
    Click and drag to zoom into a rectangle.
    Right click and drag to pan.
    Use x/y plus scroll to zoom one axis.
    Ctrl + left click resets the limits.
    Right click inside the plot for axis scale settings.
    """
end

function _render_plot_help_button!(id::AbstractString)::Nothing
    popup_id = "plot_help_popup_$id"
    button_width = ig.CalcTextSize("?").x + 18.0f0
    remaining = ig.GetContentRegionAvail().x
    ig.SameLine()
    ig.SetCursorPosX(ig.GetCursorPosX() + max(0.0f0, remaining - button_width))
    if ig.Button("?##plot_help_$id")
        ig.OpenPopup(popup_id)
    end
    if ig.BeginItemTooltip()
        ig.TextUnformatted("Help")
        ig.EndTooltip()
    end
    if ig.BeginPopup(popup_id)
        ig.PushTextWrapPos(ig.GetFontSize() * 38.0)
        ig.TextUnformatted(_plot_help_text())
        ig.PopTextWrapPos()
        ig.EndPopup()
    end
    return nothing
end

function _plot_export_path(default_name::AbstractString)::String
    path = save_file(String(default_name); filterlist="png,jpg,jpeg,svg,pdf")
    isempty(path) && return ""
    _, ext = splitext(path)
    isempty(ext) && return path * ".png"
    return path
end

function _export_plot_figure!(
    state::Dict{Symbol,Any},
    figure::Figure,
    default_name::AbstractString,
)::Nothing
    path = _plot_export_path(default_name)
    isempty(path) && return nothing

    try
        Makie.save(path, figure)
        delete!(state, :plot_export_error)
    catch err
        message = "Export failed for $path: $(sprint(showerror, err))"
        state[:plot_export_error] = message
        @error "Plot export failed" details = message
    end
    return nothing
end

function _plot_export_name(
    plot_kind::Union{Nothing,Type{<:PlotKind}},
    measurements::Vector{MeasurementInfo},
)::String
    kind_name = plot_kind === nothing ? "plot" : lowercase(String(nameof(plot_kind)))
    if length(measurements) == 1
        stem = splitext(basename(only(measurements).filepath))[1]
        return "$stem-$kind_name.png"
    end
    return "$(length(measurements))-measurements-$kind_name.png"
end

function _detach_plot_window!(
    ui_state::Dict{Symbol,Any},
    measurements::Vector{MeasurementInfo},
    plot_kind::Type{<:PlotKind},
)::Nothing
    open_plots = get!(ui_state, :open_plot_windows) do
        Vector{Dict{Symbol,Any}}()
    end
    target_id = _next_plot_window_id!(ui_state)
    title = length(measurements) == 1 ?
        only(measurements).clean_title :
        "$(length(measurements)) measurements"
    push!(open_plots, Dict{Symbol,Any}(
        :title => title,
        :target_id => target_id,
        :measurements => copy(measurements),
        :plot_kind => plot_kind,
        :live => false,
    ))
    return nothing
end

function _main_plot_measurements!(ui_state::Dict{Symbol,Any})::Vector{MeasurementInfo}
    selected = get(ui_state, :selected_measurements, MeasurementInfo[])
    if get(ui_state, :main_plot_live, true) === true
        ui_state[:main_plot_measurements] = copy(selected)
        return selected
    end

    measurements = get(ui_state, :main_plot_measurements, nothing)
    if measurements isa Vector{MeasurementInfo}
        return measurements
    end
    ui_state[:main_plot_measurements] = copy(selected)
    return selected
end

function _entry_plot_measurements!(
    entry::Dict{Symbol,Any},
    selected::Vector{MeasurementInfo},
)::Vector{MeasurementInfo}
    if get(entry, :live, false) === true
        entry[:measurements] = copy(selected)
        return selected
    end

    measurements = get(entry, :measurements, MeasurementInfo[])
    if measurements isa Vector{MeasurementInfo}
        return measurements
    end
    return MeasurementInfo[]
end

function _render_main_plot_toolbar!(
    ui_state::Dict{Symbol,Any},
    measurements::Vector{MeasurementInfo},
    plot_kind::Union{Nothing,Type{<:PlotKind}},
)::Nothing
    if _render_plot_kind_combo!(ui_state, :main_plot_kind, "##main_plot_kind")
        delete!(ui_state, :_last_plot_key)
        ui_state[:plot_error] = ""
    end

    ig.SameLine()
    if _render_live_checkbox!(ui_state, :main_plot_live, true, "Live##main_plot_live")
        if get(ui_state, :main_plot_live, true) === false
            ui_state[:main_plot_measurements] = copy(get(ui_state, :selected_measurements, MeasurementInfo[]))
        end
        delete!(ui_state, :_last_plot_key)
    end

    can_detach = plot_kind !== nothing && !isempty(measurements)
    ig.SameLine()
    !can_detach && ig.BeginDisabled()
    if ig.Button("Detach") && can_detach
        _detach_plot_window!(ui_state, measurements, plot_kind)
    end
    !can_detach && ig.EndDisabled()

    can_export = haskey(ui_state, :plot_figure)
    ig.SameLine()
    !can_export && ig.BeginDisabled()
    if ig.Button("Export") && can_export
        _export_plot_figure!(
            ui_state,
            ui_state[:plot_figure],
            _plot_export_name(plot_kind, measurements),
        )
    end
    !can_export && ig.EndDisabled()

    _render_plot_help_button!("main")
    if haskey(ui_state, :plot_export_error)
        _render_plot_error_detail!(ui_state[:plot_export_error])
    end
    return nothing
end

function _render_extra_plot_toolbar!(
    entry::Dict{Symbol,Any},
    selected_measurements::Vector{MeasurementInfo},
    id::AbstractString,
)::Nothing
    if _render_plot_kind_combo!(entry, :plot_kind, "##plot_kind_$id")
        delete!(entry, :_last_plot_key)
        delete!(entry, :plot_error)
    end

    ig.SameLine()
    if _render_live_checkbox!(entry, :live, false, "Live##plot_live_$id")
        if get(entry, :live, false) === false
            entry[:measurements] = copy(selected_measurements)
        end
        delete!(entry, :_last_plot_key)
    end

    can_export = haskey(entry, :figure)
    ig.SameLine()
    !can_export && ig.BeginDisabled()
    if ig.Button("Export##plot_export_$id") && can_export
        plot_kind = _plot_kind_from_state!(entry, :plot_kind)
        measurements = _entry_plot_measurements!(entry, selected_measurements)
        _export_plot_figure!(
            entry,
            entry[:figure],
            _plot_export_name(plot_kind, measurements),
        )
    end
    !can_export && ig.EndDisabled()

    _render_plot_help_button!(id)
    if haskey(entry, :plot_export_error)
        _render_plot_error_detail!(entry[:plot_export_error])
    end
    return nothing
end

function _plot_runtime_warmup_figure()
    fig = Figure(size=(64, 48))
    ax1 = Axis(fig[1, 1], xlabel="x", ylabel="y", title="warm")
    ax2 = Axis(fig[1, 2], xlabel="x", ylabel="y")
    lineplot = lines!(ax1, [0.0, 1.0], [0.0, 1.0], color=:blue, linewidth=2)
    scatter!(ax1, [0.5], [0.5], color=:red, markersize=8)
    lines!(ax2, [0.0, 1.0], [1.0, 0.0], color=:green, linewidth=2)
    Legend(fig[0, 1], [lineplot], ["warm"], tellwidth=false, tellheight=false)
    return fig
end

function _ensure_plot_runtime_warmed!(ui_state)
    get(ui_state, :plot_runtime_warmed, false) && return

    flags = ig.ImGuiWindowFlags_NoDecoration |
            ig.ImGuiWindowFlags_NoInputs |
            ig.ImGuiWindowFlags_NoBackground |
            ig.ImGuiWindowFlags_NoSavedSettings |
            ig.ImGuiWindowFlags_NoNav |
            ig.ImGuiWindowFlags_NoFocusOnAppearing
    ig.SetNextWindowPos((-10_000.0, -10_000.0), ig.ImGuiCond_Always)
    ig.SetNextWindowSize((8.0, 8.0), ig.ImGuiCond_Always)

    if ig.Begin("###plot_runtime_warmup", C_NULL, flags)
        fig = get(ui_state, :plot_runtime_warmup_figure, nothing)
        if fig === nothing
            fig = _plot_runtime_warmup_figure()
            ui_state[:plot_runtime_warmup_figure] = fig
        end
        MakieFigure("_plot_runtime_warmup", fig; auto_resize_x=false, auto_resize_y=false, tooltip=false, stats=false)
        ui_state[:plot_runtime_warmed] = true
    end
    ig.End()
end

function render_plot_window(ui_state::Dict{Symbol,Any})::Nothing
    proj = ui_state[:project]
    plot_measurements = _main_plot_measurements!(ui_state)
    plot_kind = _plot_kind_from_state!(ui_state, :main_plot_kind)
    plot_key = nothing
    plot_status = :empty

    if isempty(plot_measurements)
        plot_status = :empty
    elseif plot_kind === nothing
        plot_status = :needs_kind
    else
        plot_key = _plot_view_key(plot_kind, plot_measurements)
        plot_status = :ready
    end

    should_queue = plot_status == :ready &&
        get(ui_state, :_last_plot_key, nothing) != plot_key
    if should_queue
        _queue_plot_job!(ui_state, proj, plot_measurements, plot_kind; target=:main, target_id="main", plot_key)
    elseif plot_status != :ready && get(ui_state, :main_plot_live, true) === true
        delete!(ui_state, :plot_figure)
        delete!(ui_state, :_last_plot_key)
    end

    if ig.Begin("Plot Area")
        _render_main_plot_toolbar!(ui_state, plot_measurements, plot_kind)
        ig.Separator()

        if get(ui_state, :debug_plot_mode, false)
            ig.TextColored((0.2, 0.8, 0.2, 1.0), "Debug Plot Mode")
            ig.SameLine()
            _helpmarker("Debug mode is ON: plots have DEBUG flag.")
        end
        _render_plot_indicator!(ui_state)

        if haskey(ui_state, :plot_figure)
            f = ui_state[:plot_figure]
            _time!(ui_state, :makie_fig) do
                MakieFigure("measurement_plot", f; auto_resize_x=true, auto_resize_y=true, tooltip=false)
            end
        else
            if _plot_target_loading(ui_state, :main; target_id="main")
                ig.TextDisabled("Loading plot...")
            elseif plot_status == :needs_kind
                ig.TextDisabled("No plot kind selected")
            elseif plot_status == :empty
                ig.TextDisabled("No measurement selected")
            else
                ig.TextColored((1.0, 0.4, 0.4, 1.0), "Plot generation failed")
                ig.Text("Check file format and measurement type")
                _render_plot_error_detail!(get(ui_state, :plot_error, ""))
            end
        end
    end
    ig.End()
    return nothing
end

function render_additional_plot_windows(ui_state::Dict{Symbol,Any})::Nothing
    proj = ui_state[:project]
    open_plots = get(ui_state, :open_plot_windows, nothing)
    open_plots === nothing && return
    isempty(open_plots) && return
    selected_measurements = get(ui_state, :selected_measurements, MeasurementInfo[])
    to_keep = Vector{Dict{Symbol,Any}}()
    for entry in open_plots
        if !haskey(entry, :target_id)
            entry[:target_id] = _next_plot_window_id!(ui_state)
        end
        target_id = string(entry[:target_id])
        entry[:target_id] = target_id
        open_ref = Ref(true)
        window_id = replace(target_id, '/' => '_')
        title = get(entry, :title, "Plot")
        if ig.Begin("Plot: $title###plot_window_$window_id", open_ref)
            _render_extra_plot_toolbar!(entry, selected_measurements, window_id)
            ig.Separator()

            measurements = _entry_plot_measurements!(entry, selected_measurements)
            plot_kind = _plot_kind_from_state!(entry, :plot_kind)
            plot_status = isempty(measurements) ? :empty :
                plot_kind === nothing ? :needs_kind : :ready
            plot_key = plot_status == :ready ? _plot_view_key(plot_kind, measurements) : nothing
            debug_changed = get(entry, :debug, false) != get(ui_state, :debug_plot_mode, false)
            should_queue = plot_status == :ready &&
                (get(entry, :_last_plot_key, nothing) != plot_key || debug_changed)
            if should_queue
                _queue_plot_job!(ui_state, proj, measurements, plot_kind; target=:extra, target_id, plot_key)
            elseif plot_status != :ready && get(entry, :live, false) === true
                delete!(entry, :figure)
                delete!(entry, :_last_plot_key)
            end

            if haskey(entry, :figure)
                f = entry[:figure]
                _time!(ui_state, :makie_fig) do
                    MakieFigure("measurement_plot_$window_id", f; auto_resize_x=true, auto_resize_y=true, tooltip=false)
                end
            elseif _plot_target_loading(ui_state, :extra; target_id=target_id)
                ig.TextDisabled("Loading plot...")
            elseif plot_status == :needs_kind
                ig.TextDisabled("No plot kind selected")
            elseif plot_status == :empty
                ig.TextDisabled("No measurement selected")
            elseif haskey(entry, :plot_error)
                ig.TextColored((1.0, 0.4, 0.4, 1.0), "Plot generation failed")
                _render_plot_error_detail!(entry[:plot_error])
            else
                ig.Text("No plot available")
            end
        end
        ig.End()
        open_ref[] && push!(to_keep, entry)
    end
    ui_state[:open_plot_windows] = to_keep
    return nothing
end
