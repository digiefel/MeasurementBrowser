function _plot_kind_name(kind::Type{<:PlotKind})::String
    return String(nameof(kind))
end

function _selected_plot_kind!(ui_state)
    kinds = plot_kinds()
    isempty(kinds) && return nothing
    current = get(ui_state, :selected_plot_kind, nothing)
    if current isa Type && current <: PlotKind && current in kinds
        return current
    end
    ui_state[:selected_plot_kind] = first(kinds)
    return ui_state[:selected_plot_kind]
end

function _measurement_plot_window_entry(ui_state, measurement::MeasurementInfo)
    plot_kind = _selected_plot_kind!(ui_state)
    target_id = plot_kind === nothing ? measurement.unique_id :
        "$(measurement.unique_id)#$(nameof(plot_kind))"
    return Dict{Symbol,Any}(
        :target_id => target_id,
        :filepath => measurement.filepath,
        :measurement => measurement,
        :title => measurement.clean_title,
        :plot_kind => plot_kind,
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
        failed ? (entry[:plot_error] = message) : delete!(entry, :plot_error)
        entry[:debug] = job.debug
        return failed
    end
    return failed
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
                message = bt === nothing ? sprint(showerror, err) : sprint(showerror, err, bt)
                _finish_plot_job!(ui_state, job; error=message)
                ui_state[:plot_state] = :error
                @error "Plot drawing failed" exception = (err, bt)
            end
            _finalize_plot_job!(ui_state)
        elseif msg.kind == :canceled
            ui_state[:plot_state] = :canceled
            _finalize_plot_job!(ui_state)
        elseif msg.kind == :error
            job = get(ui_state, :active_plot_job, nothing)
            job isa PlotJob || continue
            message = msg.bt === nothing ? sprint(showerror, msg.error) : sprint(showerror, msg.error, msg.bt)
            _finish_plot_job!(ui_state, job; error=message)
            ui_state[:plot_state] = :error
            @error "Plot job failed" exception = (msg.error, msg.bt)
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

function _open_combined_plot_window!(ui_state, proj, measurements::Vector{MeasurementInfo}, plot_kind::Type{<:PlotKind})
    open_plots = get!(ui_state, :open_plot_windows) do
        Vector{Dict{Symbol,Any}}()
    end
    counter = get!(ui_state, :_combined_plot_counter) do
        0
    end
    ui_state[:_combined_plot_counter] = counter + 1
    target_id = "combined_$(counter + 1)"
    entry = Dict(
        :title => "Combined: $(_plot_kind_name(plot_kind))",
        :id => target_id,
        :target_id => target_id,
        :combined_kind => plot_kind,
        :measurements => measurements,
    )
    push!(open_plots, entry)
    _queue_plot_job!(ui_state, proj, measurements, plot_kind; target=:extra, target_id)
end

function render_plot_window(ui_state)
    proj = ui_state[:project]
    selected_measurements = get(ui_state, :selected_measurements, MeasurementInfo[])
    selected_kind = _selected_plot_kind!(ui_state)
    plot_measurements = MeasurementInfo[]
    plot_kind = selected_kind
    plot_key = nothing
    plot_status = :empty

    if isempty(selected_measurements)
        plot_status = :empty
    elseif plot_kind === nothing
        plot_status = :needs_kind
    else
        plot_measurements = selected_measurements
        plot_key = (
            nameof(plot_kind),
            sort([measurement.unique_id for measurement in plot_measurements]),
            _plot_parameters_key([measurement.parameters for measurement in plot_measurements]),
        )
        plot_status = :ready
    end

    generate_combined = get(ui_state, :generate_combined_plot, false)
    if generate_combined
        ui_state[:generate_combined_plot] = false
        _clear_plot_jobs!(ui_state)
    end
    should_queue = plot_status == :ready &&
        get(ui_state, :_last_plot_key, nothing) != plot_key &&
        (length(plot_measurements) == 1 || generate_combined)
    if should_queue
        _queue_plot_job!(ui_state, proj, plot_measurements, plot_kind; target=:main, target_id="main", plot_key)
    end

    if ig.Begin("Plot Area")
        if get(ui_state, :debug_plot_mode, false)
            ig.TextColored((0.2, 0.8, 0.2, 1.0), "Debug Plot Mode")
            ig.SameLine()
            _helpmarker("Debug mode is ON: plots have DEBUG flag.")
        end
        _render_plot_indicator!(ui_state)

        if haskey(ui_state, :plot_figure)
            f = ui_state[:plot_figure]
            _time!(ui_state, :makie_fig) do
                MakieFigure("measurement_plot", f; auto_resize_x=true, auto_resize_y=true)
            end
        else
            if _plot_target_loading(ui_state, :main; target_id="main")
                ig.TextDisabled(length(plot_measurements) > 1 ? "Loading combined plot..." : "Loading plot...")
            elseif plot_status == :needs_kind
                ig.TextDisabled("No plot kind selected")
            elseif plot_status == :empty
                ig.TextDisabled("No measurement selected")
            else
                ig.TextColored((1.0, 0.4, 0.4, 1.0), length(plot_measurements) > 1 ? "Combined plot generation failed" : "Plot generation failed")
                ig.Text("Check file format and measurement type")
                _render_plot_error_detail!(get(ui_state, :plot_error, ""))
            end
        end

        ig.Separator()
    end
    ig.End()
end

function render_combined_plots_window(ui_state)
    proj = ui_state[:project]
    if ig.Begin("Combined Plots")
        current_kind = _selected_plot_kind!(ui_state)
        kind_label = current_kind === nothing ? "None" : _plot_kind_name(current_kind)

        if ig.BeginCombo("Plot Kind", kind_label)
            for plot_kind in plot_kinds()
                if ig.Selectable(_plot_kind_name(plot_kind), current_kind === plot_kind)
                    ui_state[:selected_plot_kind] = plot_kind
                end
            end
            ig.EndCombo()
        end

        selected_measurements = get(ui_state, :selected_measurements, MeasurementInfo[])
        can_generate = current_kind !== nothing && !isempty(selected_measurements)
        ig.Separator()

        if !isempty(selected_measurements)
            ig.TextColored((0.6, 0.8, 1.0, 1.0), "Selected: $(length(selected_measurements)) measurements")
        else
            ig.TextDisabled("No measurements selected")
        end

        if can_generate
            ig.TextColored((0.2, 0.8, 0.2, 1.0), "Ready")
        elseif length(selected_measurements) > 1
            ig.TextColored((0.8, 0.8, 0.2, 1.0), "Select a plot kind")
        end

        ig.Separator()
        if !can_generate
            ig.BeginDisabled()
        end
        if ig.Button("Generate Combined Plot", (-1, 0))
            ui_state[:generate_combined_plot] = true
        end
        if can_generate && ig.BeginPopupContextItem()
            if ig.MenuItem("Open in New Window")
                _open_combined_plot_window!(ui_state, proj, selected_measurements, current_kind)
            end
            ig.EndPopup()
        end
        if !can_generate
            ig.EndDisabled()
        end
    end
    ig.End()
end

# Render any additional plot windows opened via right-click context menu.
function render_additional_plot_windows(ui_state)
    proj = ui_state[:project]
    open_plots = get(ui_state, :open_plot_windows, nothing)
    open_plots === nothing && return
    isempty(open_plots) && return
    to_keep = Vector{Dict{Symbol,Any}}()
    for entry in open_plots
        is_combined = haskey(entry, :combined_kind)
        filepath = get(entry, :filepath, "")
        measurements = is_combined ? get(entry, :measurements, MeasurementInfo[]) :
            (get(entry, :measurement, nothing) isa MeasurementInfo ? [entry[:measurement]] : MeasurementInfo[])
        isempty(measurements) && continue

        plot_kind = is_combined ? entry[:combined_kind] : entry[:plot_kind]
        title = get(entry, :title, is_combined ? "Combined Plot" : basename(filepath))
        target_id = string(get(entry, :target_id, get(entry, :id, filepath)))
        entry[:target_id] = target_id
        if plot_kind === nothing
            entry[:plot_error] = "No plot kind is available for this measurement."
        end

        existing_debug = get(entry, :debug, false)
        has_failure = haskey(entry, :plot_error)
        refresh = (!haskey(entry, :figure) && !has_failure) ||
            existing_debug != get(ui_state, :debug_plot_mode, false)
        if refresh
            _queue_plot_job!(ui_state, proj, measurements, plot_kind; target=:extra, target_id)
        end

        open_ref = Ref(true)
        window_id = replace(target_id, '/' => '_')
        if ig.Begin("Plot: $title###plot_window_$window_id", open_ref)
            if haskey(entry, :figure)
                f = entry[:figure]
                _time!(ui_state, :makie_fig) do
                    MakieFigure("measurement_plot_$window_id", f; auto_resize_x=true, auto_resize_y=true)
                end
            elseif _plot_target_loading(ui_state, :extra; target_id=target_id)
                ig.TextDisabled("Loading plot...")
            elseif haskey(entry, :plot_error)
                ig.TextColored((1.0, 0.4, 0.4, 1.0), "Plot generation failed")
                _render_plot_error_detail!(entry[:plot_error])
            else
                ig.Text("No plot available")
            end
            ig.Separator()
            ig.TextDisabled(is_combined ? title : basename(filepath))
        end
        ig.End()
        open_ref[] && push!(to_keep, entry)
    end
    ui_state[:open_plot_windows] = to_keep
end
