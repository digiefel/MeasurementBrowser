function _measurement_plot_window_entry(measurement::MeasurementInfo)
    return Dict{Symbol,Any}(
        :target_id => measurement.unique_id,
        :filepath => measurement.filepath,
        :measurement => measurement,
        :title => measurement.clean_title,
        :measurement_kind => measurement.measurement_kind,
        :params => _measurement_parameters(measurement),
    )
end

function _cache_plot_version(ui_state)
    try
        identity = get(ui_state, :cache_identity, nothing)
        identity isa ProjectCacheIdentity ||
            error("Plot job requires an active HDF5 cache identity")
        status = get(ui_state, :cache_status, nothing)
        status isa ProjectCacheStatus ||
            error("Plot job requires a loaded HDF5 cache status")
        if get(ui_state, :warned_about_cache_version, false)
            @debug "Successfully determined cache version for plot job after previous failure" cache_id=identity.cache_id
            ui_state[:warned_about_cache_version] = false
        end
        return (
            identity.cache_id,
            identity.cache_path,
            status.fresh_files,
            status.stale_files,
            status.new_files,
            status.deleted_files,
            status.error_files,
        )
    catch err
        if !get(ui_state, :warned_about_cache_version, false)
            @warn "Failed to cache plot job" exception=(err, catch_backtrace())
            ui_state[:warned_about_cache_version] = true
        end
        return nothing
    end
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
    plot_kind::Symbol,
    target::Symbol;
    target_id::String,
    plot_key=nothing,
)
    isempty(measurements) && error("Plot job requires at least one measurement")
    debug = get(ui_state, :debug_plot_mode, false)
    cache_identity = get(ui_state, :cache_identity, nothing)
    cache_identity = cache_identity isa ProjectCacheIdentity ? cache_identity : nothing
    cache_version = _cache_plot_version(ui_state)
    device_params = [merge(measurement.device_info.parameters, measurement.parameters) for measurement in measurements]
    return PlotJob(
        _plot_job_key(proj, target, target_id, measurements, plot_kind, device_params, cache_version, debug),
        plot_key,
        proj,
        target,
        target_id,
        measurements,
        plot_kind,
        device_params,
        cache_identity,
        cache_version,
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
            cancel_token[] && throw(PlotCancelled())
            put!(events, (kind=:data, plot_id=plot_id, data=data))
        catch err
            if err isa PlotCancelled
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

function _queue_plot_job!(ui_state, proj, measurements::Vector{MeasurementInfo}, plot_kind::Symbol;
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
        entry[:cache_version] = job.cache_version
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

function _compatible_measurements(proj, measurements::Vector{MeasurementInfo}, combined_type)
    combined_type === nothing && return MeasurementInfo[]
    return filter(m -> m.measurement_kind in compatible_kinds(proj, combined_type), measurements)
end

function _open_combined_plot_window!(ui_state, proj, measurements::Vector{MeasurementInfo}, combined_type)
    open_plots = get!(ui_state, :open_plot_windows) do
        Vector{Dict{Symbol,Any}}()
    end
    counter = get!(ui_state, :_combined_plot_counter) do
        0
    end
    ui_state[:_combined_plot_counter] = counter + 1
    target_id = "combined_$(counter + 1)"
    entry = Dict(
        :title => "Combined: $(string(combined_type))",
        :id => target_id,
        :target_id => target_id,
        :combined_kind => combined_type,
        :measurements => measurements,
    )
    push!(open_plots, entry)
    _queue_plot_job!(ui_state, proj, measurements, combined_type; target=:extra, target_id)
end

function render_plot_window(ui_state)
    proj = ui_state[:project]
    selected_measurements = get(ui_state, :selected_measurements, MeasurementInfo[])
    combined_type = get(ui_state, :combined_plot_type, nothing)
    plot_measurements = MeasurementInfo[]
    plot_kind = nothing
    plot_key = nothing
    plot_status = :empty

    if length(selected_measurements) == 1
        measurement = only(selected_measurements)
        plot_measurements = selected_measurements
        plot_kind = measurement.measurement_kind
        plot_key = (
            measurement.unique_id,
            _cache_plot_version(ui_state),
            _plot_parameters_key([measurement.parameters]),
        )
        plot_status = :ready
    elseif length(selected_measurements) > 1
        if combined_type === nothing
            plot_status = :needs_kind
        else
            plot_measurements = _compatible_measurements(proj, selected_measurements, combined_type)
            plot_kind = combined_type
            plot_key = (combined_type, sort([m.filepath for m in plot_measurements]))
            plot_status = length(plot_measurements) >= 2 ? :ready : :incompatible
        end
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
    elseif plot_status == :incompatible
        delete!(ui_state, :plot_figure)
        ui_state[:_last_plot_key] = plot_key
        ui_state[:plot_error] = ""
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
            if plot_status == :incompatible
                ig.TextColored((1.0, 0.6, 0.2, 1.0), "Not enough compatible measurements for $(combined_type)")
                ig.Text("Selected: $(length(selected_measurements)), Compatible: $(length(plot_measurements))")
                if combined_type === :tlm_analysis
                    ig.TextDisabled("TLM Analysis requires ≥2 TLM 4-point measurements")
                elseif combined_type === :pund_fatigue
                    ig.TextDisabled("PUND Fatigue requires ≥2 PUND measurements")
                end
            elseif _plot_target_loading(ui_state, :main; target_id="main")
                ig.TextDisabled(length(plot_measurements) > 1 ? "Loading combined plot..." : "Loading plot...")
            elseif plot_status == :needs_kind
                ig.TextColored((0.6, 0.8, 1.0, 1.0), "Multiple measurements selected ($(length(selected_measurements)))")
                ig.Text("Choose a combined plot type and click Generate in the Combined Plots window")
            elseif plot_status == :empty
                ig.TextDisabled("Select measurements from the Hierarchy panel")
                ig.TextDisabled("• Single measurement: regular plot")
                ig.TextDisabled("• Multiple measurements + plot type: combined plot")
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
        ig.Text("Select Combined Plot Type:")

        ig.TextDisabled("1. Select measurements from the Hierarchy panel")
        ig.TextDisabled("2. Choose a plot type below")
        ig.TextDisabled("3. Click Generate to create combined plot")
        ig.Separator()

        current_type = get(ui_state, :combined_plot_type, nothing)
        type_label = current_type === nothing ? "None" : string(current_type)

        if ig.BeginCombo("Plot Type", type_label)
            for (type_key, type_name, type_description) in combined_plot_types(proj)
                if ig.Selectable(type_name, current_type === type_key)
                    ui_state[:combined_plot_type] = type_key
                end
                if type_key !== nothing
                    ig.SameLine()
                    _helpmarker(type_description)
                end
            end
            ig.EndCombo()
        end

        selected_measurements = get(ui_state, :selected_measurements, MeasurementInfo[])
        compatible_measurements = _compatible_measurements(proj, selected_measurements, current_type)
        compatible_count = length(compatible_measurements)
        can_generate = current_type !== nothing && compatible_count >= 2
        ig.Separator()

        if !isempty(selected_measurements)
            ig.TextColored((0.6, 0.8, 1.0, 1.0), "Selected: $(length(selected_measurements)) measurements")
        else
            ig.TextDisabled("No measurements selected")
        end

        if can_generate
            ig.TextColored((0.2, 0.8, 0.2, 1.0), "Ready: $compatible_count compatible measurements")
        elseif current_type !== nothing && length(selected_measurements) > 1
            ig.TextColored((1.0, 0.6, 0.2, 1.0), "Warning: Only $compatible_count compatible (need 2+)")
        elseif current_type !== nothing
            ig.TextDisabled("Select multiple measurements for combined plots")
        elseif length(selected_measurements) > 1
            ig.TextColored((0.8, 0.8, 0.2, 1.0), "Select a plot type above to continue")
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
                _open_combined_plot_window!(ui_state, proj, compatible_measurements, current_type)
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

        plot_kind = is_combined ? entry[:combined_kind] : entry[:measurement_kind]
        title = get(entry, :title, is_combined ? "Combined Plot" : basename(filepath))
        target_id = string(get(entry, :target_id, get(entry, :id, filepath)))
        entry[:target_id] = target_id

        cache_version = _cache_plot_version(ui_state)
        existing_cache_version = get(entry, :cache_version, nothing)
        existing_debug = get(entry, :debug, false)
        has_failure = haskey(entry, :plot_error)
        refresh = (!haskey(entry, :figure) && !has_failure) ||
            existing_cache_version != cache_version ||
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
