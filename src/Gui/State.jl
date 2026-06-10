function _project_for_preference(pref::AbstractString)
    pref == "auto" && return something(DEFAULT_PROJECT[])
    for project in PROJECTS
        project_name(project) == pref && return project
    end
    error("Unknown project preference '$pref'")
end

"""
Replace the current workspace with the project selected for one source root.
"""
function _open_project_path!(ui_state, path::String; persist=true)
    norm_path = _normalize_project_path(path)
    ui_state[:project_preference] = _project_preference_for_path(ui_state, norm_path)
    proj = _project_for_preference(ui_state[:project_preference])
    previous_workspace = get(ui_state, :workspace, nothing)
    previous_root = previous_workspace isa Workspace.Workspace ?
        previous_workspace.root_path :
        ""
    _cancel_figure_script_job!(ui_state)
    previous_workspace isa Workspace.Workspace && close_workspace!(previous_workspace)
    _clear_plot_views!(ui_state)
    delete!(ui_state, :main_plot_kind)
    delete!(ui_state, :main_plot_measurement_kind)
    ui_state[:main_plot_measurement_ids] = String[]
    ui_state[:main_plot_measurements] = MeasurementInfo[]
    ui_state[:open_plot_windows] = Vector{Dict{Symbol,Any}}()
    previous_root == norm_path ||
        _reset_figure_script_state!(ui_state, norm_path)
    _load_project_view!(ui_state, norm_path, proj)
    _load_tag_state_for_root!(ui_state, norm_path)
    ui_state[:workspace] = open_workspace(proj, norm_path)
    _apply_project_view!(ui_state, ui_state[:project_view_loaded])
    _invalidate_figure_script_scan_cache!(ui_state)
    persist && _persist_preferences!(ui_state; path=norm_path)
    return nothing
end

function _project_status_text(ui_state)
    workspace = get(ui_state, :workspace, nothing)
    workspace isa Workspace.Workspace || return "No project loaded"
    return "Active: $(project_name(workspace.project))"
end

function _init_tag_state!(ui_state)
    ui_state[:show_bad] = true
    ui_state[:tag_state] = nothing
    ui_state[:tag_state_error] = ""
    ui_state[:expanded_device_paths] = String[]
end

function _init_plot_state!(ui_state)
    ui_state[:plot_error] = ""
    ui_state[:plot_runtime_warmed] = false
    ui_state[:main_plot_live] = true
    ui_state[:main_plot_measurements] = MeasurementInfo[]
    ui_state[:_plot_window_counter] = 0
    ui_state[:plot_kind_by_measurement_kind] = Dict{String,String}()
end

"""Select and reveal every measurement produced by one source file."""
function select_source_file!(
    ui_state::Dict{Symbol,Any},
    filepath::AbstractString,
)::Bool
    workspace = get(ui_state, :workspace, nothing)
    workspace isa Workspace.Workspace || return false
    path = String(filepath)
    measurements = [
        measurement
        for measurement in workspace.index.hierarchy.all_measurements
        if measurement.filepath == path
    ]
    isempty(measurements) && return false

    device_paths =
        unique([device_path_key(measurement.device_info) for measurement in measurements])
    expanded_paths = copy(get(ui_state, :expanded_device_paths, String[]))
    for device_path in device_paths
        parts = split(device_path, '/')
        for depth in 1:(length(parts) - 1)
            parent_path = join(parts[1:depth], '/')
            parent_path in expanded_paths || push!(expanded_paths, parent_path)
        end
    end

    ui_state[:expanded_device_paths] = expanded_paths
    workspace.selection.device_paths = device_paths
    workspace.selection.measurement_ids =
        [measurement.unique_id for measurement in measurements]
    ui_state[:scroll_to_device_path] = first(device_paths)
    ui_state[:scroll_to_measurement_id] = first(measurements).unique_id
    return true
end

const _FIGURE_SCRIPT_NAME_BUFFER_SIZE = 192
const _FIGURE_GROUP_NAME_BUFFER_SIZE = 192
const _FIGURE_OUTPUT_DIR_BUFFER_SIZE = 512

function _text_buffer(capacity::Int)
    capacity > 1 || error("Text buffer capacity must be greater than 1")
    return fill(UInt8(0), capacity)
end

function _buffer_string(buffer::Vector{UInt8})
    terminator = findfirst(==(0x00), buffer)
    end_index = something(terminator, length(buffer) + 1) - 1
    return String(buffer[1:end_index])
end

function _set_buffer_string!(buffer::Vector{UInt8}, value::AbstractString)
    bytes = collect(codeunits(String(value)))
    length(bytes) < length(buffer) || error("String value exceeds text buffer capacity")
    fill!(buffer, 0x00)
    buffer[1:length(bytes)] = bytes
    return buffer
end

function _init_figure_script_state!(ui_state)
    ui_state[:show_figure_script_window] = false
    ui_state[:figure_script_root_path] = ""
    ui_state[:figure_script_output_dir_buffer] = _text_buffer(_FIGURE_OUTPUT_DIR_BUFFER_SIZE)
    ui_state[:figure_script_name_buffer] = _text_buffer(_FIGURE_SCRIPT_NAME_BUFFER_SIZE)
    ui_state[:figure_script_group_name_buffer] = _text_buffer(_FIGURE_GROUP_NAME_BUFFER_SIZE)
    ui_state[:figure_script_groups] = NamedMeasurementGroup[]
    ui_state[:figure_script_group_matches] = Dict{String,Vector{MeasurementInfo}}()
    ui_state[:figure_script_group_matches_valid] = false
    ui_state[:figure_script_fact_index] = nothing
    ui_state[:figure_script_fact_index_valid] = false
    ui_state[:figure_script_selected_group] = 0
    ui_state[:figure_script_error] = ""
    ui_state[:figure_script_status] = ""
    ui_state[:figure_script_overwrite_confirm] = ""
    ui_state[:figure_script_job_state] = :idle
    ui_state[:figure_script_job_seq] = 0
    ui_state[:active_figure_script_job_id] = 0
    ui_state[:figure_script_job_events] = nothing
    ui_state[:figure_script_job_profiles] = Any[]
end

function _remember_background_task!(ui_state, task::Task)::Task
    workspace = ui_state[:workspace]::Workspace.Workspace
    return track_task!(workspace, task)
end

"""
Stop browser and workspace work before the render loop exits.
"""
function _shutdown_background_jobs!(ui_state)::Nothing
    get(ui_state, :shutdown_complete, false) && return nothing
    ui_state[:shutting_down] = true
    _cancel_figure_script_job!(ui_state)
    workspace = get(ui_state, :workspace, nothing)
    workspace isa Workspace.Workspace && close_workspace!(workspace)
    ui_state[:shutdown_complete] = true
    return nothing
end

function _figure_script_groups(ui_state)
    return get(ui_state, :figure_script_groups, NamedMeasurementGroup[])
end

function _selected_figure_script_group(ui_state)
    index = get(ui_state, :figure_script_selected_group, 0)
    groups = _figure_script_groups(ui_state)
    1 <= index <= length(groups) || return nothing
    return groups[index]
end

function _set_selected_figure_script_group!(ui_state, index::Int)
    groups = _figure_script_groups(ui_state)
    if 1 <= index <= length(groups)
        ui_state[:figure_script_selected_group] = index
        _set_buffer_string!(ui_state[:figure_script_group_name_buffer], groups[index].name)
        return
    end
    ui_state[:figure_script_selected_group] = 0
    _set_buffer_string!(ui_state[:figure_script_group_name_buffer], "")
end

function _clear_figure_script_messages!(ui_state)
    ui_state[:figure_script_error] = ""
    ui_state[:figure_script_status] = ""
end

function _set_figure_script_error!(ui_state, err::Exception)
    ui_state[:figure_script_error] = sprint(showerror, err)
    ui_state[:figure_script_status] = ""
end

function _set_figure_script_status!(ui_state, message::AbstractString)
    ui_state[:figure_script_status] = String(message)
    ui_state[:figure_script_error] = ""
end

function _append_perf_sample!(ui_state, key::Symbol, duration_ms::Float64, alloc_bytes::Int)
    timings = get!(() -> Dict{Symbol,Vector{Float64}}(), ui_state, :_timings)
    allocs = get!(() -> Dict{Symbol,Vector{Int}}(), ui_state, :_allocs)
    vec = get!(() -> Float64[], timings, key)
    push!(vec, duration_ms)
    length(vec) > 400 && popfirst!(vec)
    avec = get!(() -> Int[], allocs, key)
    push!(avec, alloc_bytes)
    length(avec) > 400 && popfirst!(avec)
    return nothing
end

function _record_figure_script_job_profile!(ui_state, operation::Symbol, profile::FigureScriptInferenceProfile)
    history = get!(ui_state, :figure_script_job_profiles) do
        Any[]
    end
    entry = (
        operation=operation,
        profile=profile,
    )
    push!(history, entry)
    length(history) > 24 && popfirst!(history)

    total_alloc = sum(section.alloc_bytes for section in profile.sections)
    _append_perf_sample!(ui_state, :figure_script_job_total, profile.total_ms, total_alloc)
    for section in profile.sections
        key = Symbol("figure_script_" * String(section.key))
        _append_perf_sample!(ui_state, key, section.duration_ms, section.alloc_bytes)
    end
    return nothing
end

function _figure_script_job_running(ui_state)
    return get(ui_state, :figure_script_job_state, :idle) == :running
end

function _finalize_figure_script_job!(ui_state)
    ui_state[:figure_script_job_state] = :idle
    ui_state[:figure_script_job_events] = nothing
end

function _cancel_figure_script_job!(ui_state)
    _figure_script_job_running(ui_state) || return
    ui_state[:active_figure_script_job_id] = get(ui_state, :active_figure_script_job_id, 0) + 1
    _clear_figure_script_messages!(ui_state)
    _finalize_figure_script_job!(ui_state)
end

function _figure_script_job_request(
    ui_state,
    operation::Symbol,
    groups::Vector{NamedMeasurementGroup},
    selected_measurements::Vector{MeasurementInfo},
    all_measurements::Vector{MeasurementInfo},
    group_name::AbstractString;
    group_index::Int=0,
    status_message::AbstractString,
    completion_message::AbstractString,
)
    operation in (:create, :replace) || error("Unsupported figure-script job operation '$operation'")
    workspace = ui_state[:workspace]::Workspace.Workspace
    return Dict{Symbol,Any}(
        :operation => operation,
        :groups => copy(groups),
        :selected_measurements => copy(selected_measurements),
        :all_measurements => all_measurements,
        :group_name => String(group_name),
        :group_index => group_index,
        :scan_id => workspace.scan.id,
        :root_path => workspace.root_path,
        :tag_state => get(ui_state, :tag_state, nothing),
        :status_message => String(status_message),
        :completion_message => String(completion_message),
    )
end

function _resolve_figure_script_job(request::Dict{Symbol,Any})
    inferred_group, profile = _infer_measurement_group_profiled(
        request[:group_name],
        request[:selected_measurements],
        request[:all_measurements];
        tag_state=get(request, :tag_state, nothing),
    )
    groups = copy(request[:groups])
    if request[:operation] == :create
        push!(groups, inferred_group)
        selected_index = length(groups)
    else
        selected_index = request[:group_index]
        1 <= selected_index <= length(groups) || error("Figure-script group index $selected_index is out of bounds")
        groups[selected_index] = inferred_group
    end
    _validate_named_measurement_groups(groups)
    return groups, selected_index, profile
end

function _launch_figure_script_job!(ui_state, request::Dict{Symbol,Any})
    _figure_script_job_running(ui_state) && throw(FigureScriptValidationError(
        "Wait for the current figure-script update to finish",
    ))

    job_id = get(ui_state, :figure_script_job_seq, 0) + 1
    ui_state[:figure_script_job_seq] = job_id
    ui_state[:active_figure_script_job_id] = job_id
    ui_state[:figure_script_job_state] = :running
    _set_figure_script_status!(ui_state, request[:status_message])

    events = Channel{NamedTuple}(1)
    ui_state[:figure_script_job_events] = events

    task = Base.Threads.@spawn begin
        try
            groups, selected_index, profile = _resolve_figure_script_job(request)
            put!(events, (
                kind=:done,
                job_id=job_id,
                operation=request[:operation],
                group_name=request[:group_name],
                profile=profile,
                groups=groups,
                selected_index=selected_index,
                scan_id=request[:scan_id],
                root_path=request[:root_path],
                completion_message=request[:completion_message],
            ))
        catch err
            profile = nothing
            cause = err
            bt = catch_backtrace()
            if err isa FigureScriptProfiledError
                profile = err.profile
                cause = err.cause
                bt = err.bt
            end
            put!(events, (kind=:error, job_id=job_id, operation=request[:operation], profile=profile, error=cause, bt=bt))
        finally
            close(events)
        end
    end
    _remember_background_task!(ui_state, task)
end

function _poll_figure_script_job_events!(ui_state)
    events = get(ui_state, :figure_script_job_events, nothing)
    events === nothing && return

    while isready(events)
        msg = take!(events)
        msg.job_id == get(ui_state, :active_figure_script_job_id, 0) || continue

        if msg.kind == :done
            _record_figure_script_job_profile!(ui_state, msg.operation, msg.profile)
            workspace = ui_state[:workspace]::Workspace.Workspace
            if msg.scan_id == workspace.scan.id &&
               msg.root_path == workspace.root_path
                ui_state[:figure_script_groups] = msg.groups
                _invalidate_figure_script_group_matches!(ui_state)
                _set_selected_figure_script_group!(ui_state, msg.selected_index)
                _set_figure_script_status!(
                    ui_state,
                    "$(msg.completion_message) ($(round(msg.profile.total_ms / 1000; digits=2)) s)",
                )
            end
            _finalize_figure_script_job!(ui_state)
        elseif msg.kind == :error
            msg.profile === nothing || _record_figure_script_job_profile!(ui_state, msg.operation, msg.profile)
            @error "Figure-script job failed" exception = (msg.error, msg.bt)
            _set_figure_script_error!(ui_state, msg.error)
            _finalize_figure_script_job!(ui_state)
        end
    end
end

function _reset_figure_script_state!(ui_state, root_path::AbstractString="")
    _init_figure_script_state!(ui_state)
    normalized_root = String(root_path)
    ui_state[:figure_script_root_path] = normalized_root
    if !isempty(normalized_root)
        _set_buffer_string!(
            ui_state[:figure_script_output_dir_buffer],
            _figure_script_output_dir_for_path(ui_state, normalized_root),
        )
    end
end

function _current_scan_measurements(ui_state)
    workspace = get(ui_state, :workspace, nothing)
    workspace isa Workspace.Workspace || return MeasurementInfo[]
    return workspace.index.hierarchy.all_measurements
end

function _invalidate_figure_script_group_matches!(ui_state)
    ui_state[:figure_script_group_matches] = Dict{String,Vector{MeasurementInfo}}()
    ui_state[:figure_script_group_matches_valid] = false
    return nothing
end

function _invalidate_figure_script_scan_cache!(ui_state)
    _invalidate_figure_script_group_matches!(ui_state)
    ui_state[:figure_script_fact_index] = nothing
    ui_state[:figure_script_fact_index_valid] = false
    return nothing
end

function _ensure_figure_script_fact_index!(ui_state)
    if get(ui_state, :figure_script_fact_index_valid, false)
        index = get(ui_state, :figure_script_fact_index, nothing)
        index isa _FigureScriptFactIndex && return index
    end

    index = _build_figure_script_fact_index(_current_scan_measurements(ui_state))
    ui_state[:figure_script_fact_index] = index
    ui_state[:figure_script_fact_index_valid] = true
    return index
end

function _refresh_figure_script_group_matches!(ui_state)
    groups = _figure_script_groups(ui_state)
    index = _ensure_figure_script_fact_index!(ui_state)
    matches = Dict{String,Vector{MeasurementInfo}}()
    for group in groups
        matches[group.name] = _matching_measurements(index, group)
    end
    ui_state[:figure_script_group_matches] = matches
    ui_state[:figure_script_group_matches_valid] = true
    return matches
end

function _ensure_figure_script_group_matches!(ui_state)
    if !get(ui_state, :figure_script_group_matches_valid, false)
        return _refresh_figure_script_group_matches!(ui_state)
    end
    return get!(ui_state, :figure_script_group_matches) do
        Dict{String,Vector{MeasurementInfo}}()
    end
end

function _selected_measurements_in_panel_order(ui_state)
    filter_meas = get(ui_state, :_imgui_text_filter_meas, nothing)
    all_measurements = _measurements_of_selected_devices(ui_state)
    visible_measurements = if filter_meas === nothing
        all_measurements
    else
        workspace = ui_state[:workspace]::Workspace.Workspace
        _visible_measurements(ui_state, workspace.project, all_measurements, filter_meas)
    end
    workspace = ui_state[:workspace]::Workspace.Workspace
    selected_ids = Set(workspace.selection.measurement_ids)
    return [measurement for measurement in visible_measurements if measurement.unique_id in selected_ids]
end

function _group_measurements_in_current_scan(ui_state, group::NamedMeasurementGroup)
    matches = _ensure_figure_script_group_matches!(ui_state)
    return get(matches, group.name, MeasurementInfo[])
end

function _create_figure_script_group_from_selection!(ui_state)
    selected_measurements = _selected_measurements_in_panel_order(ui_state)
    isempty(selected_measurements) && throw(FigureScriptValidationError("Select one or more measurements before creating a group"))
    name = strip(_buffer_string(ui_state[:figure_script_group_name_buffer]))
    isempty(name) && throw(FigureScriptValidationError("Enter a group name before creating a group"))

    request = _figure_script_job_request(
        ui_state,
        :create,
        _figure_script_groups(ui_state),
        selected_measurements,
        _current_scan_measurements(ui_state),
        name;
        status_message="Creating group '$name'...",
        completion_message="Created group '$name'",
    )
    _launch_figure_script_job!(ui_state, request)
    return nothing
end

function _add_selection_to_figure_script_group!(ui_state)
    group = _selected_figure_script_group(ui_state)
    group === nothing && throw(FigureScriptValidationError("Select a group before adding measurements"))
    selected_measurements = _selected_measurements_in_panel_order(ui_state)
    isempty(selected_measurements) && throw(FigureScriptValidationError("Select one or more measurements before adding them"))

    all_measurements = _current_scan_measurements(ui_state)
    merged_ids = Set(measurement.unique_id for measurement in _group_measurements_in_current_scan(ui_state, group))
    foreach(measurement -> push!(merged_ids, measurement.unique_id), selected_measurements)
    merged_measurements = [measurement for measurement in all_measurements if measurement.unique_id in merged_ids]
    request = _figure_script_job_request(
        ui_state,
        :replace,
        _figure_script_groups(ui_state),
        merged_measurements,
        all_measurements,
        group.name;
        group_index=get(ui_state, :figure_script_selected_group, 0),
        status_message="Updating group '$(group.name)'...",
        completion_message="Updated group '$(group.name)'",
    )
    _launch_figure_script_job!(ui_state, request)
    return nothing
end

function _remove_selection_from_figure_script_group!(ui_state)
    group = _selected_figure_script_group(ui_state)
    group === nothing && throw(FigureScriptValidationError("Select a group before removing measurements"))
    selected_measurements = _selected_measurements_in_panel_order(ui_state)
    isempty(selected_measurements) && throw(FigureScriptValidationError("Select one or more measurements before removing them"))

    remaining_ids = Set(measurement.unique_id for measurement in _group_measurements_in_current_scan(ui_state, group))
    foreach(measurement -> delete!(remaining_ids, measurement.unique_id), selected_measurements)
    isempty(remaining_ids) && throw(FigureScriptValidationError("Measurement groups cannot be empty"))
    all_measurements = _current_scan_measurements(ui_state)
    remaining_measurements = [measurement for measurement in all_measurements if measurement.unique_id in remaining_ids]
    request = _figure_script_job_request(
        ui_state,
        :replace,
        _figure_script_groups(ui_state),
        remaining_measurements,
        all_measurements,
        group.name;
        group_index=get(ui_state, :figure_script_selected_group, 0),
        status_message="Updating group '$(group.name)'...",
        completion_message="Updated group '$(group.name)'",
    )
    _launch_figure_script_job!(ui_state, request)
    return nothing
end

function _rename_selected_figure_script_group!(ui_state)
    group = _selected_figure_script_group(ui_state)
    group === nothing && throw(FigureScriptValidationError("Select a group before renaming it"))
    name = strip(_buffer_string(ui_state[:figure_script_group_name_buffer]))
    isempty(name) && throw(FigureScriptValidationError("Group name cannot be empty"))

    groups = copy(_figure_script_groups(ui_state))
    selected_index = get(ui_state, :figure_script_selected_group, 0)
    groups[selected_index] = NamedMeasurementGroup(name, group.filter)
    _validate_named_measurement_groups(groups)
    ui_state[:figure_script_groups] = groups
    _invalidate_figure_script_group_matches!(ui_state)
    _set_selected_figure_script_group!(ui_state, selected_index)
    return nothing
end

function _delete_selected_figure_script_group!(ui_state)
    group = _selected_figure_script_group(ui_state)
    group === nothing && throw(FigureScriptValidationError("Select a group before deleting it"))
    groups = copy(_figure_script_groups(ui_state))
    deleteat!(groups, get(ui_state, :figure_script_selected_group, 0))
    ui_state[:figure_script_groups] = groups
    _invalidate_figure_script_group_matches!(ui_state)
    _set_selected_figure_script_group!(ui_state, min(get(ui_state, :figure_script_selected_group, 0), length(groups)))
    return nothing
end

# Timing & allocation utilities
# usage: _time!(ui_state, :key) do ... end
function _time!(f::Function, ui_state, key::Symbol)
    timings = get!(() -> Dict{Symbol,Vector{Float64}}(), ui_state, :_timings)
    allocs = get!(() -> Dict{Symbol,Vector{Int}}(), ui_state, :_allocs)
    t0 = time_ns()
    bytes = @allocated f()
    dt_ms = (time_ns() - t0) / 1e6
    vec = get!(() -> Float64[], timings, key)
    push!(vec, dt_ms)
    length(vec) > 400 && popfirst!(vec)
    avec = get!(() -> Int[], allocs, key)
    push!(avec, bytes)
    length(avec) > 400 && popfirst!(avec)
    nothing
end

function _read_proc_int(path::String, prefix::String)
    isfile(path) || return nothing
    for line in eachline(path)
        startswith(line, prefix) || continue
        fields = split(strip(line))
        length(fields) >= 2 || return nothing
        try
            return parse(Int, fields[2])
        catch
            return nothing
        end
    end
    return nothing
end

function _memory_snapshot()
    return (
        vmrss_kb=_read_proc_int("/proc/self/status", "VmRSS:"),
        rssanon_kb=_read_proc_int("/proc/self/status", "RssAnon:"),
        vmsize_kb=_read_proc_int("/proc/self/status", "VmSize:"),
        vmpeak_kb=_read_proc_int("/proc/self/status", "VmPeak:"),
        read_bytes=_read_proc_int("/proc/self/io", "read_bytes:"),
        gc_live_bytes=Int(Base.gc_live_bytes()),
        maxrss_bytes=Int(Sys.maxrss()),
    )
end

function _fmt_kb(kb::Union{Nothing,Integer})
    kb === nothing && return "n/a"
    gib = kb / (1024^2)
    gib >= 1 && return @sprintf("%.2f GiB", gib)
    return @sprintf("%.0f MiB", kb / 1024)
end

function _fmt_bytes(bytes::Union{Nothing,Integer})
    bytes === nothing && return "n/a"
    gib = bytes / (1024^3)
    gib >= 1 && return @sprintf("%.2f GiB", gib)
    return @sprintf("%.0f MiB", bytes / (1024^2))
end

function _collect_gl_info!()
    try
        Dict(
            :vendor => unsafe_string(gl.glGetString(gl.GL_VENDOR)),
            :renderer => unsafe_string(gl.glGetString(gl.GL_RENDERER)),
            :version => unsafe_string(gl.glGetString(gl.GL_VERSION)),
            :sl => unsafe_string(gl.glGetString(gl.GL_SHADING_LANGUAGE_VERSION)),
        )
    catch err
        @warn "GL info query failed" error = err
        Dict{Symbol,String}()
    end
end

function _print_perf_summary(ui_state)
    @debug begin
        gi = ui_state[:_gl_info]
        timings = get(ui_state, :_timings, Dict{Symbol,Vector{Float64}}())
        allocs = get(ui_state, :_allocs, Dict{Symbol,Vector{Int}}())
        msg = """\n
        ==== Performance Summary ====
        GL Vendor:   $(get(gi, :vendor, "?"))
        GL Renderer: $(get(gi, :renderer, "?"))
        GL Version:  $(get(gi, :version, "?"))
        """
        if !isempty(timings)
            msg = msg * @sprintf "%-12s %5s %9s %9s %9s %12s %12s\n" "Key" "n" "Mean(ms)" "Max(ms)" "Last(ms)" "AllocMean(KB)" "AllocLast(KB)"
        end
        for k in sort(collect(keys(timings)))
            v = timings[k]
            isempty(v) && continue
            a = get(allocs, k, Int[])
            n = length(v)
            mean_ms = mean(v)
            max_ms = maximum(v)
            last_ms = v[end]
            mean_alloc = isempty(a) ? 0.0 : mean(a) / 1024
            last_alloc = isempty(a) ? 0.0 : a[end] / 1024
            msg = msg * @sprintf "%-12s %5d %9.2f %9.2f %9.2f %12.1f %12.1f\n" String(k) n mean_ms max_ms last_ms mean_alloc last_alloc
        end
        msg * "=============================="
    end
end

function _helpmarker(desc::String)
    ig.TextDisabled("(?)")
    if ig.BeginItemTooltip()
        ig.PushTextWrapPos(ig.GetFontSize() * 35.0)
        ig.TextUnformatted(desc)
        ig.PopTextWrapPos()
        ig.EndTooltip()
    end
end
