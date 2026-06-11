"""Return the project selected by the saved project preference."""
function _project_for_preference(pref::AbstractString)::AbstractProject
    pref == "auto" && return something(DEFAULT_PROJECT[])
    for project in PROJECTS
        project_name(project) == pref && return project
    end
    error("Unknown project preference '$pref'")
end

"""
Replace the current workspace with the project selected for one source root.
"""
function _open_project_path!(
    state::BrowserState,
    path::String;
    persist::Bool=true,
)::Nothing
    norm_path = _normalize_project_path(path)
    state.project_preference = _project_preference_for_path(state, norm_path)
    project = _project_for_preference(state.project_preference)
    previous_workspace = state.workspace
    previous_root = previous_workspace isa Workspace.Workspace ?
        previous_workspace.root_path :
        ""
    _cancel_figure_script_job!(state)
    previous_workspace isa Workspace.Workspace && close_workspace!(previous_workspace)
    state.plots = PlotState(debug=state.plots.debug)
    previous_root == norm_path ||
        _reset_figure_script_state!(state, norm_path)
    view = _load_project_view(norm_path)
    !isempty(view.project) && view.project != project_name(project) &&
        (view = PersistedProjectView(project=project_name(project)))
    _load_tag_state_for_root!(state, norm_path)
    state.workspace = open_workspace(project, norm_path)
    _apply_project_view!(state, view)
    state.saved_project_view = view
    _invalidate_figure_script_scan_cache!(state)
    persist && _persist_preferences!(state; path=norm_path)
    return nothing
end

"""Select and reveal every measurement produced by one source file."""
function select_source_file!(
    state::BrowserState,
    filepath::AbstractString,
)::Bool
    workspace = state.workspace
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
    expanded_paths = copy(state.expanded_device_paths)
    for device_path in device_paths
        parts = split(device_path, '/')
        for depth in 1:(length(parts) - 1)
            parent_path = join(parts[1:depth], '/')
            parent_path in expanded_paths || push!(expanded_paths, parent_path)
        end
    end

    state.expanded_device_paths = expanded_paths
    workspace.selection.device_paths = device_paths
    workspace.selection.measurement_ids =
        [measurement.unique_id for measurement in measurements]
    state.scroll_to_device_path = first(device_paths)
    state.scroll_to_measurement_id = first(measurements).unique_id
    return true
end

"""Read a null-terminated CImGui text buffer."""
function _buffer_string(buffer::Vector{UInt8})::String
    terminator = findfirst(==(0x00), buffer)
    end_index = something(terminator, length(buffer) + 1) - 1
    return String(buffer[1:end_index])
end

"""Replace the text in a fixed-size CImGui input buffer."""
function _set_buffer_string!(
    buffer::Vector{UInt8},
    value::AbstractString,
)::Vector{UInt8}
    bytes = collect(codeunits(String(value)))
    length(bytes) < length(buffer) || error("String value exceeds text buffer capacity")
    fill!(buffer, 0x00)
    buffer[1:length(bytes)] = bytes
    return buffer
end

"""
Stop browser and workspace work before the render loop exits.
"""
function _shutdown_background_jobs!(state::BrowserState)::Nothing
    state.shutdown_complete && return nothing
    state.shutting_down = true
    _cancel_figure_script_job!(state)
    workspace = state.workspace
    workspace isa Workspace.Workspace && close_workspace!(workspace)
    state.shutdown_complete = true
    return nothing
end

function _selected_figure_script_group(
    state::BrowserState,
)::Union{Nothing,NamedMeasurementGroup}
    index = state.figure_scripts.selected_group
    groups = state.figure_scripts.groups
    1 <= index <= length(groups) || return nothing
    return groups[index]
end

"""Select one figure-script group and copy its name into the editor."""
function _set_selected_figure_script_group!(
    state::BrowserState,
    index::Int,
)::Nothing
    figure_scripts = state.figure_scripts
    groups = figure_scripts.groups
    if 1 <= index <= length(groups)
        figure_scripts.selected_group = index
        _set_buffer_string!(figure_scripts.group_name_buffer, groups[index].name)
        return nothing
    end
    figure_scripts.selected_group = 0
    _set_buffer_string!(figure_scripts.group_name_buffer, "")
    return nothing
end

"""Clear the user-facing figure-script result messages."""
function _clear_figure_script_messages!(state::BrowserState)::Nothing
    state.figure_scripts.error = ""
    state.figure_scripts.status = ""
    return nothing
end

"""Show a figure-script error and clear any previous success message."""
function _set_figure_script_error!(
    state::BrowserState,
    err::Exception,
)::Nothing
    state.figure_scripts.error = sprint(showerror, err)
    state.figure_scripts.status = ""
    return nothing
end

"""Show a figure-script success message and clear any previous error."""
function _set_figure_script_status!(
    state::BrowserState,
    message::AbstractString,
)::Nothing
    state.figure_scripts.status = String(message)
    state.figure_scripts.error = ""
    return nothing
end

"""Retain a bounded history of one performance measurement."""
function _append_perf_sample!(
    state::BrowserState,
    key::Symbol,
    duration_ms::Float64,
    alloc_bytes::Int,
)::Nothing
    samples = get!(() -> Float64[], state.performance.timings, key)
    allocations = get!(() -> Int[], state.performance.allocations, key)
    push!(samples, duration_ms)
    length(samples) > 400 && popfirst!(samples)
    push!(allocations, alloc_bytes)
    length(allocations) > 400 && popfirst!(allocations)
    return nothing
end

"""Add one completed figure-script profile to diagnostics."""
function _record_figure_script_job_profile!(
    state::BrowserState,
    operation::Symbol,
    profile::FigureScriptInferenceProfile,
)::Nothing
    history = state.figure_scripts.job_profiles
    entry = (
        operation=operation,
        profile=profile,
    )
    push!(history, entry)
    length(history) > 24 && popfirst!(history)

    total_alloc = sum(section.alloc_bytes for section in profile.sections)
    _append_perf_sample!(
        state,
        :figure_script_job_total,
        profile.total_ms,
        total_alloc,
    )
    for section in profile.sections
        key = Symbol("figure_script_" * String(section.key))
        _append_perf_sample!(state, key, section.duration_ms, section.alloc_bytes)
    end
    return nothing
end

function _figure_script_job_running(state::BrowserState)::Bool
    return state.figure_scripts.job_state == :running
end

"""Forget the completed or canceled figure-script worker."""
function _finalize_figure_script_job!(state::BrowserState)::Nothing
    state.figure_scripts.job_state = :idle
    state.figure_scripts.job_events = nothing
    return nothing
end

"""Invalidate the active figure-script worker so late results are ignored."""
function _cancel_figure_script_job!(state::BrowserState)::Nothing
    _figure_script_job_running(state) || return nothing
    state.figure_scripts.active_job_id += 1
    _clear_figure_script_messages!(state)
    _finalize_figure_script_job!(state)
    return nothing
end

"""
Capture immutable inputs for one figure-script group update.

The worker receives this copy so later browser edits cannot change in-flight work.
"""
function _figure_script_job_request(
    state::BrowserState,
    operation::Symbol,
    groups::Vector{NamedMeasurementGroup},
    selected_measurements::Vector{MeasurementInfo},
    all_measurements::Vector{MeasurementInfo},
    group_name::AbstractString;
    group_index::Int=0,
    status_message::AbstractString,
    completion_message::AbstractString,
)::Dict{Symbol,Any}
    operation in (:create, :replace) || error("Unsupported figure-script job operation '$operation'")
    workspace = state.workspace::Workspace.Workspace
    return Dict{Symbol,Any}(
        :operation => operation,
        :groups => copy(groups),
        :selected_measurements => copy(selected_measurements),
        :all_measurements => all_measurements,
        :group_name => String(group_name),
        :group_index => group_index,
        :scan_id => workspace.scan.id,
        :root_path => workspace.root_path,
        :tag_state => state.tag_state,
        :status_message => String(status_message),
        :completion_message => String(completion_message),
    )
end

"""Infer and validate the groups produced by one figure-script request."""
function _resolve_figure_script_job(
    request::Dict{Symbol,Any},
)::Tuple{Vector{NamedMeasurementGroup},Int,FigureScriptInferenceProfile}
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

"""Run one figure-script group update outside the render loop."""
function _launch_figure_script_job!(
    state::BrowserState,
    request::Dict{Symbol,Any},
)::Nothing
    _figure_script_job_running(state) && throw(FigureScriptValidationError(
        "Wait for the current figure-script update to finish",
    ))

    figure_scripts = state.figure_scripts
    job_id = figure_scripts.job_sequence + 1
    figure_scripts.job_sequence = job_id
    figure_scripts.active_job_id = job_id
    figure_scripts.job_state = :running
    _set_figure_script_status!(state, request[:status_message])

    events = Channel{NamedTuple}(1)
    figure_scripts.job_events = events

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
    track_task!(state.workspace::Workspace.Workspace, task)
    return nothing
end

"""Apply completed figure-script worker messages to browser state."""
function _poll_figure_script_job_events!(state::BrowserState)::Nothing
    events = state.figure_scripts.job_events
    events === nothing && return nothing

    while isready(events)
        msg = take!(events)
        msg.job_id == state.figure_scripts.active_job_id || continue

        if msg.kind == :done
            _record_figure_script_job_profile!(state, msg.operation, msg.profile)
            workspace = state.workspace::Workspace.Workspace
            if msg.scan_id == workspace.scan.id &&
               msg.root_path == workspace.root_path
                state.figure_scripts.groups = msg.groups
                _invalidate_figure_script_group_matches!(state)
                _set_selected_figure_script_group!(state, msg.selected_index)
                _set_figure_script_status!(
                    state,
                    "$(msg.completion_message) ($(round(msg.profile.total_ms / 1000; digits=2)) s)",
                )
            end
            _finalize_figure_script_job!(state)
        elseif msg.kind == :error
            msg.profile === nothing ||
                _record_figure_script_job_profile!(state, msg.operation, msg.profile)
            @error "Figure-script job failed" exception = (msg.error, msg.bt)
            _set_figure_script_error!(state, msg.error)
            _finalize_figure_script_job!(state)
        end
    end
    return nothing
end

"""Clear figure-script state when the browser changes source roots."""
function _reset_figure_script_state!(
    state::BrowserState,
    root_path::AbstractString="",
)::Nothing
    state.figure_scripts = FigureScriptState()
    normalized_root = String(root_path)
    state.figure_scripts.root_path = normalized_root
    if !isempty(normalized_root)
        _set_buffer_string!(
            state.figure_scripts.output_dir_buffer,
            _figure_script_output_dir_for_path(state, normalized_root),
        )
    end
    return nothing
end

"""Return the measurements currently present in the workspace index."""
function _current_scan_measurements(state::BrowserState)::Vector{MeasurementInfo}
    workspace = state.workspace
    workspace isa Workspace.Workspace || return MeasurementInfo[]
    return workspace.index.hierarchy.all_measurements
end

"""Discard cached matches after group definitions change."""
function _invalidate_figure_script_group_matches!(state::BrowserState)::Nothing
    empty!(state.figure_scripts.group_matches)
    state.figure_scripts.group_matches_valid = false
    return nothing
end

"""Discard every figure-script value derived from the measurement index."""
function _invalidate_figure_script_scan_cache!(state::BrowserState)::Nothing
    _invalidate_figure_script_group_matches!(state)
    state.figure_scripts.fact_index = nothing
    state.figure_scripts.fact_index_valid = false
    return nothing
end

"""Build the figure-script search index once for the current measurement index."""
function _ensure_figure_script_fact_index!(
    state::BrowserState,
)::_FigureScriptFactIndex
    if state.figure_scripts.fact_index_valid
        index = state.figure_scripts.fact_index
        index isa _FigureScriptFactIndex && return index
    end

    index = _build_figure_script_fact_index(_current_scan_measurements(state))
    state.figure_scripts.fact_index = index
    state.figure_scripts.fact_index_valid = true
    return index
end

"""Recompute which measurements match every figure-script group."""
function _refresh_figure_script_group_matches!(
    state::BrowserState,
)::Dict{String,Vector{MeasurementInfo}}
    groups = state.figure_scripts.groups
    index = _ensure_figure_script_fact_index!(state)
    matches = Dict{String,Vector{MeasurementInfo}}()
    for group in groups
        matches[group.name] = _matching_measurements(index, group)
    end
    state.figure_scripts.group_matches = matches
    state.figure_scripts.group_matches_valid = true
    return matches
end

"""Return current figure-script group matches, computing them when needed."""
function _ensure_figure_script_group_matches!(
    state::BrowserState,
)::Dict{String,Vector{MeasurementInfo}}
    if !state.figure_scripts.group_matches_valid
        return _refresh_figure_script_group_matches!(state)
    end
    return state.figure_scripts.group_matches
end

"""Return selected measurements in the same order shown by the browser."""
function _selected_measurements_in_panel_order(
    state::BrowserState,
)::Vector{MeasurementInfo}
    filter_meas = state.measurement_filter_widget
    all_measurements = _measurements_of_selected_devices(state)
    visible_measurements = if filter_meas === nothing
        all_measurements
    else
        workspace = state.workspace::Workspace.Workspace
        _visible_measurements(state, workspace.project, all_measurements, filter_meas)
    end
    workspace = state.workspace::Workspace.Workspace
    selected_ids = Set(workspace.selection.measurement_ids)
    return [measurement for measurement in visible_measurements if measurement.unique_id in selected_ids]
end

"""Return current measurements matching one figure-script group."""
function _group_measurements_in_current_scan(
    state::BrowserState,
    group::NamedMeasurementGroup,
)::Vector{MeasurementInfo}
    matches = _ensure_figure_script_group_matches!(state)
    return get(matches, group.name, MeasurementInfo[])
end

"""Create a figure-script group from the current browser selection."""
function _create_figure_script_group_from_selection!(
    state::BrowserState,
)::Nothing
    selected_measurements = _selected_measurements_in_panel_order(state)
    isempty(selected_measurements) && throw(FigureScriptValidationError("Select one or more measurements before creating a group"))
    name = strip(_buffer_string(state.figure_scripts.group_name_buffer))
    isempty(name) && throw(FigureScriptValidationError("Enter a group name before creating a group"))

    request = _figure_script_job_request(
        state,
        :create,
        state.figure_scripts.groups,
        selected_measurements,
        _current_scan_measurements(state),
        name;
        status_message="Creating group '$name'...",
        completion_message="Created group '$name'",
    )
    _launch_figure_script_job!(state, request)
    return nothing
end

"""Add the current browser selection to the selected figure-script group."""
function _add_selection_to_figure_script_group!(state::BrowserState)::Nothing
    group = _selected_figure_script_group(state)
    group === nothing && throw(FigureScriptValidationError("Select a group before adding measurements"))
    selected_measurements = _selected_measurements_in_panel_order(state)
    isempty(selected_measurements) && throw(FigureScriptValidationError("Select one or more measurements before adding them"))

    all_measurements = _current_scan_measurements(state)
    merged_ids = Set(
        measurement.unique_id
        for measurement in _group_measurements_in_current_scan(state, group)
    )
    foreach(measurement -> push!(merged_ids, measurement.unique_id), selected_measurements)
    merged_measurements = [measurement for measurement in all_measurements if measurement.unique_id in merged_ids]
    request = _figure_script_job_request(
        state,
        :replace,
        state.figure_scripts.groups,
        merged_measurements,
        all_measurements,
        group.name;
        group_index=state.figure_scripts.selected_group,
        status_message="Updating group '$(group.name)'...",
        completion_message="Updated group '$(group.name)'",
    )
    _launch_figure_script_job!(state, request)
    return nothing
end

"""Remove the current browser selection from the selected figure-script group."""
function _remove_selection_from_figure_script_group!(state::BrowserState)::Nothing
    group = _selected_figure_script_group(state)
    group === nothing && throw(FigureScriptValidationError("Select a group before removing measurements"))
    selected_measurements = _selected_measurements_in_panel_order(state)
    isempty(selected_measurements) && throw(FigureScriptValidationError("Select one or more measurements before removing them"))

    remaining_ids = Set(
        measurement.unique_id
        for measurement in _group_measurements_in_current_scan(state, group)
    )
    foreach(measurement -> delete!(remaining_ids, measurement.unique_id), selected_measurements)
    isempty(remaining_ids) && throw(FigureScriptValidationError("Measurement groups cannot be empty"))
    all_measurements = _current_scan_measurements(state)
    remaining_measurements = [measurement for measurement in all_measurements if measurement.unique_id in remaining_ids]
    request = _figure_script_job_request(
        state,
        :replace,
        state.figure_scripts.groups,
        remaining_measurements,
        all_measurements,
        group.name;
        group_index=state.figure_scripts.selected_group,
        status_message="Updating group '$(group.name)'...",
        completion_message="Updated group '$(group.name)'",
    )
    _launch_figure_script_job!(state, request)
    return nothing
end

"""Rename the selected figure-script group after validating all group names."""
function _rename_selected_figure_script_group!(state::BrowserState)::Nothing
    group = _selected_figure_script_group(state)
    group === nothing && throw(FigureScriptValidationError("Select a group before renaming it"))
    name = strip(_buffer_string(state.figure_scripts.group_name_buffer))
    isempty(name) && throw(FigureScriptValidationError("Group name cannot be empty"))

    groups = copy(state.figure_scripts.groups)
    selected_index = state.figure_scripts.selected_group
    groups[selected_index] = NamedMeasurementGroup(name, group.filter)
    _validate_named_measurement_groups(groups)
    state.figure_scripts.groups = groups
    _invalidate_figure_script_group_matches!(state)
    _set_selected_figure_script_group!(state, selected_index)
    return nothing
end

"""Delete the selected figure-script group and select its nearest neighbor."""
function _delete_selected_figure_script_group!(state::BrowserState)::Nothing
    group = _selected_figure_script_group(state)
    group === nothing && throw(FigureScriptValidationError("Select a group before deleting it"))
    groups = copy(state.figure_scripts.groups)
    deleteat!(groups, state.figure_scripts.selected_group)
    state.figure_scripts.groups = groups
    _invalidate_figure_script_group_matches!(state)
    _set_selected_figure_script_group!(
        state,
        min(state.figure_scripts.selected_group, length(groups)),
    )
    return nothing
end

"""Measure one browser operation and retain its latest timing and allocation samples."""
function _time!(
    f::Function,
    state::BrowserState,
    key::Symbol,
)::Nothing
    t0 = time_ns()
    bytes = @allocated f()
    dt_ms = (time_ns() - t0) / 1e6
    _append_perf_sample!(state, key, dt_ms, bytes)
    return nothing
end

"""Read one integer field from a Linux process-information file."""
function _read_proc_int(
    path::String,
    prefix::String,
)::Union{Nothing,Int}
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

"""Collect the process-memory values shown in the performance window."""
function _memory_snapshot()::NamedTuple
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

"""Format a kibibyte count for the performance window."""
function _fmt_kb(kb::Union{Nothing,Integer})::String
    kb === nothing && return "n/a"
    gib = kb / (1024^2)
    gib >= 1 && return @sprintf("%.2f GiB", gib)
    return @sprintf("%.0f MiB", kb / 1024)
end

"""Format a byte count for the performance window."""
function _fmt_bytes(bytes::Union{Nothing,Integer})::String
    bytes === nothing && return "n/a"
    gib = bytes / (1024^3)
    gib >= 1 && return @sprintf("%.2f GiB", gib)
    return @sprintf("%.0f MiB", bytes / (1024^2))
end

"""Read identifying strings from the active OpenGL context."""
function _gl_info()::Dict{Symbol,String}
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

"""Write the collected performance samples to the debug log at shutdown."""
function _print_perf_summary(state::BrowserState)::Nothing
    @debug begin
        gi = state.performance.gl_info
        timings = state.performance.timings
        allocs = state.performance.allocations
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
    return nothing
end

"""Render a small question mark with explanatory hover text."""
function _helpmarker(desc::String)::Nothing
    ig.TextDisabled("(?)")
    if ig.BeginItemTooltip()
        ig.PushTextWrapPos(ig.GetFontSize() * 35.0)
        ig.TextUnformatted(desc)
        ig.PopTextWrapPos()
        ig.EndTooltip()
    end
    return nothing
end
