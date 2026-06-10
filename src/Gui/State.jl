function _project_for_preference(pref::AbstractString)
    pref == "auto" && return something(_default_project[])
    for project in KNOWN_PROJECTS
        project_name(project) == pref && return project
    end
    error("Unknown project preference '$pref'")
end

function _open_project_path!(ui_state, path::String; persist=true)
    norm_path = _normalize_project_path(path)
    ui_state[:project_preference] = _project_preference_for_path(ui_state, norm_path)
    proj = _project_for_preference(ui_state[:project_preference])
    cache_id = _cache_id_for_path!(ui_state, norm_path)
    _load_project_view!(ui_state, norm_path, proj)
    _load_tag_state_for_root!(ui_state, norm_path)
    _launch_project_reload_job!(ui_state, norm_path, proj, cache_id; persist)
    _launch_source_scan_job!(ui_state, norm_path, proj; persist=false, replace=true)
end

function _project_status_text(ui_state)
    if haskey(ui_state, :project)
        return "Active: $(project_name(ui_state[:project]))"
    end
    return "No project loaded"
end

function _new_scan_progress()
    return Dict{Symbol,Any}(
        :phase => :idle,
        :total_csv => 0,
        :processed_csv => 0,
        :loaded_measurements => 0,
        :skipped_csv => 0,
        :current_path => "",
    )
end

function _init_scan_state!(ui_state)
    ui_state[:source_scan_state] = :idle
    ui_state[:source_scan_progress] = _new_scan_progress()
    ui_state[:source_scan_error] = ""
    ui_state[:source_scan_seq] = 0
    ui_state[:active_source_scan_id] = 0
    ui_state[:source_scan_events] = nothing
    ui_state[:source_scan_cancel_token] = nothing
    ui_state[:source_scan_path] = ""
    ui_state[:source_scan_result] = nothing
end

function _init_cache_state!(ui_state)
    ui_state[:cache_state] = :idle
    ui_state[:cache_progress] = _new_scan_progress()
    ui_state[:cache_error] = ""
    ui_state[:cache_id] = ""
    ui_state[:cache_identity] = nothing
    ui_state[:cache_status] = nothing
    ui_state[:cache_source_checked] = false
    ui_state[:cache_errors] = Pair{String,String}[]
    ui_state[:cache_seq] = 0
    ui_state[:active_cache_id] = 0
    ui_state[:cache_events] = nothing
    ui_state[:cache_cancel_token] = nothing
    set_active_project_cache!(nothing)
end

function _init_tag_state!(ui_state)
    ui_state[:show_bad] = true
    ui_state[:tag_state] = nothing
    ui_state[:tag_state_error] = ""
    ui_state[:expanded_device_paths] = String[]
    ui_state[:selected_device_paths] = String[]
    ui_state[:selected_measurement_ids] = String[]
    ui_state[:selected_devices] = HierarchyNode[]
    ui_state[:selected_measurements] = MeasurementInfo[]
    ui_state[:selected_all_measurements] = MeasurementInfo[]
    ui_state[:selected_measurement_id_set] = Set{String}()
    ui_state[:selected_path] = String[]
    ui_state[:measurement_index] = Dict{String,MeasurementInfo}()
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
    path = String(filepath)
    measurements = [
        measurement
        for measurement in get(ui_state, :all_measurements, MeasurementInfo[])
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
    ui_state[:selected_device_paths] = device_paths
    ui_state[:selected_measurement_ids] = [measurement.unique_id for measurement in measurements]
    ui_state[:scroll_to_device_path] = first(device_paths)
    ui_state[:scroll_to_measurement_id] = first(measurements).unique_id
    _apply_visible_selection!(ui_state)
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
    tasks = get!(ui_state, :background_tasks) do
        Task[]
    end
    filter!(!istaskdone, tasks)
    push!(tasks, task)
    return task
end

function _request_cancel_token!(ui_state, key::Symbol)::Nothing
    token = get(ui_state, key, nothing)
    token === nothing && return nothing
    Base.Threads.atomic_xchg!(token, true)
    return nothing
end

function _request_background_jobs_cancel!(ui_state)::Nothing
    ui_state[:shutting_down] = true
    _request_cancel_token!(ui_state, :source_scan_cancel_token)
    _request_cancel_token!(ui_state, :cache_cancel_token)
    _source_scan_running(ui_state) && (ui_state[:source_scan_state] = :canceling)
    _cache_action_blocked(ui_state) && (ui_state[:cache_state] = :canceling)
    return nothing
end

function _wait_for_background_jobs!(ui_state)::Nothing
    tasks = get(ui_state, :background_tasks, Task[])
    tasks isa Vector{Task} || return nothing
    while true
        filter!(!istaskdone, tasks)
        isempty(tasks) && return nothing
        sleep(0.02)
    end
end

function _shutdown_background_jobs!(ui_state)::Nothing
    get(ui_state, :shutdown_complete, false) && return nothing
    _request_background_jobs_cancel!(ui_state)
    _wait_for_background_jobs!(ui_state)
    ui_state[:shutdown_complete] = true
    return nothing
end

function _clear_selection!(ui_state)
    ui_state[:selected_device_paths] = String[]
    ui_state[:selected_measurement_ids] = String[]
    ui_state[:selected_devices] = HierarchyNode[]
    ui_state[:selected_measurements] = MeasurementInfo[]
    ui_state[:selected_all_measurements] = MeasurementInfo[]
    ui_state[:selected_measurement_id_set] = Set{String}()
    ui_state[:selected_path] = String[]
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
    return Dict{Symbol,Any}(
        :operation => operation,
        :groups => copy(groups),
        :selected_measurements => copy(selected_measurements),
        :all_measurements => all_measurements,
        :group_name => String(group_name),
        :group_index => group_index,
        :scan_id => get(ui_state, :source_scan_seq, 0),
        :root_path => get(ui_state, :root_path, ""),
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
            if msg.scan_id == get(ui_state, :source_scan_seq, 0) &&
               msg.root_path == get(ui_state, :root_path, "")
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
    return get(ui_state, :all_measurements, MeasurementInfo[])
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
    all_measurements = _selected_measurements(ui_state)
    visible_measurements = if filter_meas === nothing
        all_measurements
    else
        proj = ui_state[:project]
        _visible_measurements(ui_state, proj, all_measurements, filter_meas)
    end
    selected_ids = Set(get(ui_state, :selected_measurement_ids, String[]))
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

function _scan_running(ui_state)
    return _source_scan_running(ui_state)
end

function _cache_action_blocked(ui_state)
    return get(ui_state, :cache_state, :idle) in (:loading, :writing, :canceling)
end

function _source_scan_running(ui_state)
    return get(ui_state, :source_scan_state, :idle) in (:counting, :discovering, :scanning, :analyzing, :canceling)
end

function _request_cache_cancel!(ui_state)
    token = get(ui_state, :cache_cancel_token, nothing)
    token === nothing && return
    Base.Threads.atomic_xchg!(token, true)
    ui_state[:cache_state] = :canceling
end

function _request_source_scan_cancel!(ui_state)
    token = get(ui_state, :source_scan_cancel_token, nothing)
    token === nothing && return
    Base.Threads.atomic_xchg!(token, true)
    ui_state[:source_scan_state] = :canceling
end

function _begin_scan!(
    ui_state,
    path::String,
    proj::AbstractProject,
    has_device_metadata::Bool=_has_device_metadata(path),
)
    current_project = get(ui_state, :project, nothing)
    same_project =
        get(ui_state, :root_path, "") == path &&
        current_project isa AbstractProject &&
        project_name(current_project) == project_name(proj)

    _clear_plot_views!(ui_state)
    if !same_project
        _clear_selection!(ui_state)
        delete!(ui_state, :main_plot_kind)
        delete!(ui_state, :main_plot_measurement_kind)
        ui_state[:main_plot_measurement_ids] = String[]
        ui_state[:main_plot_measurements] = MeasurementInfo[]
        ui_state[:open_plot_windows] = Vector{Dict{Symbol,Any}}()
    end
    if get(ui_state, :figure_script_root_path, "") != path
        _reset_figure_script_state!(ui_state, path)
    else
        ui_state[:figure_script_root_path] = path
    end

    hierarchy = MeasurementHierarchy(path, has_device_metadata, proj, 0)
    ui_state[:scan_hierarchy] = hierarchy
    ui_state[:hierarchy_root] = hierarchy.root
    ui_state[:all_measurements] = hierarchy.all_measurements
    ui_state[:root_path] = path
    ui_state[:has_device_metadata] = has_device_metadata
    ui_state[:project] = proj
    ui_state[:skipped_count] = 0
    ui_state[:device_metadata_keys] = Symbol[]
    ui_state[:measurement_index] = Dict{String,MeasurementInfo}()
    if !same_project
        view = get(ui_state, :project_view_loaded, PersistedProjectView(project=project_name(proj)))
        _apply_project_view!(ui_state, view)
    else
        _apply_visible_selection!(ui_state)
        _refresh_plot_measurement_refs!(ui_state)
    end
    _invalidate_figure_script_scan_cache!(ui_state)
end

function _build_measurement_index(measurements::Vector{MeasurementInfo})
    measurement_index = Dict{String,MeasurementInfo}()
    sizehint!(measurement_index, length(measurements))
    for measurement in measurements
        haskey(measurement_index, measurement.unique_id) && error(
            "Duplicate measurement id generated during scan: $(measurement.unique_id)",
        )
        measurement_index[measurement.unique_id] = measurement
    end
    return measurement_index
end

function _device_metadata_keys(measurements::Vector{MeasurementInfo})
    metadata_keys = Set{Symbol}()
    for measurement in measurements
        for key in keys(measurement.device_info.parameters)
            push!(metadata_keys, key)
        end
    end
    return sort!(collect(metadata_keys); by=String)
end

function _apply_scan_snapshot!(ui_state, snapshot)
    hierarchy = snapshot.hierarchy
    ui_state[:scan_hierarchy] = hierarchy
    ui_state[:hierarchy_root] = hierarchy.root
    ui_state[:all_measurements] = hierarchy.all_measurements
    ui_state[:measurement_index] = snapshot.measurement_index
    ui_state[:device_metadata_keys] = snapshot.device_metadata_keys
    ui_state[:skipped_count] = snapshot.skipped_count
    ui_state[:has_device_metadata] = hierarchy.has_device_metadata || _has_device_metadata(hierarchy.root_path)
    _apply_visible_selection!(ui_state)
    _refresh_plot_measurement_refs!(ui_state)
    _invalidate_figure_script_scan_cache!(ui_state)
end

function _apply_cache_snapshot!(ui_state, index::ProjectCacheIndex)
    identity = index.identity
    hierarchy = index.source.hierarchy
    _apply_scan_snapshot!(ui_state, (
        hierarchy,
        measurement_index=_build_measurement_index(hierarchy.all_measurements),
        device_metadata_keys=_device_metadata_keys(hierarchy.all_measurements),
        skipped_count=hierarchy.skipped_count,
    ))
    ui_state[:root_path] = identity.root_path
    ui_state[:project] = index.source.project
    ui_state[:cache_id] = identity.cache_id
    ui_state[:cache_identity] = identity
    set_active_project_cache!(index)
    source = get(ui_state, :source_scan_result, nothing)
    if source isa SourceScan && source.root_path == identity.root_path
        ui_state[:cache_status] = cache_status(index, source)
        ui_state[:cache_source_checked] = true
    else
        cached_files = length(index.files)
        ui_state[:cache_status] = ProjectCacheStatus(
            cached_files,
            cached_files,
            0,
            0,
            0,
            0,
            length(index.analysis_errors),
        )
        ui_state[:cache_source_checked] = false
    end
    ui_state[:cache_errors] = sort!(collect(index.analysis_errors); by=first)
    ui_state[:cache_state] = :ready
    ui_state[:cache_error] = ""
    return nothing
end

function _apply_source_scan!(ui_state, source::SourceScan)
    ui_state[:source_scan_result] = source
    _apply_scan_snapshot!(ui_state, (
        hierarchy=source.hierarchy,
        measurement_index=_build_measurement_index(source.hierarchy.all_measurements),
        device_metadata_keys=_device_metadata_keys(source.hierarchy.all_measurements),
        skipped_count=source.hierarchy.skipped_count,
    ))
    ui_state[:root_path] = source.root_path
    ui_state[:project] = source.project
    ui_state[:has_device_metadata] = source.hierarchy.has_device_metadata
    identity = get(ui_state, :cache_identity, nothing)
    if identity isa ProjectCacheIdentity && identity.root_path == source.root_path
        ui_state[:cache_errors] =
            sort!(collect(ProjectCacheIndex(identity, source).analysis_errors); by=first)
        cache_state = get(ui_state, :cache_state, :idle)
        if cache_state == :ready
            ui_state[:cache_status] = cache_status(project_cache_index(identity), source)
            ui_state[:cache_source_checked] = true
        elseif !isfile(identity.cache_path) || cache_state == :missing
            ui_state[:cache_status] = ProjectCacheStatus(
                length(source.files),
                0,
                0,
                0,
                length(source.files),
                0,
                0,
            )
            ui_state[:cache_source_checked] = true
        else
            ui_state[:cache_source_checked] = false
        end
    else
        set_active_project_cache!(nothing)
    end
    return nothing
end

function _cache_status_needs_update(status::ProjectCacheStatus)::Bool
    return status.stale_files > 0 || status.new_files > 0 || status.deleted_files > 0
end

function _maybe_update_cache_from_source!(ui_state)::Nothing
    _cache_action_blocked(ui_state) && return nothing
    source = get(ui_state, :source_scan_result, nothing)
    source isa SourceScan || return nothing
    identity = get(ui_state, :cache_identity, nothing)
    identity isa ProjectCacheIdentity || return nothing
    identity.root_path == source.root_path || return nothing
    project_name(source.project) == identity.project_name || return nothing

    cache_state = get(ui_state, :cache_state, :idle)
    if cache_state == :missing
        _launch_cache_update_job!(ui_state; full_rebuild=true)
        return nothing
    end

    status = get(ui_state, :cache_status, nothing)
    status isa ProjectCacheStatus || return nothing
    get(ui_state, :cache_source_checked, false) || return nothing
    _cache_status_needs_update(status) || return nothing
    _launch_cache_update_job!(ui_state)
    return nothing
end

function _append_measurements!(ui_state, measurements::Vector{MeasurementInfo})
    isempty(measurements) && return nothing
    hierarchy = get(ui_state, :scan_hierarchy, nothing)
    hierarchy isa MeasurementHierarchy || return nothing
    measurement_index = get!(ui_state, :measurement_index) do
        Dict{String,MeasurementInfo}()
    end
    metadata_keys = Set(get(ui_state, :device_metadata_keys, Symbol[]))
    appended = false
    for measurement in measurements
        haskey(measurement_index, measurement.unique_id) && continue
        insert_measurement!(hierarchy, measurement)
        measurement_index[measurement.unique_id] = measurement
        foreach(key -> push!(metadata_keys, key), keys(measurement.device_info.parameters))
        appended = true
    end
    appended || return nothing
    sort!(hierarchy)
    ui_state[:scan_hierarchy] = hierarchy
    ui_state[:hierarchy_root] = hierarchy.root
    ui_state[:all_measurements] = hierarchy.all_measurements
    ui_state[:measurement_index] = measurement_index
    ui_state[:device_metadata_keys] = sort!(collect(metadata_keys); by=String)
    _apply_visible_selection!(ui_state)
    _refresh_plot_measurement_refs!(ui_state)
    _invalidate_figure_script_scan_cache!(ui_state)
    return nothing
end

function _launch_source_scan_job!(
    ui_state,
    path::String,
    proj::AbstractProject;
    persist::Bool=true,
    replace::Bool=false,
)
    if _source_scan_running(ui_state)
        _request_source_scan_cancel!(ui_state)
        replace || return
    end

    norm_path = _normalize_project_path(path)
    cache_id = get(ui_state, :cache_id, "")
    if isempty(cache_id)
        cache_id = _cache_id_for_path!(ui_state, norm_path)
    end
    identity = project_cache_identity(cache_id, proj, norm_path)
    ui_state[:cache_id] = cache_id
    ui_state[:cache_identity] = identity
    _load_tag_state_for_root!(ui_state, norm_path)
    ui_state[:source_scan_seq] = get(ui_state, :source_scan_seq, 0) + 1
    scan_id = ui_state[:source_scan_seq]
    ui_state[:active_source_scan_id] = scan_id
    ui_state[:source_scan_path] = norm_path
    ui_state[:source_scan_state] = :discovering
    ui_state[:source_scan_progress] = _new_scan_progress()
    ui_state[:source_scan_error] = ""
    _begin_scan!(ui_state, norm_path, proj, _has_device_metadata(norm_path))

    events = Channel{NamedTuple}(Inf)
    cancel_token = Base.Threads.Atomic{Bool}(false)
    ui_state[:source_scan_events] = events
    ui_state[:source_scan_cancel_token] = cancel_token

    task = Base.Threads.@spawn begin
        try
            source = _with_cancel(() -> cancel_token[]) do
                scan_source(
                    norm_path;
                    project=proj,
                    on_progress=(progress) -> put!(events, (
                        kind=:progress,
                        scan_id=scan_id,
                        progress=progress,
                    )),
                    on_measurements=(measurements) -> put!(events, (
                        kind=:measurements,
                        scan_id=scan_id,
                        measurements=measurements,
                    )),
                )
            end
            put!(events, (
                kind=:source_scan_result,
                scan_id=scan_id,
                identity=identity,
                source=source,
                persist=persist,
            ))
        catch err
            if _is_job_cancelled(err)
                put!(events, (kind=:canceled, scan_id=scan_id))
            else
                put!(events, (kind=:error, scan_id=scan_id, error=err, bt=catch_backtrace()))
            end
        finally
            close(events)
        end
    end
    _remember_background_task!(ui_state, task)
    return nothing
end

function _launch_project_reload_job!(
    ui_state,
    path::String,
    proj::AbstractProject,
    cache_id::String;
    persist::Bool=true,
)
    if _source_scan_running(ui_state)
        _request_source_scan_cancel!(ui_state)
    end
    if _cache_action_blocked(ui_state)
        _request_cache_cancel!(ui_state)
    end

    identity = project_cache_identity(cache_id, proj, path)
    current_root = get(ui_state, :root_path, "")
    if current_root != identity.root_path || !haskey(ui_state, :scan_hierarchy)
        _begin_scan!(ui_state, identity.root_path, proj, _has_device_metadata(identity.root_path))
    else
        ui_state[:root_path] = identity.root_path
        ui_state[:project] = proj
    end

    cache_id_num = get(ui_state, :cache_seq, 0) + 1
    ui_state[:cache_seq] = cache_id_num
    ui_state[:active_cache_id] = cache_id_num
    ui_state[:cache_progress] = _new_scan_progress()
    ui_state[:cache_id] = cache_id
    ui_state[:cache_identity] = identity
    set_active_project_cache!(nothing)
    ui_state[:cache_state] = :loading
    ui_state[:cache_error] = ""
    ui_state[:cache_status] = nothing
    ui_state[:cache_errors] = Pair{String,String}[]

    events = Channel{NamedTuple}(Inf)
    cancel_token = Base.Threads.Atomic{Bool}(false)
    ui_state[:cache_events] = events
    ui_state[:cache_cancel_token] = cancel_token

    task = Base.Threads.@spawn begin
        try
            snapshot = _with_cancel(() -> cancel_token[]) do
                load_project_cache(
                    identity;
                    on_progress=(progress) -> put!(events, (
                        kind=:progress,
                        cache_id=cache_id_num,
                        progress=progress,
                    )),
                )
            end
            put!(events, (kind=:cache_result, cache_id=cache_id_num, snapshot=snapshot, persist=persist))
        catch err
            if _is_job_cancelled(err)
                put!(events, (kind=:canceled, cache_id=cache_id_num))
            elseif err isa ProjectCacheError
                put!(events, (
                    kind=:cache_missing,
                    cache_id=cache_id_num,
                    identity=identity,
                    error=sprint(showerror, err),
                    persist=persist,
                ))
            else
                put!(events, (kind=:error, cache_id=cache_id_num, error=err, bt=catch_backtrace()))
            end
        finally
            close(events)
        end
    end
    _remember_background_task!(ui_state, task)
    return nothing
end

function _launch_cache_update_job!(ui_state; full_rebuild::Bool=false)
    identity = get(ui_state, :cache_identity, nothing)
    identity isa ProjectCacheIdentity || error("No project cache is bound to the current project")
    project = get(ui_state, :project, nothing)
    project isa AbstractProject || error("No project is open")
    previous_cache_state = get(ui_state, :cache_state, :idle)
    _cancel_figure_script_job!(ui_state)
    _clear_plot_views!(ui_state)

    cache_id_num = get(ui_state, :cache_seq, 0) + 1
    ui_state[:cache_seq] = cache_id_num
    ui_state[:active_cache_id] = cache_id_num
    ui_state[:cache_progress] = _new_scan_progress()
    ui_state[:cache_state] = :writing
    operation = full_rebuild ? (previous_cache_state == :missing ? :build : :rebuild) : :update
    ui_state[:cache_operation] = operation
    ui_state[:cache_error] = ""

    events = Channel{NamedTuple}(Inf)
    cancel_token = Base.Threads.Atomic{Bool}(false)
    ui_state[:cache_events] = events
    ui_state[:cache_cancel_token] = cancel_token
    existing_source = get(ui_state, :source_scan_result, nothing)
    if !(existing_source isa SourceScan) ||
       existing_source.root_path != identity.root_path ||
       project_name(existing_source.project) != identity.project_name
        existing_source = nothing
    end

    task = Base.Threads.@spawn begin
        try
            source = existing_source
            if source === nothing
                source = _with_cancel(() -> cancel_token[]) do
                    scan_source(
                        identity.root_path;
                        project,
                        on_progress=(progress) -> put!(events, (
                            kind=:source_progress,
                            cache_id=cache_id_num,
                            progress=progress,
                        )),
                        on_measurements=(measurements) -> put!(events, (
                            kind=:source_measurements,
                            cache_id=cache_id_num,
                            measurements=measurements,
                        )),
                    )
                end
                put!(events, (kind=:source_scan_result, cache_id=cache_id_num, source=source))
            end
            snapshot = _with_cancel(() -> cancel_token[]) do
                write_project_cache!(
                    identity,
                    source;
                    replace=full_rebuild,
                    on_progress=(progress) -> put!(events, (
                        kind=:progress,
                        cache_id=cache_id_num,
                        progress=progress,
                    )),
                )
            end
            put!(events, (kind=:cache_result, cache_id=cache_id_num, snapshot=snapshot, persist=true))
        catch err
            if _is_job_cancelled(err)
                put!(events, (kind=:canceled, cache_id=cache_id_num))
            else
                put!(events, (kind=:error, cache_id=cache_id_num, error=err, bt=catch_backtrace()))
            end
        finally
            close(events)
        end
    end
    _remember_background_task!(ui_state, task)
end

function _queue_cache_update!(ui_state; full_rebuild::Bool=false)
    if _cache_action_blocked(ui_state)
        _request_cache_cancel!(ui_state)
        return
    end
    _launch_cache_update_job!(ui_state; full_rebuild)
end

function _finalize_source_scan!(ui_state)
    ui_state[:source_scan_events] = nothing
    ui_state[:source_scan_cancel_token] = nothing
end

function _finalize_cache!(ui_state)
    ui_state[:cache_events] = nothing
    ui_state[:cache_cancel_token] = nothing
end

function _poll_cache_events!(ui_state)
    events = get(ui_state, :cache_events, nothing)
    events === nothing && return

    pending_measurements = MeasurementInfo[]
    while isready(events)
        msg = take!(events)
        msg.cache_id == get(ui_state, :active_cache_id, 0) || continue

        if msg.kind == :progress
            p = msg.progress
            progress_dict = Dict{Symbol,Any}(
                :phase => p.phase,
                :total_csv => p.total_csv,
                :processed_csv => p.processed_csv,
                :loaded_measurements => p.loaded_measurements,
                :skipped_csv => p.skipped_csv,
                :current_path => p.current_path,
            )
            ui_state[:cache_progress] = progress_dict
        elseif msg.kind == :source_progress
            p = msg.progress
            ui_state[:source_scan_progress] = Dict{Symbol,Any}(
                :phase => p.phase,
                :total_csv => p.total_csv,
                :processed_csv => p.processed_csv,
                :loaded_measurements => p.loaded_measurements,
                :skipped_csv => p.skipped_csv,
                :current_path => p.current_path,
            )
            ui_state[:source_scan_state] = p.phase
        elseif msg.kind == :source_measurements
            append!(pending_measurements, msg.measurements)
        elseif msg.kind == :source_scan_result
            _apply_source_scan!(ui_state, msg.source)
            ui_state[:source_scan_state] = :done
            _maybe_update_cache_from_source!(ui_state)
        elseif msg.kind == :cache_result
            _apply_cache_snapshot!(ui_state, msg.snapshot)
            get(msg, :persist, false) &&
                _persist_preferences!(ui_state; path=msg.snapshot.identity.root_path)
            _finalize_cache!(ui_state)
            _maybe_update_cache_from_source!(ui_state)
        elseif msg.kind == :cache_missing
            identity = msg.identity
            proj = get(ui_state, :project, nothing)
            proj isa AbstractProject || error("No project is open")
            current_root = get(ui_state, :root_path, "")
            if current_root != identity.root_path || !haskey(ui_state, :scan_hierarchy)
                _begin_scan!(ui_state, identity.root_path, proj, _has_device_metadata(identity.root_path))
            end
            ui_state[:cache_id] = identity.cache_id
            ui_state[:cache_identity] = identity
            ui_state[:cache_state] = :missing
            ui_state[:cache_error] = msg.error
            msg.persist && _persist_preferences!(ui_state; path=identity.root_path)
            _finalize_cache!(ui_state)
            _maybe_update_cache_from_source!(ui_state)
        elseif msg.kind == :canceled
            ui_state[:cache_state] = :canceled
            _finalize_cache!(ui_state)
        elseif msg.kind == :error
            ui_state[:cache_state] = :error
            ui_state[:cache_error] = "Cache job failed. See the console for full details."
            identity = get(ui_state, :cache_identity, nothing)
            @error(
                "Cache job failed",
                cache_path=identity isa ProjectCacheIdentity ? identity.cache_path : "",
                operation=get(ui_state, :cache_operation, :load),
                exception=(msg.error, msg.bt),
            )
            _finalize_cache!(ui_state)
        end
    end
    _append_measurements!(ui_state, pending_measurements)
end

function _poll_source_scan_events!(ui_state)
    events = get(ui_state, :source_scan_events, nothing)
    events === nothing && return

    pending_measurements = MeasurementInfo[]
    while isready(events)
        msg = take!(events)
        msg.scan_id == get(ui_state, :active_source_scan_id, 0) || continue

        if msg.kind == :progress
            p = msg.progress
            progress_dict = Dict{Symbol,Any}(
                :phase => p.phase,
                :total_csv => p.total_csv,
                :processed_csv => p.processed_csv,
                :loaded_measurements => p.loaded_measurements,
                :skipped_csv => p.skipped_csv,
                :current_path => p.current_path,
            )
            ui_state[:source_scan_progress] = progress_dict
            ui_state[:source_scan_state] = p.phase
        elseif msg.kind == :measurements
            append!(pending_measurements, msg.measurements)
        elseif msg.kind == :source_scan_result
            ui_state[:cache_identity] = msg.identity
            ui_state[:cache_id] = msg.identity.cache_id
            _apply_source_scan!(ui_state, msg.source)
            ui_state[:source_scan_state] = :done
            msg.persist && _persist_preferences!(ui_state; path=msg.identity.root_path)
            _finalize_source_scan!(ui_state)
            _maybe_update_cache_from_source!(ui_state)
        elseif msg.kind == :canceled
            ui_state[:source_scan_state] = :canceled
            _finalize_source_scan!(ui_state)
        elseif msg.kind == :error
            ui_state[:source_scan_state] = :error
            ui_state[:source_scan_error] =
                "Source scan failed. See the console for full details."
            progress = get(ui_state, :source_scan_progress, _new_scan_progress())
            project = get(ui_state, :project, nothing)
            @error(
                "Source scan job failed",
                project=project isa AbstractProject ? project_name(project) : "",
                root=get(ui_state, :root_path, ""),
                current_file=get(progress, :current_path, ""),
                exception=(msg.error, msg.bt),
            )
            _finalize_source_scan!(ui_state)
        end
    end
    _append_measurements!(ui_state, pending_measurements)
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
