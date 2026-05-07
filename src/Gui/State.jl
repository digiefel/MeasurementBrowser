# ---------------------------------------------------------------------------
# Preferences (persistent project selection + recent folders)
# ---------------------------------------------------------------------------

function _prefs_path()
    return joinpath(homedir(), ".config", "MeasurementBrowser", "prefs.toml")
end

function _load_prefs()
    path = _prefs_path()
    isfile(path) || return Dict{String,Any}()
    return TOML.parsefile(path)
end

function _save_prefs(data::Dict)
    path = _prefs_path()
    mkpath(dirname(path))
    open(path, "w") do io
        TOML.print(io, data)
    end
end

const _MAX_RECENT_PROJECTS = 12

function _normalize_project_path(path::AbstractString)
    norm_path = abspath(expanduser(String(path)))
    return ispath(norm_path) ? realpath(norm_path) : norm_path
end

function _sanitize_project_preference(pref::String)
    pref == "auto" && return "auto"
    for p in KNOWN_PROJECTS
        project_name(p) == pref && return pref
    end
    return "auto"
end

function _sanitize_figure_script_output_dir(value)
    value isa AbstractString || return ""
    return String(strip(String(value)))
end

function _sanitize_cache_id(value)
    value isa AbstractString || return ""
    stripped = strip(String(value))
    isempty(stripped) && return ""
    occursin(r"^[A-Za-z0-9_.-]+$", stripped) ||
        error("Invalid cache id '$stripped'; expected letters, numbers, '.', '_' or '-'")
    return stripped
end

function _parse_recent_projects(prefs::Dict{String,Any})
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
        cache_id = _sanitize_cache_id(get(entry, "cache_id", ""))
        push!(recents, Dict{String,String}(
            "path" => _normalize_project_path(path),
            "project_preference" => _sanitize_project_preference(pref),
            "figure_script_output_dir" => figure_script_output_dir,
            "cache_id" => cache_id,
        ))
    end

    return recents
end

function _update_recent_projects(
    recents::Vector{Dict{String,String}},
    path::AbstractString,
    pref::AbstractString,
    figure_script_output_dir::AbstractString,
    cache_id::AbstractString,
)
    norm_path = _normalize_project_path(path)
    filter!(entry -> get(entry, "path", "") != norm_path, recents)
    pushfirst!(recents, Dict{String,String}(
        "path" => norm_path,
        "project_preference" => String(pref),
        "figure_script_output_dir" => _sanitize_figure_script_output_dir(figure_script_output_dir),
        "cache_id" => _sanitize_cache_id(cache_id),
    ))
    length(recents) > _MAX_RECENT_PROJECTS && resize!(recents, _MAX_RECENT_PROJECTS)
    return recents
end

function _current_figure_script_output_dir(ui_state)
    haskey(ui_state, :figure_script_output_dir_buffer) || return ""
    return _sanitize_figure_script_output_dir(_buffer_string(ui_state[:figure_script_output_dir_buffer]))
end

function _persist_preferences!(ui_state; path::Union{Nothing,String}=nothing)
    prefs = _load_prefs()
    pref = _sanitize_project_preference(string(get(ui_state, :project_preference, "auto")))
    ui_state[:project_preference] = pref
    prefs["project"] = pref

    recents = _parse_recent_projects(prefs)
    if path !== nothing && !isempty(path)
        cache_id = string(get(ui_state, :cache_id, ""))
        _update_recent_projects(recents, path, pref, _current_figure_script_output_dir(ui_state), cache_id)
        prefs["recent_projects"] = recents
    end

    _save_prefs(prefs)
    ui_state[:recent_projects] = recents
end

function _recent_project_entry_for_path(ui_state, path::String)
    recents = get(ui_state, :recent_projects, Dict{String,String}[])
    for entry in recents
        get(entry, "path", "") == path && return entry
    end
    return nothing
end

function _project_preference_for_path(ui_state, path::String)
    entry = _recent_project_entry_for_path(ui_state, path)
    if entry !== nothing
        pref = get(entry, "project_preference", "auto")
        return _sanitize_project_preference(pref)
    end
    pref = string(get(ui_state, :project_preference, "auto"))
    return _sanitize_project_preference(pref)
end

function _figure_script_output_dir_for_path(ui_state, path::String)
    entry = _recent_project_entry_for_path(ui_state, path)
    entry === nothing && return ""
    return _sanitize_figure_script_output_dir(get(entry, "figure_script_output_dir", ""))
end

function _cache_id_for_path!(ui_state, path::String)
    cache_id = project_cache_id(path)
    pref = _project_preference_for_path(ui_state, path)
    recents = get!(ui_state, :recent_projects) do
        Dict{String,String}[]
    end
    _update_recent_projects(
        recents,
        path,
        pref,
        _figure_script_output_dir_for_path(ui_state, path),
        cache_id,
    )
    prefs = _load_prefs()
    prefs["recent_projects"] = recents
    prefs["project"] = pref
    _save_prefs(prefs)
    return cache_id
end

function _persist_current_project_preferences!(ui_state)
    current_root = get(ui_state, :root_path, "")
    isempty(current_root) && return
    _persist_preferences!(ui_state; path=current_root)
end

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
    _load_bad_registry_for_root!(ui_state, norm_path)
    _launch_project_reload_job!(ui_state, norm_path, proj, cache_id; persist)
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
    ui_state[:scan_state] = :idle
    ui_state[:scan_progress] = _new_scan_progress()
    ui_state[:scan_error] = ""
    ui_state[:scan_seq] = 0
    ui_state[:active_scan_id] = 0
    ui_state[:scan_events] = nothing
    ui_state[:scan_cancel_token] = nothing
    ui_state[:scan_path] = ""
    ui_state[:source_scan_state] = :idle
    ui_state[:source_scan_progress] = _new_scan_progress()
    ui_state[:source_scan_error] = ""
    ui_state[:source_scan_seq] = 0
    ui_state[:active_source_scan_id] = 0
    ui_state[:source_scan_events] = nothing
    ui_state[:source_scan_cancel_token] = nothing
    ui_state[:source_scan_path] = ""
end

function _init_cache_state!(ui_state)
    ui_state[:cache_state] = :idle
    ui_state[:cache_error] = ""
    ui_state[:cache_id] = ""
    ui_state[:cache_identity] = nothing
    ui_state[:cache_status] = nothing
    ui_state[:cache_semantic_fields] = Dict{Symbol,Vector{Symbol}}()
    ui_state[:cache_errors] = ProjectCacheFileError[]
end

function _init_bad_state!(ui_state)
    ui_state[:show_bad] = true
    ui_state[:bad_registry] = nothing
    ui_state[:bad_registry_error] = ""
    ui_state[:selected_device_paths] = String[]
    ui_state[:selected_measurement_ids] = String[]
    ui_state[:selected_devices] = HierarchyNode[]
    ui_state[:selected_measurements] = MeasurementInfo[]
    ui_state[:selected_all_measurements] = MeasurementInfo[]
    ui_state[:selected_measurement_id_set] = Set{String}()
    ui_state[:selected_path] = String[]
    ui_state[:measurement_index] = Dict{String,MeasurementInfo}()
    ui_state[:computed_stats_cache] = Dict{Tuple{String,Int},Dict{Symbol,Any}}()
end

function _init_plot_state!(ui_state)
    ui_state[:plot_state] = :idle
    ui_state[:plot_error] = ""
    ui_state[:plot_seq] = 0
    ui_state[:active_plot_id] = 0
    ui_state[:plot_events] = nothing
    ui_state[:plot_cancel_token] = nothing
    ui_state[:active_plot_job] = nothing
    ui_state[:pending_plot_job] = nothing
    ui_state[:plot_runtime_warmed] = false
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
        :scan_id => get(ui_state, :active_scan_id, 0),
        :root_path => get(ui_state, :root_path, ""),
        :status_message => String(status_message),
        :completion_message => String(completion_message),
    )
end

function _resolve_figure_script_job(request::Dict{Symbol,Any})
    inferred_group, profile = _infer_measurement_group_profiled(
        request[:group_name],
        request[:selected_measurements],
        request[:all_measurements],
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

    Base.Threads.@spawn begin
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
end

function _poll_figure_script_job_events!(ui_state)
    events = get(ui_state, :figure_script_job_events, nothing)
    events === nothing && return

    while isready(events)
        msg = take!(events)
        msg.job_id == get(ui_state, :active_figure_script_job_id, 0) || continue

        if msg.kind == :done
            _record_figure_script_job_profile!(ui_state, msg.operation, msg.profile)
            if msg.scan_id == get(ui_state, :active_scan_id, 0) &&
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
    return [measurement for measurement in visible_measurements if measurement.id in selected_ids]
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
    merged_ids = Set(measurement.id for measurement in _group_measurements_in_current_scan(ui_state, group))
    foreach(measurement -> push!(merged_ids, measurement.id), selected_measurements)
    merged_measurements = [measurement for measurement in all_measurements if measurement.id in merged_ids]
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

    remaining_ids = Set(measurement.id for measurement in _group_measurements_in_current_scan(ui_state, group))
    foreach(measurement -> delete!(remaining_ids, measurement.id), selected_measurements)
    isempty(remaining_ids) && throw(FigureScriptValidationError("Measurement groups cannot be empty"))
    all_measurements = _current_scan_measurements(ui_state)
    remaining_measurements = [measurement for measurement in all_measurements if measurement.id in remaining_ids]
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
    return get(ui_state, :scan_state, :idle) in (
        :cache_reload,
        :cache_load,
        :cache_check,
        :cache_discovery,
        :cache_update,
        :canceling,
    )
end

function _cache_action_blocked(ui_state)
    return get(ui_state, :scan_state, :idle) in (
        :cache_reload,
        :cache_load,
        :cache_discovery,
        :cache_update,
        :canceling,
    )
end

function _source_scan_running(ui_state)
    return get(ui_state, :source_scan_state, :idle) in (:counting, :cache_check, :canceling)
end

function _request_scan_cancel!(ui_state)
    token = get(ui_state, :scan_cancel_token, nothing)
    token === nothing && return
    Base.Threads.atomic_xchg!(token, true)
    ui_state[:scan_state] = :canceling
end

function _request_source_scan_cancel!(ui_state)
    token = get(ui_state, :source_scan_cancel_token, nothing)
    token === nothing && return
    Base.Threads.atomic_xchg!(token, true)
    ui_state[:source_scan_state] = :canceling
end

function _begin_scan!(ui_state, path::String, proj::AbstractProject, has_device_metadata::Bool)
    _clear_plot_jobs!(ui_state)
    _clear_selection!(ui_state)
    delete!(ui_state, :plot_figure)
    delete!(ui_state, :_last_plot_key)
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
    _invalidate_figure_script_scan_cache!(ui_state)
end

function _build_measurement_index(measurements::Vector{MeasurementInfo})
    measurement_index = Dict{String,MeasurementInfo}()
    sizehint!(measurement_index, length(measurements))
    for measurement in measurements
        haskey(measurement_index, measurement.id) && error(
            "Duplicate measurement id generated during scan: $(measurement.id)",
        )
        measurement_index[measurement.id] = measurement
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
    _apply_visible_selection!(ui_state)
    _invalidate_figure_script_scan_cache!(ui_state)
end

function _apply_cache_snapshot!(ui_state, snapshot::ProjectCacheSnapshot)
    _apply_scan_snapshot!(ui_state, (
        hierarchy=snapshot.hierarchy,
        measurement_index=_build_measurement_index(snapshot.hierarchy.all_measurements),
        device_metadata_keys=_device_metadata_keys(snapshot.hierarchy.all_measurements),
        skipped_count=snapshot.hierarchy.skipped_count,
    ))
    ui_state[:root_path] = snapshot.identity.root_path
    ui_state[:project] = _project_by_name(snapshot.identity.project_name)
    ui_state[:cache_id] = snapshot.identity.cache_id
    ui_state[:cache_identity] = snapshot.identity
    ui_state[:cache_status] = snapshot.status
    ui_state[:cache_semantic_fields] = snapshot.semantic_fields
    ui_state[:cache_errors] = snapshot.errors
    ui_state[:cache_state] = :ready
    ui_state[:cache_error] = ""
    ui_state[:scan_state] = :done
    return nothing
end

function _apply_cache_metadata!(ui_state, metadata)
    identity = metadata.identity
    ui_state[:root_path] = identity.root_path
    ui_state[:project] = _project_by_name(identity.project_name)
    ui_state[:cache_id] = identity.cache_id
    ui_state[:cache_identity] = identity
    ui_state[:cache_status] = metadata.status
    ui_state[:cache_semantic_fields] = metadata.semantic_fields
    ui_state[:cache_errors] = metadata.errors
    ui_state[:skipped_count] = metadata.skipped_count
    ui_state[:has_device_metadata] = metadata.has_device_metadata
    ui_state[:cache_error] = ""
    return nothing
end

function _append_cache_measurements!(ui_state, measurements::Vector{MeasurementInfo})
    isempty(measurements) && return nothing
    hierarchy = get(ui_state, :scan_hierarchy, nothing)
    hierarchy isa MeasurementHierarchy || return nothing
    measurement_index = get!(ui_state, :measurement_index) do
        Dict{String,MeasurementInfo}()
    end
    metadata_keys = Set(get(ui_state, :device_metadata_keys, Symbol[]))
    appended = false
    for measurement in measurements
        haskey(measurement_index, measurement.id) && continue
        insert_measurement!(hierarchy, measurement)
        measurement_index[measurement.id] = measurement
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
    _invalidate_figure_script_scan_cache!(ui_state)
    return nothing
end

function _launch_source_scan_job!(
    ui_state,
    path::String,
    proj::AbstractProject;
    persist::Bool=true,
)
    if _source_scan_running(ui_state)
        _request_source_scan_cancel!(ui_state)
        return
    end

    norm_path = _normalize_project_path(path)
    cache_id = get(ui_state, :cache_id, "")
    if isempty(cache_id)
        cache_id = _cache_id_for_path!(ui_state, norm_path)
    end
    identity = project_cache_identity(cache_id, proj, norm_path)
    ui_state[:cache_id] = cache_id
    ui_state[:cache_identity] = identity
    _load_bad_registry_for_root!(ui_state, norm_path)
    ui_state[:source_scan_seq] = get(ui_state, :source_scan_seq, 0) + 1
    scan_id = ui_state[:source_scan_seq]
    ui_state[:active_source_scan_id] = scan_id
    ui_state[:source_scan_path] = norm_path
    ui_state[:source_scan_state] = :cache_check
    ui_state[:source_scan_progress] = _new_scan_progress()
    ui_state[:source_scan_error] = ""

    events = Channel{NamedTuple}(Inf)
    cancel_token = Base.Threads.Atomic{Bool}(false)
    ui_state[:source_scan_events] = events
    ui_state[:source_scan_cancel_token] = cancel_token

    Base.Threads.@spawn begin
        try
            cached_fingerprints, cached_statuses = try
                _cached_cache_index(identity)
            catch err
                err isa ProjectCacheMissingError || rethrow()
                (Dict{String,FileFingerprint}(), Dict{String,String}())
            end
            raw = collect_csv_fingerprints(
                norm_path;
                should_cancel=() -> cancel_token[],
                on_progress=(progress) -> put!(events, (
                    kind=:progress,
                    scan_id=scan_id,
                    progress=progress,
                )),
                phase=:cache_check,
            )
            status = _cache_status_from_fingerprints(raw, cached_fingerprints, cached_statuses)
            put!(events, (
                kind=:source_scan_result,
                scan_id=scan_id,
                identity=identity,
                status=status,
                persist=persist,
            ))
        catch err
            if _is_scan_cancel_error(err) || err isa PlotCancelled
                put!(events, (kind=:canceled, scan_id=scan_id))
            else
                put!(events, (kind=:error, scan_id=scan_id, error=err, bt=catch_backtrace()))
            end
        finally
            close(events)
        end
    end
    return nothing
end

function _launch_project_reload_job!(
    ui_state,
    path::String,
    proj::AbstractProject,
    cache_id::String;
    persist::Bool=true,
)
    if _scan_running(ui_state)
        _request_scan_cancel!(ui_state)
    end
    if _source_scan_running(ui_state)
        _request_source_scan_cancel!(ui_state)
    end

    identity = project_cache_identity(cache_id, proj, path)
    current_root = get(ui_state, :root_path, "")
    if current_root != identity.root_path || !haskey(ui_state, :scan_hierarchy)
        _begin_scan!(ui_state, identity.root_path, proj, isfile(joinpath(identity.root_path, "device_info.txt")))
    else
        ui_state[:root_path] = identity.root_path
        ui_state[:project] = proj
    end

    scan_id = get(ui_state, :scan_seq, 0) + 1
    ui_state[:scan_seq] = scan_id
    ui_state[:active_scan_id] = scan_id
    ui_state[:scan_path] = identity.root_path
    ui_state[:scan_state] = :cache_reload
    ui_state[:scan_progress] = _new_scan_progress()
    ui_state[:scan_error] = ""
    ui_state[:cache_id] = cache_id
    ui_state[:cache_identity] = identity
    ui_state[:cache_state] = :loading
    ui_state[:cache_error] = ""
    ui_state[:cache_status] = nothing
    ui_state[:cache_errors] = ProjectCacheFileError[]
    delete!(ui_state, :cache_load_progress)
    delete!(ui_state, :source_check_progress)

    events = Channel{NamedTuple}(Inf)
    cancel_token = Base.Threads.Atomic{Bool}(false)
    ui_state[:scan_events] = events
    ui_state[:scan_cancel_token] = cancel_token

    Base.Threads.@spawn begin
        try
            cached_fingerprints, cached_statuses = _cached_cache_index(identity)
            cache_metadata_ref = Ref{Any}(nothing)
            status_ref = Ref{Union{Nothing,ProjectCacheStatus}}(nothing)
            errors = Channel{Any}(Inf)

            @sync begin
                Base.Threads.@spawn try
                    metadata = _stream_project_cache_contents(
                        identity.root_path,
                        proj,
                        cache_id;
                        should_cancel=() -> cancel_token[],
                        on_progress=(progress) -> put!(events, (
                            kind=:progress,
                            scan_id=scan_id,
                            progress=progress,
                        )),
                        on_file_loaded=(measurements) -> put!(events, (
                            kind=:cache_measurements,
                            scan_id=scan_id,
                            measurements=measurements,
                        )),
                    )
                    cache_metadata_ref[] = metadata
                    put!(events, (
                        kind=:cache_loaded,
                        scan_id=scan_id,
                        metadata=metadata,
                    ))
                catch err
                    put!(errors, (err, catch_backtrace()))
                end

                Base.Threads.@spawn try
                    raw = collect_csv_fingerprints(
                        identity.root_path;
                        should_cancel=() -> cancel_token[],
                        on_progress=(progress) -> put!(events, (
                            kind=:progress,
                            scan_id=scan_id,
                            progress=progress,
                        )),
                        phase=:cache_check,
                    )
                    status_ref[] = _cache_status_from_fingerprints(raw, cached_fingerprints, cached_statuses)
                catch err
                    put!(errors, (err, catch_backtrace()))
                end
            end

            if isready(errors)
                err, bt = take!(errors)
                throw(err)
            end
            metadata = cache_metadata_ref[]
            status = status_ref[]
            metadata === nothing && error("Cache reload finished without loading cache contents")
            status === nothing && error("Cache reload finished without checking source files")
            checked_metadata = (
                identity=metadata.identity,
                status=status,
                semantic_fields=metadata.semantic_fields,
                errors=metadata.errors,
                skipped_count=metadata.skipped_count,
                has_device_metadata=metadata.has_device_metadata,
            )
            put!(events, (
                kind=:reload_result,
                scan_id=scan_id,
                metadata=checked_metadata,
                persist=persist,
            ))
        catch err
            if _is_scan_cancel_error(err) || err isa PlotCancelled
                put!(events, (kind=:canceled, scan_id=scan_id))
            elseif err isa ProjectCacheMissingError
                put!(events, (
                    kind=:cache_missing,
                    scan_id=scan_id,
                    identity=identity,
                    error=sprint(showerror, err),
                    persist=persist,
                ))
            else
                put!(events, (kind=:error, scan_id=scan_id, error=err, bt=catch_backtrace()))
            end
        finally
            close(events)
        end
    end
    return nothing
end

function _launch_cache_update_job!(ui_state; full_rebuild::Bool=false)
    identity = get(ui_state, :cache_identity, nothing)
    identity isa ProjectCacheIdentity || error("No project cache is bound to the current project")
    project = _project_by_name(identity.project_name)
    previous_cache_state = get(ui_state, :cache_state, :idle)
    _cancel_figure_script_job!(ui_state)
    _clear_plot_jobs!(ui_state)

    scan_id = get(ui_state, :scan_seq, 0) + 1
    ui_state[:scan_seq] = scan_id
    ui_state[:active_scan_id] = scan_id
    ui_state[:scan_path] = identity.root_path
    ui_state[:scan_state] = :cache_update
    ui_state[:scan_progress] = _new_scan_progress()
    ui_state[:scan_error] = ""
    ui_state[:cache_state] = :updating
    ui_state[:cache_operation] = full_rebuild ? :rebuild :
        previous_cache_state == :missing ? :create : :update
    ui_state[:cache_error] = ""
    delete!(ui_state, :cache_load_progress)
    delete!(ui_state, :source_check_progress)

    events = Channel{NamedTuple}(Inf)
    cancel_token = Base.Threads.Atomic{Bool}(false)
    ui_state[:scan_events] = events
    ui_state[:scan_cancel_token] = cancel_token

    Base.Threads.@spawn begin
        try
            snapshot = build_project_cache!(
                identity.root_path,
                project,
                identity.cache_id;
                full_rebuild,
                should_cancel=() -> cancel_token[],
                on_progress=(progress) -> put!(events, (
                    kind=:progress,
                    scan_id=scan_id,
                    progress=progress,
                )),
            )
            put!(events, (kind=:cache_result, scan_id=scan_id, snapshot=snapshot))
        catch err
            if _is_scan_cancel_error(err) || err isa PlotCancelled
                put!(events, (kind=:canceled, scan_id=scan_id))
            else
                put!(events, (kind=:error, scan_id=scan_id, error=err, bt=catch_backtrace()))
            end
        finally
            close(events)
        end
    end
end

function _queue_cache_update!(ui_state; full_rebuild::Bool=false)
    if _cache_action_blocked(ui_state)
        _request_scan_cancel!(ui_state)
        return
    end
    _launch_cache_update_job!(ui_state; full_rebuild)
end

function _finalize_scan!(ui_state)
    ui_state[:scan_events] = nothing
    ui_state[:scan_cancel_token] = nothing
end

function _finalize_source_scan!(ui_state)
    ui_state[:source_scan_events] = nothing
    ui_state[:source_scan_cancel_token] = nothing
end

function _poll_scan_events!(ui_state)
    events = get(ui_state, :scan_events, nothing)
    events === nothing && return

    while isready(events)
        msg = take!(events)
        msg.scan_id == get(ui_state, :active_scan_id, 0) || continue

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
            ui_state[:scan_progress] = progress_dict
            p.phase == :cache_load && (ui_state[:cache_load_progress] = progress_dict)
            p.phase == :cache_check && (ui_state[:source_check_progress] = progress_dict)
            p.phase in (:scanning, :cache_update) && (ui_state[:skipped_count] = p.skipped_csv)
            if !(get(ui_state, :cache_state, :idle) == :loading && p.phase == :cache_check)
                ui_state[:scan_state] = p.phase
            end
        elseif msg.kind == :cache_result
            _apply_cache_snapshot!(ui_state, msg.snapshot)
            _persist_preferences!(ui_state; path=msg.snapshot.identity.root_path)
            _finalize_scan!(ui_state)
        elseif msg.kind == :cache_measurements
            _append_cache_measurements!(ui_state, msg.measurements)
        elseif msg.kind == :cache_loaded
            _apply_cache_metadata!(ui_state, msg.metadata)
            delete!(ui_state, :cache_load_progress)
            ui_state[:cache_state] = :checking
            ui_state[:scan_state] = :cache_check
        elseif msg.kind == :reload_result
            _apply_cache_metadata!(ui_state, msg.metadata)
            ui_state[:cache_state] = :ready
            ui_state[:scan_state] = :done
            msg.persist && _persist_preferences!(ui_state; path=msg.metadata.identity.root_path)
            _finalize_scan!(ui_state)
        elseif msg.kind == :cache_missing
            identity = msg.identity
            proj = _project_by_name(identity.project_name)
            _begin_scan!(ui_state, identity.root_path, proj, isfile(joinpath(identity.root_path, "device_info.txt")))
            ui_state[:cache_id] = identity.cache_id
            ui_state[:cache_identity] = identity
            ui_state[:cache_state] = :missing
            ui_state[:cache_error] = msg.error
            ui_state[:scan_state] = :cache_missing
            msg.persist && _persist_preferences!(ui_state; path=identity.root_path)
            _finalize_scan!(ui_state)
        elseif msg.kind == :canceled
            ui_state[:scan_state] = :canceled
            ui_state[:cache_state] = :canceled
            _finalize_scan!(ui_state)
        elseif msg.kind == :error
            ui_state[:scan_state] = :error
            ui_state[:scan_error] = sprint(showerror, msg.error, msg.bt)
            ui_state[:cache_state] = :error
            ui_state[:cache_error] = ui_state[:scan_error]
            @error "Scan job failed" exception = (msg.error, msg.bt)
            _finalize_scan!(ui_state)
        end
    end
end

function _poll_source_scan_events!(ui_state)
    events = get(ui_state, :source_scan_events, nothing)
    events === nothing && return

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
        elseif msg.kind == :source_scan_result
            ui_state[:cache_identity] = msg.identity
            ui_state[:cache_id] = msg.identity.cache_id
            ui_state[:cache_status] = msg.status
            ui_state[:source_scan_state] = :done
            msg.persist && _persist_preferences!(ui_state; path=msg.identity.root_path)
            _finalize_source_scan!(ui_state)
        elseif msg.kind == :canceled
            ui_state[:source_scan_state] = :canceled
            _finalize_source_scan!(ui_state)
        elseif msg.kind == :error
            ui_state[:source_scan_state] = :error
            ui_state[:source_scan_error] = sprint(showerror, msg.error, msg.bt)
            @error "Source scan job failed" exception = (msg.error, msg.bt)
            _finalize_source_scan!(ui_state)
        end
    end
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
