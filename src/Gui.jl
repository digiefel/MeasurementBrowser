import GLFW
using GLMakie
import GLMakie.Makie as Makie
import CImGui as ig
import CImGui.CSyntax: @c
import ModernGL as gl
using Printf

using Statistics: mean

include("MakieIntegration.jl")
using .MakieImguiIntegration

using TOML
using NativeFileDialog: pick_folder

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
    occursin(r"^\d{8}_\d{6}$", stripped) ||
        error("Invalid cache id '$stripped'; expected YYYYMMDD_hhmmss")
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
    entry = _recent_project_entry_for_path(ui_state, path)
    cache_id = entry === nothing ? "" : get(entry, "cache_id", "")
    if isempty(cache_id)
        cache_id = new_project_cache_id()
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
    end
    return cache_id
end

function _persist_current_project_preferences!(ui_state)
    current_root = get(ui_state, :root_path, "")
    isempty(current_root) && return
    _persist_preferences!(ui_state; path=current_root)
end

function _project_for_preference(pref::AbstractString)
    pref == "auto" && return something(_default_project[], RUO2_PROJECT)
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
    ui_state[:active_plot_request] = nothing
    ui_state[:pending_plot_request] = nothing
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

function _copy_bad_registry(registry::BadRegistry)
    return BadRegistry(copy(registry.devices), copy(registry.measurements))
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

function _load_bad_registry_for_root!(ui_state, root_path::String)
    if isempty(root_path)
        ui_state[:bad_registry] = nothing
        ui_state[:bad_registry_error] = ""
        return
    end

    try
        ui_state[:bad_registry] = load_bad_registry(root_path)
        ui_state[:bad_registry_error] = ""
    catch err
        if err isa BadRegistryParseError || err isa BadRegistryIOError
            ui_state[:bad_registry] = nothing
            ui_state[:bad_registry_error] = sprint(showerror, err)
            ui_state[:show_bad] = true
            return
        end
        rethrow()
    end
end

function _bad_registry_ready(ui_state)
    return get(ui_state, :bad_registry, nothing) isa BadRegistry && isempty(get(ui_state, :bad_registry_error, ""))
end

function _bad_registry_or_error(ui_state)
    registry = get(ui_state, :bad_registry, nothing)
    registry isa BadRegistry || error("Bad registry unavailable: $(get(ui_state, :bad_registry_error, ""))")
    return registry
end

function _device_path_key(node::HierarchyNode)
    isempty(node.measurements) && error("Leaf node '$(node.name)' has no measurements")
    return device_path_key(first(node.measurements).device_info)
end

function _device_location(node::HierarchyNode)
    isempty(node.measurements) && error("Leaf node '$(node.name)' has no measurements")
    return copy(first(node.measurements).device_info.location)
end

function _device_is_explicitly_bad(ui_state, device_key::String)
    return device_key in _bad_registry_or_error(ui_state).devices
end

function _measurement_is_explicitly_bad(ui_state, measurement_id::String)
    return measurement_id in _bad_registry_or_error(ui_state).measurements
end

function _measurement_is_bad(ui_state, measurement::MeasurementInfo)
    return _measurement_is_explicitly_bad(ui_state, measurement.id) ||
           _device_is_explicitly_bad(ui_state, device_path_key(measurement.device_info))
end

function _assert_bad_registry_visibility_available(ui_state)
    _bad_registry_ready(ui_state) && return
    error("Cannot hide bad items while bad registry is unavailable: $(get(ui_state, :bad_registry_error, ""))")
end

function _device_is_visible(ui_state, device_key::String)
    get(ui_state, :show_bad, true) && return true
    _assert_bad_registry_visibility_available(ui_state)
    return !_device_is_explicitly_bad(ui_state, device_key)
end

function _measurement_is_visible(ui_state, measurement::MeasurementInfo)
    get(ui_state, :show_bad, true) && return true
    _assert_bad_registry_visibility_available(ui_state)
    return !_measurement_is_bad(ui_state, measurement)
end

function _project_visible_selection(ui_state)
    hierarchy = get(ui_state, :scan_hierarchy, nothing)
    if hierarchy === nothing
        return HierarchyNode[], MeasurementInfo[], String[]
    end

    selected_devices = HierarchyNode[]
    for path_key in get(ui_state, :selected_device_paths, String[])
        node = get(hierarchy.index, device_path_tuple(path_key), nothing)
        node === nothing && continue
        node.kind == :leaf || error("Selected device path '$path_key' does not point to a leaf device")
        !_device_is_visible(ui_state, path_key) && continue
        push!(selected_devices, node)
    end

    visible_device_keys = Set(_device_path_key(node) for node in selected_devices)
    measurement_index = get(ui_state, :measurement_index, Dict{String,MeasurementInfo}())
    selected_measurements = MeasurementInfo[]
    for measurement_id in get(ui_state, :selected_measurement_ids, String[])
        measurement = get(measurement_index, measurement_id, nothing)
        measurement === nothing && continue
        device_path_key(measurement.device_info) in visible_device_keys || continue
        !_measurement_is_visible(ui_state, measurement) && continue
        push!(selected_measurements, measurement)
    end
    selected_path = length(selected_devices) == 1 ? _device_location(selected_devices[1]) : String[]
    return selected_devices, selected_measurements, selected_path
end

function _apply_visible_selection!(ui_state)
    selected_devices, selected_measurements, selected_path = _project_visible_selection(ui_state)

    all_measurements = MeasurementInfo[]
    sizehint!(all_measurements, sum(length(device.measurements) for device in selected_devices; init=0))
    for device in selected_devices
        for measurement in device.measurements
            push!(all_measurements, measurement)
        end
    end

    ui_state[:selected_devices] = selected_devices
    ui_state[:selected_measurements] = selected_measurements
    ui_state[:selected_all_measurements] = all_measurements
    ui_state[:selected_measurement_id_set] = Set(get(ui_state, :selected_measurement_ids, String[]))
    ui_state[:selected_path] = selected_path
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
    proj = get(ui_state, :project, RUO2_PROJECT)
    filter_meas = get(ui_state, :_imgui_text_filter_meas, nothing)
    all_measurements = _selected_measurements(ui_state)
    visible_measurements = filter_meas === nothing ?
        all_measurements :
        _visible_measurements(ui_state, proj, all_measurements, filter_meas)
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

function _set_devices_bad!(ui_state, device_keys::Vector{String}, bad::Bool)
    unique_keys = unique(copy(device_keys))
    isempty(unique_keys) && return false
    _bad_registry_ready(ui_state) || return false

    root_path = get(ui_state, :root_path, "")
    isempty(root_path) && error("Cannot update bad registry without an active project root")

    registry = get(ui_state, :bad_registry, nothing)
    registry isa BadRegistry || error("Bad registry is unavailable for editing")

    updated = _copy_bad_registry(registry)
    for device_key in unique_keys
        bad ? push!(updated.devices, device_key) : delete!(updated.devices, device_key)
    end

    try
        save_bad_registry(root_path, updated)
    catch err
        if err isa BadRegistryIOError
            ui_state[:bad_registry_error] = sprint(showerror, err)
            return false
        end
        rethrow()
    end

    ui_state[:bad_registry] = updated
    ui_state[:bad_registry_error] = ""
    _apply_visible_selection!(ui_state)
    return true
end

function _set_measurements_bad!(ui_state, measurement_ids::Vector{String}, bad::Bool)
    unique_ids = unique(copy(measurement_ids))
    isempty(unique_ids) && return false
    _bad_registry_ready(ui_state) || return false

    root_path = get(ui_state, :root_path, "")
    isempty(root_path) && error("Cannot update bad registry without an active project root")

    registry = get(ui_state, :bad_registry, nothing)
    registry isa BadRegistry || error("Bad registry is unavailable for editing")

    updated = _copy_bad_registry(registry)
    for measurement_id in unique_ids
        bad ? push!(updated.measurements, measurement_id) : delete!(updated.measurements, measurement_id)
    end

    try
        save_bad_registry(root_path, updated)
    catch err
        if err isa BadRegistryIOError
            ui_state[:bad_registry_error] = sprint(showerror, err)
            return false
        end
        rethrow()
    end

    ui_state[:bad_registry] = updated
    ui_state[:bad_registry_error] = ""
    _apply_visible_selection!(ui_state)
    return true
end

function _selection_targets(selected_items::Vector{T}, clicked_item::T) where {T}
    if clicked_item in selected_items
        return selected_items
    end
    return T[clicked_item]
end

function _render_bad_registry_error!(ui_state)
    message = get(ui_state, :bad_registry_error, "")
    isempty(message) && return
    ig.TextColored((1.0, 0.5, 0.5, 1.0), "bad_measurements error")
    if ig.BeginItemTooltip()
        ig.PushTextWrapPos(ig.GetFontSize() * 35.0)
        ig.TextUnformatted(message)
        ig.PopTextWrapPos()
        ig.EndTooltip()
    end
end

function _push_bad_text_style!(bad::Bool)
    bad || return false
    ig.PushStyleColor(ig.ImGuiCol_Text, (0.82, 0.35, 0.35, 1.0))
    return true
end

function _measurement_plot_window_entry(measurement::MeasurementInfo)
    return Dict{Symbol,Any}(
        :target_id => measurement.id,
        :filepath => measurement.filepath,
        :measurement => measurement,
        :title => measurement.clean_title,
        :measurement_kind => measurement.measurement_kind,
        :params => _measurement_parameters(measurement),
    )
end

function _cache_plot_version(ui_state)
    identity = get(ui_state, :cache_identity, nothing)
    identity isa ProjectCacheIdentity ||
        error("Plot request requires an active HDF5 cache identity")
    status = get(ui_state, :cache_status, nothing)
    status isa ProjectCacheStatus ||
        error("Plot request requires a loaded HDF5 cache status")
    return (
        identity.cache_id,
        identity.cache_path,
        status.fresh_files,
        status.stale_files,
        status.new_files,
        status.deleted_files,
        status.error_files,
    )
end

function _extra_plot_window_request(ui_state, proj, entry::Dict{Symbol,Any})
    filepath = get(entry, :filepath, "")
    isempty(filepath) && error("Extra plot window entry is missing filepath")

    target_id = string(get(entry, :target_id, ""))
    isempty(target_id) && error("Extra plot window entry is missing target_id")

    measurement = get(entry, :measurement, nothing)
    measurement isa MeasurementInfo ||
        error("Extra plot window entry for '$filepath' is missing cached measurement identity")

    measurement_kind = get(entry, :measurement_kind, nothing)
    measurement_kind isa Symbol || error("Extra plot window entry for '$filepath' is missing measurement_kind")

    device_params = get(entry, :params, nothing)
    device_params isa Dict{Symbol,Any} || error("Extra plot window entry for '$filepath' has invalid params")

    return _single_plot_job_request(
        ui_state,
        proj,
        filepath,
        measurement_kind,
        device_params,
        :extra;
        target_id=target_id,
        measurement=measurement,
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
    ui_state[:active_plot_request] = nothing
    ui_state[:pending_plot_request] = nothing
    ui_state[:active_plot_id] = get(ui_state, :active_plot_id, 0) + 1
end

function _request_plot_cancel!(ui_state)
    _plot_running(ui_state) || return
    token = get(ui_state, :plot_cancel_token, nothing)
    token !== nothing && Base.Threads.atomic_xchg!(token, true)
    ui_state[:plot_state] = :canceling
end

function _plot_params_key(params::Dict{Symbol,Any})
    pairs = sort(collect(params); by=x -> String(first(x)))
    return Tuple((k, repr(v)) for (k, v) in pairs)
end

function _single_plot_job_request(
    ui_state,
    proj,
    filepath::String,
    measurement_kind::Symbol,
    device_params::Dict{Symbol,Any},
    target::Symbol;
    target_id::String,
    plot_key=nothing,
    measurement::Union{Nothing,MeasurementInfo}=nothing,
)
    measurement isa MeasurementInfo ||
        error("Single plot request for '$filepath' is missing cached measurement identity")
    debug_plot_mode = get(ui_state, :debug_plot_mode, false)
    cache_identity = get(ui_state, :cache_identity, nothing)
    cache_identity isa ProjectCacheIdentity ||
        error("Single plot request for '$filepath' is missing cache identity")
    cache_version = _cache_plot_version(ui_state)
    job_key = (
        project_name(proj),
        target,
        target_id,
        measurement.id,
        cache_version,
        measurement_kind,
        _plot_params_key(device_params),
        debug_plot_mode,
    )
    return Dict{Symbol,Any}(
        :kind => :single_file,
        :job_key => job_key,
        :plot_key => plot_key,
        :project => proj,
        :filepath => filepath,
        :measurement => measurement,
        :cache_identity => cache_identity,
        :cache_version => cache_version,
        :measurement_kind => measurement_kind,
        :device_params => device_params,
        :debug_plot_mode => debug_plot_mode,
        :target => target,
        :target_id => target_id,
    )
end

function _combined_plot_job_request(
    ui_state,
    proj,
    combined_kind::Symbol,
    measurements::Vector{MeasurementInfo},
    target::Symbol;
    target_id::String,
    plot_key=nothing,
)
    debug_plot_mode = get(ui_state, :debug_plot_mode, false)
    paths = [m.filepath for m in measurements]
    device_params_list = [merge(m.device_info.parameters, m.parameters) for m in measurements]
    job_key = (
        project_name(proj),
        target,
        target_id,
        combined_kind,
        sort(paths),
        debug_plot_mode,
    )
    return Dict{Symbol,Any}(
        :kind => :combined_files,
        :job_key => job_key,
        :plot_key => plot_key,
        :project => proj,
        :combined_kind => combined_kind,
        :paths => paths,
        :device_params_list => device_params_list,
        :debug_plot_mode => debug_plot_mode,
        :target => target,
        :target_id => target_id,
    )
end

function _launch_plot_job!(ui_state, request::Dict{Symbol,Any})
    plot_id = get(ui_state, :plot_seq, 0) + 1
    ui_state[:plot_seq] = plot_id
    ui_state[:active_plot_id] = plot_id
    ui_state[:active_plot_request] = request
    ui_state[:plot_state] = :loading
    ui_state[:plot_error] = ""
    ui_state[:pending_plot_request] = nothing

    events = Channel{NamedTuple}(8)
    ui_state[:plot_events] = events
    cancel_token = Base.Threads.Atomic{Bool}(false)
    ui_state[:plot_cancel_token] = cancel_token

    Base.Threads.@spawn begin
        try
            if request[:kind] == :single_file
                analyzed = _measurement_group_for_cached_plot(
                    request[:cache_identity],
                    request[:measurement],
                )
                cancel_token[] && throw(PlotCancelled())
                put!(events, (kind=:loaded, plot_id=plot_id))
                put!(events, (kind=:analyzed, plot_id=plot_id, analyzed=analyzed))
            else
                error("Combined HDF5 plot requests are not implemented yet")
            end
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

function _queue_plot_job!(ui_state, request::Dict{Symbol,Any})
    request_job_key = get(request, :job_key, nothing)
    if _plot_running(ui_state)
        active_request = get(ui_state, :active_plot_request, nothing)
        pending_request = get(ui_state, :pending_plot_request, nothing)
        if active_request !== nothing && get(active_request, :job_key, nothing) == request_job_key
            return
        end
        if pending_request !== nothing && get(pending_request, :job_key, nothing) == request_job_key
            return
        end
        ui_state[:pending_plot_request] = request
        _request_plot_cancel!(ui_state)
        return
    end
    _launch_plot_job!(ui_state, request)
end

function _finalize_plot_job!(ui_state)
    ui_state[:plot_events] = nothing
    ui_state[:plot_cancel_token] = nothing
    ui_state[:active_plot_request] = nothing
    pending_request = get(ui_state, :pending_plot_request, nothing)
    ui_state[:pending_plot_request] = nothing
    if pending_request !== nothing
        _launch_plot_job!(ui_state, pending_request)
        return
    end
end

function _apply_plot_result!(ui_state, request::Dict{Symbol,Any}, fig)
    target = get(request, :target, :main)
    if target == :main
        _set_main_plot_figure!(ui_state, fig, get(request, :plot_key, nothing))
        return
    end

    target_id = get(request, :target_id, "")
    open_plots = get(ui_state, :open_plot_windows, nothing)
    open_plots === nothing && return

    for entry in open_plots
        entry_target_id = get(entry, :target_id, get(entry, :filepath, ""))
        entry_target_id == target_id || continue
        if fig === nothing
            delete!(entry, :figure)
        else
            entry[:figure] = fig
        end
        entry[:cache_version] = get(request, :cache_version, nothing)
        return
    end
end

function _poll_plot_events!(ui_state)
    events = get(ui_state, :plot_events, nothing)
    events === nothing && return

    while isready(events)
        msg = take!(events)
        msg.plot_id == get(ui_state, :active_plot_id, 0) || continue

        if msg.kind == :loaded
            ui_state[:plot_state] = :analyzing
        elseif msg.kind == :analyzed
            request = get(ui_state, :active_plot_request, nothing)
            request === nothing && continue
            ui_state[:plot_state] = :drawing
            fig = if request[:kind] == :single_file
                draw_plot_for_file(
                    request[:project],
                    request[:measurement_kind],
                    msg.analyzed;
                    DEBUG=request[:debug_plot_mode],
                    device_params=request[:device_params],
                )
            else
                draw_plot_for_files(
                    request[:project],
                    request[:combined_kind],
                    msg.analyzed;
                    DEBUG=request[:debug_plot_mode],
                )
            end
            _apply_plot_result!(ui_state, request, fig)
            ui_state[:plot_state] = :done
            _finalize_plot_job!(ui_state)
        elseif msg.kind == :canceled
            ui_state[:plot_state] = :canceled
            _finalize_plot_job!(ui_state)
        elseif msg.kind == :error
            ui_state[:plot_state] = :error
            ui_state[:plot_error] = sprint(showerror, msg.error, msg.bt)
            @error "Plot job failed" exception = (msg.error, msg.bt)
            _finalize_plot_job!(ui_state)
        end
    end
end

function _plot_target_loading(ui_state, target::Symbol; target_id::String="")
    _plot_running(ui_state) || return false
    active_request = get(ui_state, :active_plot_request, nothing)
    pending_request = get(ui_state, :pending_plot_request, nothing)
    for req in (active_request, pending_request)
        req === nothing && continue
        get(req, :target, :main) == target || continue
        get(req, :target_id, "") == target_id && return true
    end
    return false
end

function _cache_activity_model(ui_state)
    state = get(ui_state, :scan_state, :idle)
    progress = get(ui_state, :scan_progress, _new_scan_progress())
    total = get(progress, :total_csv, 0)
    processed = get(progress, :processed_csv, 0)
    loaded = get(progress, :loaded_measurements, 0)
    skipped = get(progress, :skipped_csv, 0)
    fraction = total > 0 ? Float32(clamp(processed / total, 0, 1)) : 0.0f0

    if state == :counting
        return (
            title="Source: Counting",
            detail="Counting source CSV files before measurement discovery.",
            progress="Found $processed CSV files",
            fraction,
            cancel_label="Cancel Source Scan",
        )
    elseif state == :scanning
        progress_text = total > 0 ?
            @sprintf("Scanned %d/%d source files, loaded %d measurements, skipped %d", processed, total, loaded, skipped) :
            @sprintf("Scanned %d source files, loaded %d measurements, skipped %d", processed, loaded, skipped)
        return (
            title="Source: Scanning",
            detail="Reading source CSV files to build the in-memory measurement list.",
            progress=progress_text,
            fraction,
            cancel_label="Cancel Source Scan",
        )
    elseif state == :cache_discovery
        progress_text = total > 0 ?
            "Checked $processed/$total source CSV files" :
            "Found $processed source CSV files"
        return (
            title="Cache: Preparing Build",
            detail="Finding source CSV files before writing cache entries.",
            progress=progress_text,
            fraction,
            cancel_label="Cancel Cache Build",
        )
    elseif state == :cache_load
        progress_text = total > 0 ?
            @sprintf("Read %d/%d cached files, loaded %d measurements", processed, total, loaded) :
            @sprintf("Loaded %d measurements", loaded)
        return (
            title="Cache: Loading",
            detail="Reading cached measurements from the HDF5 file.",
            progress=progress_text,
            fraction,
            cancel_label="Cancel Reload",
        )
    elseif state == :cache_check
        progress_text = total > 0 ?
            "Checked $processed/$total source CSV files" :
            "Found $processed source CSV files"
        return (
            title="Source: Checking",
            detail="Comparing source CSV fingerprints with the loaded cache.",
            progress=progress_text,
            fraction,
            cancel_label="Cancel Source Check",
        )
    elseif state == :cache_reload
        return (
            title="Cache: Reloading",
            detail="Starting cache reload.",
            progress="Starting...",
            fraction,
            cancel_label="Cancel Reload",
        )
    elseif state == :cache_update
        progress_text = total > 0 ?
            @sprintf("Processed %d/%d source files, cached %d measurements", processed, total, loaded) :
            @sprintf("Cached %d measurements", loaded)
        operation = get(ui_state, :cache_operation, :update)
        title = operation == :create ? "Cache: Creating" :
            operation == :rebuild ? "Cache: Rebuilding" : "Cache: Updating"
        return (
            title,
            detail="Writing source measurements into the HDF5 cache.",
            progress=progress_text,
            fraction,
            cancel_label="Cancel Cache Build",
        )
    elseif state == :cache_missing
        return (
            title="Cache: Missing",
            detail="No HDF5 cache exists for this project.",
            progress="Build required",
            fraction,
            cancel_label="Cancel",
        )
    elseif state == :canceling
        return (
            title="Canceling",
            detail="Waiting for the active background job to stop.",
            progress="Canceling...",
            fraction,
            cancel_label="Cancel",
        )
    elseif state == :canceled
        return (
            title="Canceled",
            detail="The last background job was canceled.",
            progress="Canceled",
            fraction,
            cancel_label="Cancel",
        )
    elseif state == :error
        return (
            title="Error",
            detail=get(ui_state, :scan_error, "Background job failed."),
            progress="Failed",
            fraction,
            cancel_label="Cancel",
        )
    elseif state == :done
        return (
            title="Ready",
            detail="No background cache or source job is running.",
            progress="Complete",
            fraction=1.0f0,
            cancel_label="Cancel",
        )
    end
    return (
        title="Idle",
        detail="No background cache or source job is running.",
        progress="Idle",
        fraction=0.0f0,
        cancel_label="Cancel",
    )
end

function _progress_fraction(progress)
    total = get(progress, :total_csv, 0)
    processed = get(progress, :processed_csv, 0)
    total <= 0 && return 0.0f0
    return Float32(clamp(processed / total, 0, 1))
end

function _indeterminate_progress_fraction(ui_state)
    frame = get(ui_state, :_frame, 0)
    return Float32(0.12 + 0.76 * (0.5 + 0.5 * sin(frame / 18)))
end

function _render_progress_indicator!(ui_state, item)
    fraction = get(item, :indeterminate, false) ?
        _indeterminate_progress_fraction(ui_state) :
        item.fraction
    overlay = get(item, :indeterminate, false) ? "working..." : ""
    ig.ProgressBar(fraction, (-1, 0), overlay)
end

function _source_rescan_progress_model(ui_state)
    source_state = get(ui_state, :source_scan_state, :idle)
    if !(source_state in (:counting, :cache_check, :canceling, :done, :canceled, :error))
        return nothing
    end
    if source_state == :error
        return (
            title="Rescan: Error",
            progress=get(ui_state, :source_scan_error, "Source scan failed."),
            fraction=0.0f0,
            indeterminate=false,
        )
    end

    progress = get(ui_state, :source_scan_progress, _new_scan_progress())
    total = get(progress, :total_csv, 0)
    processed = get(progress, :processed_csv, 0)
    text = if source_state == :canceling
        "Canceling source rescan..."
    elseif source_state == :canceled
        "Source rescan canceled"
    elseif source_state == :done
        total > 0 ?
            "Source rescan complete: checked $processed/$total CSV files" :
            "Source rescan complete"
    else
        total > 0 ?
            "Checking $processed/$total source CSV files" :
            "Finding source CSV files: $processed found"
    end
    title = if source_state == :done
        "Rescan: Complete"
    elseif source_state == :canceled
        "Rescan: Canceled"
    elseif source_state == :canceling
        "Rescan: Canceling"
    else
        "Rescan: Checking Source"
    end
    return (
        title,
        progress=text,
        fraction=source_state == :done ? 1.0f0 : _progress_fraction(progress),
        indeterminate=source_state in (:counting, :cache_check) && total <= 0,
    )
end

function _cache_progress_models(ui_state)
    models = NamedTuple[]
    load_progress = get(ui_state, :cache_load_progress, nothing)
    if load_progress isa Dict
        total = get(load_progress, :total_csv, 0)
        processed = get(load_progress, :processed_csv, 0)
        loaded = get(load_progress, :loaded_measurements, 0)
        text = total > 0 ?
            "Read $processed/$total cached files, loaded $loaded measurements" :
            "Loaded $loaded measurements"
        push!(models, (
            title="Cache: Loading",
            progress=text,
            fraction=_progress_fraction(load_progress),
            indeterminate=total <= 0,
        ))
    end
    state = get(ui_state, :scan_state, :idle)
    if state in (:cache_discovery, :cache_update, :cache_reload)
        activity = _cache_activity_model(ui_state)
        push!(models, (
            title=activity.title,
            progress=activity.progress,
            fraction=activity.fraction,
            indeterminate=get(get(ui_state, :scan_progress, _new_scan_progress()), :total_csv, 0) <= 0,
        ))
    end
    return models
end

function _source_progress_models(ui_state)
    models = NamedTuple[]
    rescan_model = _source_rescan_progress_model(ui_state)
    if rescan_model !== nothing
        push!(models, rescan_model)
    end

    state = get(ui_state, :scan_state, :idle)
    source_check_progress = get(ui_state, :source_check_progress, nothing)
    if source_check_progress isa Dict
        total = get(source_check_progress, :total_csv, 0)
        processed = get(source_check_progress, :processed_csv, 0)
        text = total > 0 ?
            "Checked $processed/$total source CSV files" :
            "Found $processed source CSV files"
        push!(models, (
            title="Source: Checking",
            progress=text,
            fraction=_progress_fraction(source_check_progress),
            indeterminate=total <= 0,
        ))
    elseif state == :cache_check
        check_progress = get(ui_state, :scan_progress, _new_scan_progress())
        total = get(check_progress, :total_csv, 0)
        processed = get(check_progress, :processed_csv, 0)
        text = total > 0 ?
            "Checked $processed/$total source CSV files" :
            "Found $processed source CSV files"
        push!(models, (
            title="Source: Checking",
            progress=text,
            fraction=_progress_fraction(check_progress),
            indeterminate=total <= 0,
        ))
    end
    return models
end

function _cache_toolbar_model(ui_state)
    identity = get(ui_state, :cache_identity, nothing)
    status = get(ui_state, :cache_status, nothing)
    cache_state = get(ui_state, :cache_state, :idle)
    activity = _cache_activity_model(ui_state)

    if identity === nothing
        return (
            label="Cache: No Project",
            color=(0.28, 0.28, 0.28, 1.0),
            detail="Open a project folder before building a cache.",
        )
    elseif cache_state == :updating
        return (
            label=activity.title,
            color=(0.18, 0.42, 0.78, 1.0),
            detail=activity.detail,
        )
    elseif cache_state == :loading
        return (
            label="Cache: Loading",
            color=(0.18, 0.42, 0.78, 1.0),
            detail="Reading cached measurements from the HDF5 file.",
        )
    elseif cache_state == :checking
        return (
            label="Cache: Loaded",
            color=(0.18, 0.42, 0.78, 1.0),
            detail="Loaded cached measurements.",
        )
    elseif cache_state == :missing
        return (
            label="Cache: Missing",
            color=(0.82, 0.48, 0.12, 1.0),
            detail="No HDF5 cache has been built for this project.",
        )
    elseif cache_state == :error
        return (
            label="Cache: Error",
            color=(0.72, 0.18, 0.18, 1.0),
            detail=get(ui_state, :cache_error, "Cache error"),
        )
    elseif cache_state == :canceled
        return (
            label="Cache: Canceled",
            color=(0.45, 0.45, 0.45, 1.0),
            detail="Cache update was canceled.",
        )
    elseif status isa ProjectCacheStatus
        if status.error_files > 0
            return (
                label="Cache: Errors",
                color=(0.72, 0.18, 0.18, 1.0),
                detail="$(status.error_files) cached file(s) failed to transform",
            )
        end
        has_changes = status.stale_files > 0 || status.new_files > 0 || status.deleted_files > 0
        if has_changes
            return (
                label="Cache: Stale",
                color=(0.78, 0.58, 0.12, 1.0),
                detail="$(status.stale_files) stale, $(status.new_files) new, $(status.deleted_files) deleted",
            )
        end
        return (
            label="Cache: Fresh",
            color=(0.20, 0.58, 0.30, 1.0),
            detail="$(status.fresh_files) files cached",
        )
    end

    return (
        label="Cache: Idle",
        color=(0.36, 0.36, 0.36, 1.0),
        detail="No cache loaded.",
    )
end

function _cache_button_hover_color(color)
    return (
        min(color[1] + 0.10, 1.0),
        min(color[2] + 0.10, 1.0),
        min(color[3] + 0.10, 1.0),
        color[4],
    )
end

function _cache_button_active_color(color)
    return (
        max(color[1] - 0.08, 0.0),
        max(color[2] - 0.08, 0.0),
        max(color[3] - 0.08, 0.0),
        color[4],
    )
end

function _render_cache_toolbar_button!(ui_state)
    model = _cache_toolbar_model(ui_state)
    ig.PushStyleColor(ig.ImGuiCol_Button, model.color)
    ig.PushStyleColor(ig.ImGuiCol_ButtonHovered, _cache_button_hover_color(model.color))
    ig.PushStyleColor(ig.ImGuiCol_ButtonActive, _cache_button_active_color(model.color))
    if ig.Button(model.label)
        ig.OpenPopup("cache_toolbar_popup")
    end
    ig.PopStyleColor()
    ig.PopStyleColor()
    ig.PopStyleColor()
    if ig.BeginItemTooltip()
        ig.TextUnformatted(model.detail)
        ig.EndTooltip()
    end
    _render_cache_toolbar_popup!(ui_state)
end

function _render_cache_toolbar_popup!(ui_state)
    ig.SetNextWindowSize((960, 560), ig.ImGuiCond_Always)
    if ig.BeginPopup("cache_toolbar_popup")
        _render_cache_controls!(ui_state; compact=false)
        ig.EndPopup()
    end
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
    if _scan_running(ui_state)
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
            _load_bad_registry_for_root!(ui_state, msg.metadata.identity.root_path)
            ui_state[:cache_state] = :checking
            ui_state[:scan_state] = :cache_check
        elseif msg.kind == :reload_result
            _apply_cache_metadata!(ui_state, msg.metadata)
            ui_state[:cache_state] = :ready
            ui_state[:scan_state] = :done
            _load_bad_registry_for_root!(ui_state, msg.metadata.identity.root_path)
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
            _load_bad_registry_for_root!(ui_state, identity.root_path)
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

function render_perf_window(ui_state)
    if !get(ui_state, :show_performance_window, false)
        return
    end

    if ig.Begin("Performance")
        raw_io = ig.GetIO()
        fps = unsafe_load(raw_io.Framerate)
        if fps > 0
            ig.Text(
                "FPS: $(round(fps; digits=1))  Frame: " *
                "$(round(1000 / fps; digits=2)) ms"
            )
        end

        if haskey(ui_state, :_gl_info)
            gi = ui_state[:_gl_info]
            for k in (:vendor, :renderer, :version)
                haskey(gi, k) && ig.Text("GL $(k): $(gi[k])")
            end
        end

        if haskey(ui_state, :_node_count)
            ig.Text("Tree nodes rendered: $(ui_state[:_node_count])")
        end
        rows_visible = get(ui_state, :_measurement_rows_visible, 0)
        rows_rendered = get(ui_state, :_measurement_rows_rendered, 0)
        ig.Text("Measurements visible/rendered: $rows_visible / $rows_rendered")

        mem = _memory_snapshot()
        if mem.vmrss_kb !== nothing
            start_rss = get!(ui_state, :_mem_start_rss_kb, mem.vmrss_kb)
            peak_rss = max(get(ui_state, :_mem_peak_rss_kb, mem.vmrss_kb), mem.vmrss_kb)
            ui_state[:_mem_peak_rss_kb] = peak_rss

            read_bytes = mem.read_bytes === nothing ? 0 : mem.read_bytes
            start_read = get!(ui_state, :_mem_start_read_bytes, read_bytes)

            ig.Separator()
            ig.Text("Process memory")
            ig.BulletText(
                "RSS=$( _fmt_kb(mem.vmrss_kb) )  Δstart=$( _fmt_kb(mem.vmrss_kb - start_rss) )  peak=$( _fmt_kb(peak_rss) )",
            )
            ig.BulletText("Anon=$( _fmt_kb(mem.rssanon_kb) )  VSize=$( _fmt_kb(mem.vmsize_kb) )  VPeak=$( _fmt_kb(mem.vmpeak_kb) )")
            ig.BulletText("GC live=$( _fmt_bytes(mem.gc_live_bytes) )  maxrss=$( _fmt_bytes(mem.maxrss_bytes) )")
            ig.BulletText("read_bytes Δstart=$( _fmt_bytes(read_bytes - start_read) )")
        end

        timings = get(ui_state, :_timings,
            Dict{Symbol,Vector{Float64}}())
        allocs = get(ui_state, :_allocs,
            Dict{Symbol,Vector{Int}}())

        for (k, v) in timings
            isempty(v) && continue
            a = get(allocs, k, Int[])
            last_ms = round(v[end]; digits=2)
            mean_ms = round(mean(v); digits=2)
            last_alloc = isempty(a) ? 0.0 : round(a[end] / 1024; digits=1)
            mean_alloc = isempty(a) ? 0.0 : round(mean(a) / 1024; digits=1)
            msg = @sprintf "%s: last=%.2f ms  mean=%.2f ms  alloc_last=%.1f KB  alloc_mean=%.1f KB" String(k) last_ms mean_ms last_alloc mean_alloc

            ig.BulletText(msg)
        end

        profiles = get(ui_state, :figure_script_job_profiles, Any[])
        if ig.CollapsingHeader("Figure-Script Jobs", ig.ImGuiTreeNodeFlags_DefaultOpen)
            if isempty(profiles)
                ig.TextDisabled("No figure-script jobs profiled yet")
            else
                latest = profiles[end]
                latest_profile = latest.profile
                ig.Text(
                    "Last: $(latest.operation) $(repr(latest_profile.group_name))  " *
                    "$(round(latest_profile.total_ms; digits=1)) ms  " *
                    "selected=$(latest_profile.selected_count) / total=$(latest_profile.measurement_count)",
                )
                if !isempty(latest_profile.sections)
                    ig.Text("Sections")
                    for section in latest_profile.sections
                        ig.BulletText(
                            "$(section.key): calls=$(section.calls)  " *
                            "time=$(round(section.duration_ms; digits=2)) ms  " *
                            "alloc=$(round(section.alloc_bytes / 1024; digits=1)) KB",
                        )
                    end
                end
                if !isempty(latest_profile.counters)
                    ig.Text("Counters")
                    for (key, value) in sort!(collect(latest_profile.counters); by=entry -> String(first(entry)))
                        ig.BulletText("$(key): $(value)")
                    end
                end

                if length(profiles) > 1
                    ig.Text("History")
                    for entry in Iterators.reverse(profiles[max(1, end - 5):end-1])
                        profile = entry.profile
                        ig.BulletText(
                            "$(entry.operation) $(repr(profile.group_name)): " *
                            "$(round(profile.total_ms; digits=1)) ms  " *
                            "selected=$(profile.selected_count)",
                        )
                    end
                end
            end
        end

        if ig.Button("Clear timings")
            empty!(get!(() -> Dict{Symbol,Vector{Float64}}(), ui_state, :_timings))
            empty!(get!(() -> Dict{Symbol,Vector{Int}}(), ui_state, :_allocs))
            empty!(get!(ui_state, :figure_script_job_profiles) do
                Any[]
            end)
        end
    end
    ig.End()
end

function render_menu_bar(ui_state)
    if ig.BeginMenuBar()
        if ig.BeginMenu("Project")
            ig.TextDisabled(_project_status_text(ui_state))
            ig.Separator()

            if ig.MenuItem("Open Folder...")
                path = pick_folder()
                if !isnothing(path) && !isempty(path)
                    @info "Selected path: $path"
                    _open_project_path!(ui_state, path)
                end
            end

            recents = get(ui_state, :recent_projects, Dict{String,String}[])
            if ig.BeginMenu("Recent Projects")
                if isempty(recents)
                    ig.TextDisabled("No recent projects")
                else
                    for (idx, entry) in enumerate(recents)
                        path = get(entry, "path", "")
                        isempty(path) && continue
                        pref = get(entry, "project_preference", "auto")
                        label = "$(basename(path)) [$pref]###recent_project_$idx"
                        if ig.MenuItem(label)
                            _open_project_path!(ui_state, path)
                        end
                        if ig.BeginItemTooltip()
                            ig.TextUnformatted(path)
                            ig.EndTooltip()
                        end
                    end
                end
                ig.EndMenu()
            end

            has_root = haskey(ui_state, :root_path) && !isempty(ui_state[:root_path])
            source_rescan_running = _source_scan_running(ui_state)
            rescan_label = source_rescan_running ? "Cancel Rescan" : "Rescan"
            can_rescan = has_root
            !can_rescan && ig.BeginDisabled()
            if ig.MenuItem(rescan_label)
                if source_rescan_running
                    _request_source_scan_cancel!(ui_state)
                else
                    proj = haskey(ui_state, :project) ?
                        ui_state[:project] :
                        _project_for_preference(get(ui_state, :project_preference, "auto"))
                    @info "Rescanning path: $(ui_state[:root_path])"
                    _launch_source_scan_job!(ui_state, ui_state[:root_path], proj)
                end
            end
            !can_rescan && ig.EndDisabled()

            rescan_progress = _source_rescan_progress_model(ui_state)
            if rescan_progress !== nothing
                ig.Spacing()
                ig.TextDisabled(rescan_progress.title)
                ig.TextDisabled(rescan_progress.progress)
                _render_progress_indicator!(ui_state, rescan_progress)
            end

            ig.Separator()
            if ig.MenuItem("Project Settings", C_NULL, get(ui_state, :show_project_window, false))
                ui_state[:show_project_window] = !get(ui_state, :show_project_window, false)
            end
            ig.EndMenu()
        end
        _render_cache_toolbar_button!(ui_state)
        bad_visibility_toggle_enabled = isempty(get(ui_state, :bad_registry_error, "")) || get(ui_state, :show_bad, true)
        !bad_visibility_toggle_enabled && ig.BeginDisabled()
        if ig.MenuItem("Show Bad", C_NULL, get(ui_state, :show_bad, true))
            ui_state[:show_bad] = !get(ui_state, :show_bad, true)
            _apply_visible_selection!(ui_state)
        end
        if !bad_visibility_toggle_enabled
            ig.EndDisabled()
            if ig.BeginItemTooltip()
                ig.TextUnformatted("Fix bad_measurements and rescan before hiding bad items")
                ig.EndTooltip()
            end
        end
        has_measurement_selection = !isempty(get(ui_state, :selected_measurements, MeasurementInfo[]))
        can_open_figure_scripts = has_measurement_selection ||
                                  !isempty(_figure_script_groups(ui_state)) ||
                                  get(ui_state, :show_figure_script_window, false)
        !can_open_figure_scripts && ig.BeginDisabled()
        if ig.MenuItem("Figure Script", C_NULL, get(ui_state, :show_figure_script_window, false))
            ui_state[:show_figure_script_window] = !get(ui_state, :show_figure_script_window, false)
            ui_state[:figure_script_root_path] = get(ui_state, :root_path, "")
        end
        if !can_open_figure_scripts
            ig.EndDisabled()
            if ig.BeginItemTooltip()
                ig.TextUnformatted("Select one or more measurements to build a figure script")
                ig.EndTooltip()
            end
        end
        if ig.BeginMenu("Debug")
            if ig.MenuItem("Performance Window", C_NULL, get(ui_state, :show_performance_window, false))
                ui_state[:show_performance_window] = !get(ui_state, :show_performance_window, false)
            end
            if ig.MenuItem("Debug Plot Mode", C_NULL, get(ui_state, :debug_plot_mode, false))
                ui_state[:debug_plot_mode] = !get(ui_state, :debug_plot_mode, false)
                # Invalidate cached figures when toggled
                _clear_plot_jobs!(ui_state)
                delete!(ui_state, :plot_figure)
                delete!(ui_state, :_last_plot_key)
            end
            ig.EndMenu()
        end
        ig.EndMenuBar()
    end
end

# Left panel (hierarchy tree) rendering
function _tree_node_matches_filter(filter_tree, node::HierarchyNode)
    return ig.ImGuiTextFilter_PassFilter(filter_tree, node.name, C_NULL)
end

function _measurement_filter_text(proj, measurement::MeasurementInfo)
    return string(
        display_label(proj, measurement), "\n",
        measurement.clean_title, "\n",
        kind_label(proj, measurement.measurement_kind),
    )
end

function _measurement_matches_filter(proj, measurement::MeasurementInfo, filter_obj)
    if !ig.ImGuiTextFilter_IsActive(filter_obj)
        return true
    end
    return ig.ImGuiTextFilter_PassFilter(filter_obj, _measurement_filter_text(proj, measurement), C_NULL)
end

function _visible_measurements(ui_state, proj, measurements, filter_meas)
    if get(ui_state, :show_bad, true) && !ig.ImGuiTextFilter_IsActive(filter_meas)
        return measurements
    end

    visible = MeasurementInfo[]
    sizehint!(visible, length(measurements))
    for measurement in measurements
        _measurement_is_visible(ui_state, measurement) || continue
        _measurement_matches_filter(proj, measurement, filter_meas) || continue
        push!(visible, measurement)
    end
    return visible
end

function _render_hierarchy_tree_panel(ui_state, filter_tree)
    ig.BeginChild("Tree", (0, 0), true)
    ig.SeparatorText("Device Selection")
    _render_bad_registry_error!(ui_state)

    root = get(ui_state, :hierarchy_root, nothing)
    meta_keys = get(ui_state, :device_metadata_keys, Symbol[])

    visible_devices = HierarchyNode[]
    visible_device_keys_ref = Ref{Union{Nothing,Vector{String}}}(nothing)
    selected_device_path_set = Set(get(ui_state, :selected_device_paths, String[]))
    all_device_count = 0

    if root !== nothing
        has_visible_leaf_cache = IdDict{HierarchyNode,Bool}()
        subtree_matches_cache = IdDict{HierarchyNode,Bool}()
        device_key_cache = IdDict{HierarchyNode,String}()

        device_key(node::HierarchyNode) = get!(device_key_cache, node) do
            _device_path_key(node)
        end

        function has_visible_leaf(node::HierarchyNode)
            return get!(has_visible_leaf_cache, node) do
                if isempty(children(node))
                    return _device_is_visible(ui_state, device_key(node))
                end
                for child in children(node)
                    has_visible_leaf(child) && return true
                end
                return false
            end
        end

        function subtree_matches(node::HierarchyNode)
            return get!(subtree_matches_cache, node) do
                _tree_node_matches_filter(filter_tree, node) && return true
                for child in children(node)
                    subtree_matches(child) && return true
                end
                return false
            end
        end

        function collect_visible_devices!(node::HierarchyNode, force_show::Bool=false)
            has_visible_leaf(node) || return
            force_show || subtree_matches(node) || return

            direct_match = force_show || _tree_node_matches_filter(filter_tree, node)
            if isempty(children(node))
                push!(visible_devices, node)
                return
            end
            for child in children(node)
                collect_visible_devices!(child, direct_match)
            end
        end

        function count_leaf_nodes(node::HierarchyNode)
            if isempty(children(node))
                return 1
            end
            total = 0
            for child in children(node)
                total += count_leaf_nodes(child)
            end
            return total
        end

        all_device_count = count_leaf_nodes(root)
        for child in children(root)
            collect_visible_devices!(child, false)
        end

        visible_device_keys = function ()
            keys = visible_device_keys_ref[]
            if keys === nothing
                keys = [device_key(node) for node in visible_devices]
                visible_device_keys_ref[] = keys
            end
            return keys
        end

        _render_selection_toolbar!(
            length(get(ui_state, :selected_devices, HierarchyNode[])),
            length(visible_devices),
            all_device_count,
            filter_tree,
            () -> begin
                ui_state[:selected_device_paths] = copy(visible_device_keys())
                _apply_visible_selection!(ui_state)
            end;
            item_label="devices", filter_id="##tree_filter"
        )

        node_seed = UInt64(0x9e3779b97f4a7c15)
        next_node_id(parent_id::UInt64, node::HierarchyNode) = hash(node.name, parent_id)
        to_imgui_id(node_id::UInt64) = Int32(node_id % UInt64(typemax(Int32)))

        function render_node(node::HierarchyNode, node_id::UInt64, force_show::Bool=false)
            has_visible_leaf(node) || return
            force_show || subtree_matches(node) || return

            ui_state[:_node_count] += 1
            ig.TableNextRow()

            direct_match = force_show || _tree_node_matches_filter(filter_tree, node)
            ig.PushID(to_imgui_id(node_id))

            is_leaf = isempty(children(node))
            leaf_device_key = is_leaf ? device_key(node) : nothing
            selected = is_leaf && (leaf_device_key in selected_device_path_set)

            flags = (
                ig.ImGuiTreeNodeFlags_OpenOnArrow |
                ig.ImGuiTreeNodeFlags_OpenOnDoubleClick |
                ig.ImGuiTreeNodeFlags_NavLeftJumpsToParent |
                ig.ImGuiTreeNodeFlags_SpanFullWidth |
                ig.ImGuiTreeNodeFlags_DrawLinesToNodes |
                ig.ImGuiTreeNodeFlags_SpanAllColumns
            )
            if is_leaf
                flags |= (
                    ig.ImGuiTreeNodeFlags_Leaf |
                    ig.ImGuiTreeNodeFlags_Bullet |
                    ig.ImGuiTreeNodeFlags_NoTreePushOnOpen
                )
            end
            if selected
                flags |= ig.ImGuiTreeNodeFlags_Selected
            end
            if ig.ImGuiTextFilter_IsActive(filter_tree) && direct_match && !is_leaf
                flags |= ig.ImGuiTreeNodeFlags_DefaultOpen
            end

            ig.TableSetColumnIndex(0)
            bad_text_pushed = is_leaf && leaf_device_key !== nothing &&
                              _bad_registry_ready(ui_state) &&
                              _push_bad_text_style!(_device_is_explicitly_bad(ui_state, leaf_device_key))
            opened = ig.TreeNodeEx(is_leaf ? "" : node.name, flags, node.name)
            bad_text_pushed && ig.PopStyleColor()

            if ig.IsItemClicked()
                if is_leaf
                    io = ig.GetIO()
                    shift_held = unsafe_load(io.KeyShift)
                    ctrl_held = unsafe_load(io.KeyCtrl)
                    selected_device_paths = copy(get(ui_state, :selected_device_paths, String[]))
                    _update_multi_selection!(selected_device_paths, leaf_device_key, visible_device_keys(), shift_held, ctrl_held)
                    ui_state[:selected_device_paths] = selected_device_paths
                    _apply_visible_selection!(ui_state)
                elseif !isempty(get(ui_state, :selected_device_paths, String[]))
                    ui_state[:selected_device_paths] = String[]
                    _apply_visible_selection!(ui_state)
                end
            end

            if is_leaf && ig.BeginPopupContextItem()
                target_nodes = _selection_targets(get(ui_state, :selected_devices, HierarchyNode[]), node)
                target_keys = [device_key(target) for target in target_nodes]
                selected_count = length(target_keys)
                selected_count > 1 && ig.TextDisabled("Apply to $selected_count devices")

                editable = _bad_registry_ready(ui_state)
                if !editable
                    ig.TextDisabled("Fix bad_measurements and rescan to edit")
                    ig.Separator()
                end

                !editable && ig.BeginDisabled()
                if ig.MenuItem("Mark Bad")
                    _set_devices_bad!(ui_state, target_keys, true)
                end
                if ig.MenuItem("Unmark Bad")
                    _set_devices_bad!(ui_state, target_keys, false)
                end
                !editable && ig.EndDisabled()
                ig.EndPopup()
            end

            dev_meta = nothing
            if is_leaf && !isempty(node.measurements)
                dev_meta = first(node.measurements).device_info.parameters
            end
            for (i, k) in enumerate(meta_keys)
                ig.TableSetColumnIndex(i)
                if dev_meta !== nothing && haskey(dev_meta, k)
                    ig.Text(string(dev_meta[k]))
                elseif is_leaf
                    ig.TextDisabled("--")
                end
            end

            if opened && !is_leaf
                for child in children(node)
                    render_node(child, next_node_id(node_id, child), direct_match)
                end
                ig.TreePop()
            end
            ig.PopID()
        end

        local table_flags = ig.ImGuiTableFlags_BordersV | ig.ImGuiTableFlags_BordersOuterH |
                            ig.ImGuiTableFlags_Resizable | ig.ImGuiTableFlags_RowBg |
                            ig.ImGuiTableFlags_Reorderable | ig.ImGuiTableFlags_Hideable
        ncols = 1 + length(meta_keys) + 1
        if ig.BeginTable("tree_table", ncols, table_flags)
            local index_flags = ig.ImGuiTableColumnFlags_NoHide | ig.ImGuiTableColumnFlags_NoReorder |
                                ig.ImGuiTableColumnFlags_NoSort | ig.ImGuiTableColumnFlags_WidthStretch
            ig.TableSetupColumn("Device", index_flags, 5.0)
            for k in meta_keys
                ig.TableSetupColumn(String(k), ig.ImGuiTableColumnFlags_AngledHeader | ig.ImGuiTableFlags_SizingFixedFit)
            end
            ig.TableSetupColumn("")
            ig.TableAngledHeadersRow()
            ig.TableHeadersRow()
            for child in children(root)
                render_node(child, next_node_id(node_seed, child), false)
            end
            ig.EndTable()
        end
    else
        _render_selection_toolbar!(
            0,
            0,
            0,
            filter_tree,
            () -> nothing;
            item_label="devices", filter_id="##tree_filter"
        )
        ig.Text("No data loaded")
    end

    isempty(visible_devices) && all_device_count > 0 && ig.TextDisabled("No devices match filter")
    ig.EndChild()
end

# Shared multi-select utility functions

"""
    _update_multi_selection!(selected_items, item, all_items, shift_held, ctrl_held)

Common multi-select logic for both device and measurement panels.
Handles range selection (Shift), toggle selection (Ctrl), and single selection.
Modifies selected_items in place.
"""
function _update_multi_selection!(selected_items::Vector{T}, item::T, all_items::Vector{T}, shift_held::Bool, ctrl_held::Bool) where {T}
    if shift_held && !isempty(selected_items)
        # Range selection: select from last selected to current
        last_item = selected_items[end]
        start_idx = findfirst(x -> x == last_item, all_items)
        end_idx = findfirst(x -> x == item, all_items)

        if start_idx !== nothing && end_idx !== nothing
            if start_idx > end_idx
                start_idx, end_idx = end_idx, start_idx
            end
            selected_range = all_items[start_idx:end_idx]
            # Merge with existing selection
            append!(selected_items, selected_range)
            unique!(selected_items)
        end
    elseif ctrl_held
        # Toggle selection
        if item in selected_items
            filter!(x -> x != item, selected_items)
        else
            push!(selected_items, item)
        end
    else
        # Single selection (replace existing)
        empty!(selected_items)
        push!(selected_items, item)
    end
end

"""
    _render_selection_status!(selected_count, filtered_count, total_count, item_type)

Renders consistent selection status display for both device and measurement panels.
Shows "X/Y items selected" with total count if filtering is active.
Uses appropriate colors: green for multi-select, blue for single, disabled for none.
"""
function _render_selection_status!(selected_count::Int, filtered_count::Int, total_count::Int, item_type::String)
    if selected_count > 1
        ig.TextColored((0.2, 0.8, 0.2, 1.0), "$selected_count/$filtered_count $item_type selected")
    elseif selected_count == 1
        ig.TextColored((0.6, 0.8, 1.0, 1.0), "1/$filtered_count $item_type selected")
    else
        ig.TextDisabled("0/$filtered_count $item_type selected")
    end

    if filtered_count != total_count
        ig.SameLine()
        ig.TextDisabled("($total_count total)")
    end
    ig.SameLine()
    _helpmarker("Multi-select: Shift+click=range, Ctrl+click=toggle, Ctrl+A=all")
end

function _render_selection_toolbar!(
    selected_count::Int,
    visible_count::Int,
    total_count::Int,
    filter_obj,
    on_select_all!::Function;
    item_label::String,
    filter_id::String,
)
    _render_selection_status!(selected_count, visible_count, total_count, item_label)

    ig.Text("Filter")
    ig.SameLine()
    _helpmarker("incl,-excl")
    ig.SameLine()
    ig.SetNextItemShortcut(
        ig.ImGuiMod_Ctrl | ig.ImGuiKey_F,
        ig.ImGuiInputFlags_Tooltip
    )
    ig.ImGuiTextFilter_Draw(filter_obj, filter_id, -1)
    if ig.IsKeyPressed(ig.ImGuiKey_A) && ig.IsWindowFocused()
        io = ig.GetIO()
        unsafe_load(io.KeyCtrl) && on_select_all!()
    end
end

function _collect_leaf_nodes!(devices::Vector{HierarchyNode}, node::HierarchyNode)
    if isempty(children(node))
        push!(devices, node)
        return
    end
    for child in children(node)
        _collect_leaf_nodes!(devices, child)
    end
end

function _all_devices(ui_state)
    root = get(ui_state, :hierarchy_root, nothing)
    root === nothing && return HierarchyNode[]

    devices = HierarchyNode[]
    _collect_leaf_nodes!(devices, root)
    return devices
end

function _selected_measurements(ui_state)
    haskey(ui_state, :selected_all_measurements) ||
        error("UI selection state is missing selected_all_measurements")
    return ui_state[:selected_all_measurements]
end

# Right panel (measurements list) rendering
function _render_measurements_panel(ui_state, filter_meas)
    proj = get(ui_state, :project, RUO2_PROJECT)
    ig.BeginChild("Measurements", (0, 0), true)
    ig.SeparatorText("Measurement Selection")

    selected_devices = get(ui_state, :selected_devices, HierarchyNode[])
    all_measurements = _selected_measurements(ui_state)
    visible_measurements = _visible_measurements(ui_state, proj, all_measurements, filter_meas)
    haskey(ui_state, :selected_measurement_id_set) ||
        error("UI selection state is missing selected_measurement_id_set")
    selected_measurement_id_set = ui_state[:selected_measurement_id_set]
    registry_ready = _bad_registry_ready(ui_state)

    _render_selection_toolbar!(
        length(get(ui_state, :selected_measurements, MeasurementInfo[])),
        length(visible_measurements),
        length(all_measurements),
        filter_meas,
        () -> begin
            ui_state[:selected_measurement_ids] = [measurement.id for measurement in visible_measurements]
            _apply_visible_selection!(ui_state)
        end;
        item_label="measurements", filter_id="##measurements_filter"
    )

    if !isempty(selected_devices)
        if length(selected_devices) == 1
            sel_name = join(get(ui_state, :selected_path, [""]), "/")
            ig.Text("Measurements for $sel_name")
            ig.Separator()
        else
            shown = min(3, length(selected_devices))
            first_names = join((selected_devices[i].name for i in 1:shown), ", ")
            ig.Text("Measurements from $(length(selected_devices)) devices: $first_names$(length(selected_devices) > 3 ? "..." : "")")
            ig.Separator()
        end
    end

    if isempty(selected_devices)
        ig.Text("Select one or more devices to view their measurements")
        ig.EndChild()
        return
    end

    ui_state[:_measurement_rows_visible] = length(visible_measurements)
    ui_state[:_measurement_rows_rendered] = 0

    if !isempty(visible_measurements)
        table_flags = ig.ImGuiTableFlags_BordersOuterH | ig.ImGuiTableFlags_BordersOuterV |
                      ig.ImGuiTableFlags_RowBg | ig.ImGuiTableFlags_SizingStretchProp
        if ig.BeginTable("measurements_table", 1, table_flags)
            ig.TableSetupColumn("Measurement", ig.ImGuiTableColumnFlags_WidthStretch)

            rows_rendered = 0
            clipper = ig.lib.ImGuiListClipper()
            try
                ig.Begin(clipper, length(visible_measurements), ig.GetTextLineHeightWithSpacing())
                while ig.Step(clipper)
                    start_idx = Int(unsafe_load(clipper.DisplayStart)) + 1
                    end_idx = Int(unsafe_load(clipper.DisplayEnd))
                    for idx in start_idx:end_idx
                        measurement = visible_measurements[idx]
                        ig.PushID(measurement.id)
                        ig.TableNextRow()
                        ig.TableSetColumnIndex(0)
                        rows_rendered += 1

                        is_selected = measurement.id in selected_measurement_id_set
                        bad_text_pushed = registry_ready &&
                                          _push_bad_text_style!(_measurement_is_bad(ui_state, measurement))
                        if ig.Selectable(display_label(proj, measurement), is_selected, ig.ImGuiSelectableFlags_SpanAllColumns)
                            io = ig.GetIO()
                            shift_held = unsafe_load(io.KeyShift)
                            ctrl_held = unsafe_load(io.KeyCtrl)
                            selected_measurement_ids = copy(get(ui_state, :selected_measurement_ids, String[]))
                            visible_measurement_ids = [visible.id for visible in visible_measurements]
                            _update_multi_selection!(
                                selected_measurement_ids,
                                measurement.id,
                                visible_measurement_ids,
                                shift_held,
                                ctrl_held,
                            )
                            ui_state[:selected_measurement_ids] = selected_measurement_ids
                            _apply_visible_selection!(ui_state)
                        end
                        bad_text_pushed && ig.PopStyleColor()

                        if ig.BeginPopupContextItem()
                            target_measurements = _selection_targets(get(ui_state, :selected_measurements, MeasurementInfo[]), measurement)
                            target_ids = [target.id for target in target_measurements]
                            selected_count = length(target_ids)
                            selected_count > 1 && ig.TextDisabled("Apply to $selected_count measurements")

                            if ig.MenuItem("Open Plot in New Window")
                                open_plots = get!(ui_state, :open_plot_windows) do
                                    Vector{Dict{Symbol,Any}}()
                                end
                                push!(open_plots, _measurement_plot_window_entry(measurement))
                            end

                            editable = _bad_registry_ready(ui_state)
                            if !editable
                                ig.TextDisabled("Fix bad_measurements and rescan to edit")
                                ig.Separator()
                            end

                            !editable && ig.BeginDisabled()
                            if ig.MenuItem("Mark Bad")
                                _set_measurements_bad!(ui_state, target_ids, true)
                            end
                            if ig.MenuItem("Unmark Bad")
                                _set_measurements_bad!(ui_state, target_ids, false)
                            end
                            !editable && ig.EndDisabled()
                            ig.EndPopup()
                        end
                        ig.PopID()
                    end
                end
            finally
                ig.Destroy(clipper)
            end

            ui_state[:_measurement_rows_rendered] = rows_rendered
            ui_state[:_measurement_rows_visible] = length(visible_measurements)
            ig.EndTable()
        end
    else
        ig.TextDisabled("No measurements match filter")
    end
    ig.EndChild()
end

function render_selection_window(ui_state)
    ui_state[:_node_count] = 0
    filter_tree = get!(ui_state, :_imgui_text_filter_tree) do
        ig.ImGuiTextFilter_ImGuiTextFilter(C_NULL)
    end
    filter_meas = get!(ui_state, :_imgui_text_filter_meas) do
        ig.ImGuiTextFilter_ImGuiTextFilter(C_NULL)
    end

    if ig.Begin("Hierarchy", C_NULL, ig.ImGuiWindowFlags_MenuBar)
        render_menu_bar(ui_state)
        ig.Columns(2, "main_layout")
        _time!(ui_state, :device_tree) do
            _render_hierarchy_tree_panel(ui_state, filter_tree)
        end
        ig.NextColumn()
        _time!(ui_state, :measurement_panel) do
            _render_measurements_panel(ui_state, filter_meas)
        end
    end
    ig.End()
end

function _set_main_plot_figure!(ui_state, fig, key)
    if fig === nothing
        delete!(ui_state, :plot_figure)
        delete!(ui_state, :_last_plot_key)
        return
    end
    ui_state[:plot_figure] = fig
    ui_state[:_last_plot_key] = key
end

function _queue_main_plot_request!(ui_state, proj, measurement::MeasurementInfo)
    filepath = measurement.filepath
    cache_version = _cache_plot_version(ui_state)
    plot_key = (measurement.id, cache_version, measurement.parameters)
    last_plot_key = get(ui_state, :_last_plot_key, nothing)
    plot_key == last_plot_key && return

    device_params = merge(measurement.device_info.parameters, measurement.parameters)
    request = _single_plot_job_request(
        ui_state,
        proj,
        filepath,
        measurement.measurement_kind,
        device_params,
        :main;
        target_id="main",
        plot_key=plot_key,
        measurement=measurement,
    )
    _queue_plot_job!(ui_state, request)
end

function _compatible_measurements(proj, measurements::Vector{MeasurementInfo}, combined_type)
    combined_type === nothing && return MeasurementInfo[]
    return filter(m -> m.measurement_kind in compatible_kinds(proj, combined_type), measurements)
end

function _queue_combined_plot_request!(ui_state, proj, measurements::Vector{MeasurementInfo}, combined_type;
                                       target::Symbol, target_id::String, plot_key=nothing)
    compatible = _compatible_measurements(proj, measurements, combined_type)
    length(compatible) < 2 && return compatible
    request = _combined_plot_job_request(
        ui_state,
        proj,
        combined_type,
        compatible,
        target;
        target_id=target_id,
        plot_key=plot_key,
    )
    _queue_plot_job!(ui_state, request)
    return compatible
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
    )
    push!(open_plots, entry)
    _queue_combined_plot_request!(ui_state, proj, measurements, combined_type; target=:extra, target_id)
end

function render_plot_window(ui_state)
    proj = get(ui_state, :project, RUO2_PROJECT)
    selected_measurements = get(ui_state, :selected_measurements, MeasurementInfo[])
    combined_type = get(ui_state, :combined_plot_type, nothing)

    # Check if user requested combined plot generation
    if get(ui_state, :generate_combined_plot, false) &&
       length(selected_measurements) > 1 && combined_type !== nothing

        ui_state[:generate_combined_plot] = false  # Reset flag
        _clear_plot_jobs!(ui_state)

        compatible = _compatible_measurements(proj, selected_measurements, combined_type)
        key = (combined_type, sort([m.filepath for m in compatible]))
        if length(compatible) >= 2
            _queue_combined_plot_request!(ui_state, proj, compatible, combined_type; target=:main, target_id="main", plot_key=key)
        else
            _set_main_plot_figure!(ui_state, nothing, key)
        end
    elseif length(selected_measurements) == 1
        _queue_main_plot_request!(ui_state, proj, selected_measurements[1])
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
            if length(selected_measurements) > 1 && combined_type !== nothing
                compatible = _compatible_measurements(proj, selected_measurements, combined_type)
                if length(compatible) < 2
                    ig.TextColored((1.0, 0.6, 0.2, 1.0), "Not enough compatible measurements for $(combined_type)")
                    ig.Text("Selected: $(length(selected_measurements)), Compatible: $(length(compatible))")
                    if combined_type === :tlm_analysis
                        ig.TextDisabled("TLM Analysis requires ≥2 TLM 4-point measurements")
                    elseif combined_type === :pund_fatigue
                        ig.TextDisabled("PUND Fatigue requires ≥2 PUND measurements")
                    end
                elseif _plot_target_loading(ui_state, :main; target_id="main")
                    ig.TextDisabled("Loading combined plot...")
                else
                    ig.TextColored((1.0, 0.4, 0.4, 1.0), "Combined plot generation failed")
                    ig.Text("Check file formats and data quality")
                end
            elseif length(selected_measurements) > 1 && combined_type === nothing
                ig.TextColored((0.6, 0.8, 1.0, 1.0), "Multiple measurements selected ($(length(selected_measurements)))")
                ig.Text("Choose a combined plot type and click Generate in the Combined Plots window")
            elseif length(selected_measurements) == 0
                ig.TextDisabled("Select measurements from the Hierarchy panel")
                ig.TextDisabled("• Single measurement: regular plot")
                ig.TextDisabled("• Multiple measurements + plot type: combined plot")
            elseif _plot_target_loading(ui_state, :main; target_id="main")
                ig.TextDisabled("Loading plot...")
            else
                ig.TextColored((1.0, 0.4, 0.4, 1.0), "Plot generation failed")
                ig.Text("Check file format and measurement type")
            end
        end

        ig.Separator()
    end
    ig.End()
end



function render_combined_plots_window(ui_state)
    proj = get(ui_state, :project, RUO2_PROJECT)
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




function _render_cache_controls!(ui_state; compact::Bool)
    identity = get(ui_state, :cache_identity, nothing)
    status = get(ui_state, :cache_status, nothing)
    cache_state = get(ui_state, :cache_state, :idle)
    model = _cache_toolbar_model(ui_state)
    running = _scan_running(ui_state)
    activity = _cache_activity_model(ui_state)

    if compact
        ig.Text("Cache")
    else
        ig.TextColored(model.color, model.label)
    end
    ig.TextWrapped(model.detail)
    if running
        cache_progress = _cache_progress_models(ui_state)
        if !isempty(cache_progress)
            for item in cache_progress
                ig.TextDisabled(item.title)
                ig.TextDisabled(item.progress)
                _render_progress_indicator!(ui_state, item)
            end
        elseif get(ui_state, :scan_state, :idle) != :cache_check
            ig.TextDisabled(activity.progress)
            ig.ProgressBar(activity.fraction, (-1, 0))
        end
    end

    if identity isa ProjectCacheIdentity
        ig.Separator()
        ig.Text("Project: $(identity.project_name)")
        ig.TextWrapped("Root: $(identity.root_path)")
        ig.TextWrapped("File: $(identity.cache_path)")
    end

    if status isa ProjectCacheStatus
        ig.Separator()
        ig.Text("Files")
        ig.Text("Total raw CSV: $(status.total_files)")
        ig.Text("Cached: $(status.cached_files)")
        ig.Text("Fresh: $(status.fresh_files)")
        ig.Text("Stale: $(status.stale_files)")
        ig.Text("New: $(status.new_files)")
        ig.Text("Deleted: $(status.deleted_files)")
        status.error_files > 0 && ig.TextColored(
            (1.0, 0.35, 0.35, 1.0),
            "Errors: $(status.error_files)",
        )
    elseif cache_state == :error
        message = get(ui_state, :cache_error, "")
        !isempty(message) && ig.TextWrapped(message)
    end

    cache_errors = get(ui_state, :cache_errors, ProjectCacheFileError[])
    if !isempty(cache_errors)
        ig.Separator()
        ig.TextColored((1.0, 0.35, 0.35, 1.0), "File Errors")
        for (index, file_error) in enumerate(cache_errors)
            index > 20 && break
            ig.TextWrapped(file_error.path)
            ig.TextColored((1.0, 0.55, 0.35, 1.0), file_error.message)
            index < min(length(cache_errors), 20) && ig.Separator()
        end
        length(cache_errors) > 20 && ig.TextDisabled("$(length(cache_errors) - 20) more file errors")
    end

    ig.Separator()
    has_identity = identity isa ProjectCacheIdentity
    has_root = haskey(ui_state, :root_path) && !isempty(get(ui_state, :root_path, ""))

    (!has_identity || running) && ig.BeginDisabled()
    update_label = cache_state == :missing ? "Build Cache" : "Update Cache"
    if ig.Button(update_label, (-1, 0))
        _queue_cache_update!(ui_state)
    end
    (!has_identity || running) && ig.EndDisabled()

    (!has_identity || running) && ig.BeginDisabled()
    if ig.Button("Full Rebuild Cache", (-1, 0))
        _queue_cache_update!(ui_state; full_rebuild=true)
    end
    (!has_identity || running) && ig.EndDisabled()

    (!has_root || running) && ig.BeginDisabled()
    if ig.Button("Reload Cache", (-1, 0))
        _open_project_path!(ui_state, ui_state[:root_path])
    end
    (!has_root || running) && ig.EndDisabled()

    if running && get(ui_state, :scan_state, :idle) != :cache_check
        if ig.Button(activity.cancel_label, (-1, 0))
            _request_scan_cancel!(ui_state)
        end
    end
end

function render_project_window(ui_state)
    get(ui_state, :show_project_window, false) || return

    open_ref = Ref(true)
    if ig.Begin("Project Settings###project_window", open_ref, ig.ImGuiWindowFlags_AlwaysAutoResize)
        # ── Status ──────────────────────────────────────────────────────────
        if haskey(ui_state, :project)
            proj = ui_state[:project]
            n_meas = length(get(ui_state, :all_measurements, MeasurementInfo[]))
            skipped = get(ui_state, :skipped_count, 0)
            ig.TextColored((0.4, 0.8, 1.0, 1.0), "Active: $(project_name(proj))")
            ig.SameLine()
            ig.TextDisabled("— $(project_description(proj))")
            ig.Text("$n_meas measurements loaded")
            if skipped > 0
                ig.TextColored((1.0, 0.6, 0.2, 1.0), "⚠ $skipped CSV file(s) skipped")
                ig.SameLine()
                _helpmarker("These CSV files were skipped because they don't match the active project filter.\nTry selecting a different project below.")
            else
                ig.TextDisabled("0 files skipped")
            end
        else
            ig.TextDisabled("No folder loaded yet")
        end

        ig.Separator()
        ig.Text("Select project:")
        ig.Spacing()

        # ── Radio buttons ────────────────────────────────────────────────────
        pref = get(ui_state, :project_preference, "auto")
        changed = false

        default_project = something(_default_project[], RUO2_PROJECT)
        default_label = "Default ($(project_name(default_project)))"

        if ig.RadioButton(default_label, pref == "auto")
            ui_state[:project_preference] = "auto"
            changed = true
        end
        ig.SameLine()
        _helpmarker("Use the default project without trying to infer the project from the files in the folder.")

        for p in KNOWN_PROJECTS
            pn = project_name(p)
            if ig.RadioButton(pn, pref == pn)
                ui_state[:project_preference] = pn
                changed = true
            end
            ig.SameLine()
            _helpmarker(project_description(p))
        end

        source_progress = _source_progress_models(ui_state)
        if !isempty(source_progress)
            ig.Spacing()
            for item in source_progress
                ig.TextDisabled(item.title)
                ig.TextDisabled(item.progress)
                _render_progress_indicator!(ui_state, item)
            end
            if _source_scan_running(ui_state)
                if ig.Button("Cancel Rescan", (-1, 0))
                    _request_source_scan_cancel!(ui_state)
                end
            elseif _scan_running(ui_state) && get(ui_state, :scan_state, :idle) == :cache_check
                source_activity = _cache_activity_model(ui_state)
                if ig.Button(source_activity.cancel_label, (-1, 0))
                    _request_scan_cancel!(ui_state)
                end
            end
        end

        # ── Persist + rescan on change ───────────────────────────────────────
        if changed
            current_root = get(ui_state, :root_path, "")
            _persist_preferences!(ui_state; path=isempty(current_root) ? nothing : current_root)
            if !isempty(current_root)
                @info "Project preference changed to '$(ui_state[:project_preference])' - reloading cache"
                _open_project_path!(ui_state, current_root)
            end
        end

        ig.Separator()
        _render_cache_controls!(ui_state; compact=true)
    end
    open_ref[] || (ui_state[:show_project_window] = false)
    ig.End()
end

function render_info_window(ui_state)
    proj = get(ui_state, :project, RUO2_PROJECT)
    selected_devices = get(ui_state, :selected_devices, HierarchyNode[])
    selected_measurements = get(ui_state, :selected_measurements, MeasurementInfo[])
    if ig.Begin("Information Panel")
        flags = ig.ImGuiTableFlags_Borders | ig.ImGuiTableFlags_RowBg | ig.ImGuiTableFlags_ScrollY
        ig.BeginTable("info_cols", 2, flags)
        ig.TableSetupColumn("Device")
        ig.TableSetupColumn("Measurement")
        ig.TableHeadersRow()
        ig.TableNextRow()
        ig.TableNextColumn()

        if length(selected_devices) == 1
            meas_vec = selected_devices[1].measurements
            sel_name = join(get(ui_state, :selected_path, [""]), "/")
            ig.Text("Location: $sel_name")
            ig.Separator()
            stats = begin
                try
                    get_measurements_stats(meas_vec, proj)
                catch err
                    @warn "Failed to compute stats" error = err
                    Dict{Symbol,Any}()
                end
            end
            if !isempty(stats)
                ig.Text("Stats")
                ig.BulletText("Total: $(stats[:total_measurements])")
                ig.BulletText("Types: $(join(stats[:measurement_types], ", "))")
                if haskey(stats, :first_measurement)
                    ig.BulletText("First: $(stats[:first_measurement]) ")
                    ig.BulletText("Last:  $(stats[:last_measurement]) ")
                end
            else
                ig.TextDisabled("No stats available")
            end
            ig.Separator()
            # Device-level metadata
            if !isempty(meas_vec)
                dev_meta = first(meas_vec).device_info.parameters
                if !isempty(dev_meta)
                    ig.Text("Device metadata")
                    for (k, v) in dev_meta
                        ig.BulletText("$(k): $(v)")
                    end
                else
                    ig.TextDisabled("No metadata parameters found")
                end
            end
        elseif isempty(selected_devices)
            ig.TextDisabled("Select a device to see details")
        else
            ig.TextDisabled("Select a single device to see details")
        end

        ig.TableNextColumn()
        if length(selected_measurements) == 1
            m = selected_measurements[1]
            ig.Text("Title: $(m.clean_title)")
            ig.Separator()
            ig.BulletText("Type: $(kind_label(proj, m.measurement_kind))")
            ig.BulletText("Timestamp: $(m.timestamp)")
            ig.BulletText("Filename:")
            ig.SameLine()
            ig.TextLinkOpenURL(m.filename, m.filepath)
            ig.Separator()
            if !isempty(m.parameters)
                ig.Text("Parameters")
                for (k, v) in m.parameters
                    ig.BulletText("$(k) = $(v)")
                end
            else
                ig.TextDisabled("No parameters extracted")
            end

            # ---- data-derived statistics (PUND measurements only) ---------
            if m.measurement_kind == :pund
                ig.Separator()
                ig.Text("Statistics")
                cache = get(ui_state, :computed_stats_cache, Dict{Tuple{String,Int},Dict{Symbol,Any}}())
                cache_key = (m.filepath, get(m.parameters, :fatigue_cycle, 0))
                stats = get!(cache, cache_key) do
                    compute_pund_stats(m.filepath, m.parameters, m.device_info.parameters)
                end
                if !isempty(stats)
                    order = [:voltage_max_V, :voltage_baseline_V, :voltage_min_V, :frequency_kHz, :Pr_max_uCcm2]
                    for key in order
                        haskey(stats, key) || continue
                        v = stats[key]
                        label = replace(string(key), "_" => " ")
                        if v === nothing
                            ig.BulletText("$label: —")
                        else
                            ig.BulletText("$label = $v")
                        end
                    end
                else
                    ig.TextDisabled("Could not compute statistics")
                end
            end

        elseif isempty(selected_measurements)
            ig.TextDisabled("Select a measurement to view details")
        else
            ig.TextDisabled("Select a single measurement to view details")
        end
        ig.EndTable()
    end
    ig.End()
end

function _figure_script_output_path(ui_state)
    output_directory = _buffer_string(ui_state[:figure_script_output_dir_buffer])
    script_name = _buffer_string(ui_state[:figure_script_name_buffer])
    try
        return figure_script_path(output_directory, script_name)
    catch err
        err isa FigureScriptValidationError || rethrow()
        return nothing
    end
end

function _write_figure_script_from_ui!(ui_state; overwrite::Bool=false)
    root_path = get(ui_state, :root_path, "")
    isempty(root_path) && throw(FigureScriptValidationError("Open a project folder before generating a figure script"))
    project = get(ui_state, :project, nothing)
    project isa AbstractProject || error("Figure script generation requires an active project")
    path = write_figure_script(
        _buffer_string(ui_state[:figure_script_output_dir_buffer]),
        root_path,
        project,
        _buffer_string(ui_state[:figure_script_name_buffer]),
        copy(_figure_script_groups(ui_state)),
        _current_scan_measurements(ui_state);
        overwrite=overwrite,
    )
    ui_state[:figure_script_overwrite_confirm] = ""
    _set_figure_script_status!(ui_state, "Wrote $(basename(path))")
    return path
end

function _render_figure_script_group_tooltip(proj, preview_measurements::Vector{MeasurementInfo})
    ig.BeginItemTooltip() || return
    for measurement in preview_measurements[1:min(6, end)]
        ig.BulletText("$(display_label(proj, measurement))")
    end
    length(preview_measurements) > 6 && ig.TextDisabled("...")
    ig.EndTooltip()
end

function render_figure_script_window(ui_state)
    get(ui_state, :show_figure_script_window, false) || return

    proj = get(ui_state, :project, RUO2_PROJECT)
    selected_measurements = _selected_measurements_in_panel_order(ui_state)
    selected_count = length(selected_measurements)
    job_running = _figure_script_job_running(ui_state)
    groups = _figure_script_groups(ui_state)
    group_matches = isempty(groups) ? Dict{String,Vector{MeasurementInfo}}() : _ensure_figure_script_group_matches!(ui_state)
    output_path = _figure_script_output_path(ui_state)
    overwrite_path = get(ui_state, :figure_script_overwrite_confirm, "")
    if output_path === nothing || output_path != overwrite_path
        ui_state[:figure_script_overwrite_confirm] = ""
        overwrite_path = ""
    end

    open_ref = Ref(true)
    ig.SetNextWindowSize((520, 430), ig.ImGuiCond_FirstUseEver)
    if ig.Begin("Figure Script", open_ref, ig.ImGuiWindowFlags_NoDocking)
        ig.TextDisabled("$selected_count measurements selected")
        ig.Separator()

        ig.Text("Output Directory")
        output_dir_changed = ig.InputText(
            "##figure_script_output_dir",
            ui_state[:figure_script_output_dir_buffer],
            length(ui_state[:figure_script_output_dir_buffer]),
        )
        output_dir_changed && _persist_current_project_preferences!(ui_state)
        ig.SameLine()
        if ig.Button("Choose...")
            selected_dir = pick_folder()
            if !isnothing(selected_dir) && !isempty(selected_dir)
                _set_buffer_string!(ui_state[:figure_script_output_dir_buffer], _normalize_project_path(selected_dir))
                _persist_current_project_preferences!(ui_state)
            end
        end
        ig.SameLine()
        _helpmarker("Scripts are written to this directory.")

        ig.Spacing()
        ig.Text("Script Name")
        ig.InputText("##figure_script_name", ui_state[:figure_script_name_buffer], length(ui_state[:figure_script_name_buffer]))
        ig.TextDisabled(output_path === nothing ? "<choose output directory and enter script name>" : output_path)

        ig.Spacing()
        ig.Text("Group Name")
        ig.InputText("##figure_group_name", ui_state[:figure_script_group_name_buffer], length(ui_state[:figure_script_group_name_buffer]))

        ig.Spacing()
        ig.Text("Groups")
        ig.SameLine()
        _helpmarker("Each group becomes one entry in the generated script: data[\"group_name\"]::Vector{MeasurementBrowser.FigureMeasurement}")
        if ig.BeginChild("figure_script_groups", (0, 190), true)
            if isempty(groups)
                ig.TextDisabled("No groups yet")
            else
                for (index, group) in enumerate(groups)
                    preview_measurements = get(group_matches, group.name, MeasurementInfo[])
                    match_count = length(preview_measurements)
                    label = "$(group.name) ($(match_count))###figure_group_$index"
                    if ig.Selectable(label, get(ui_state, :figure_script_selected_group, 0) == index)
                        _set_selected_figure_script_group!(ui_state, index)
                    end
                    ig.IsItemHovered() && _render_figure_script_group_tooltip(proj, preview_measurements)
                end
            end
        end
        ig.EndChild()

        group_selected = _selected_figure_script_group(ui_state) !== nothing
        can_apply_selection = selected_count > 0

        (!can_apply_selection || job_running) && ig.BeginDisabled()
        if ig.Button("New Group From Selection")
            _clear_figure_script_messages!(ui_state)
            try
                _create_figure_script_group_from_selection!(ui_state)
            catch err
                if err isa FigureScriptValidationError
                    _set_figure_script_error!(ui_state, err)
                else
                    rethrow()
                end
            end
        end
        (!can_apply_selection || job_running) && ig.EndDisabled()

        (!group_selected || job_running) && ig.BeginDisabled()
        ig.SameLine()
        if ig.Button("Rename Group")
            _clear_figure_script_messages!(ui_state)
            try
                _rename_selected_figure_script_group!(ui_state)
            catch err
                if err isa FigureScriptValidationError
                    _set_figure_script_error!(ui_state, err)
                else
                    rethrow()
                end
            end
        end
        ig.SameLine()
        if ig.Button("Delete Group")
            _clear_figure_script_messages!(ui_state)
            try
                _delete_selected_figure_script_group!(ui_state)
            catch err
                if err isa FigureScriptValidationError
                    _set_figure_script_error!(ui_state, err)
                else
                    rethrow()
                end
            end
        end
        (!group_selected || job_running) && ig.EndDisabled()

        (!group_selected || !can_apply_selection || job_running) && ig.BeginDisabled()
        if group_selected || can_apply_selection
            ig.SameLine()
        end
        if ig.Button("Add Selection")
            _clear_figure_script_messages!(ui_state)
            try
                _add_selection_to_figure_script_group!(ui_state)
            catch err
                if err isa FigureScriptValidationError
                    _set_figure_script_error!(ui_state, err)
                else
                    rethrow()
                end
            end
        end
        ig.SameLine()
        if ig.Button("Remove Selection")
            _clear_figure_script_messages!(ui_state)
            try
                _remove_selection_from_figure_script_group!(ui_state)
            catch err
                if err isa FigureScriptValidationError
                    _set_figure_script_error!(ui_state, err)
                else
                    rethrow()
                end
            end
        end
        (!group_selected || !can_apply_selection || job_running) && ig.EndDisabled()

        ig.Separator()
        job_running && ig.BeginDisabled()
        if ig.Button("Generate Script")
            _clear_figure_script_messages!(ui_state)
            try
                _write_figure_script_from_ui!(ui_state)
            catch err
                if err isa FigureScriptExistsError
                    ui_state[:figure_script_overwrite_confirm] = err.path
                    _set_figure_script_error!(ui_state, err)
                elseif err isa FigureScriptValidationError || err isa FigureScriptIOError
                    _set_figure_script_error!(ui_state, err)
                else
                    rethrow()
                end
            end
        end
        job_running && ig.EndDisabled()

        if !isempty(overwrite_path)
            job_running && ig.BeginDisabled()
            ig.SameLine()
            if ig.Button("Overwrite Existing")
                _clear_figure_script_messages!(ui_state)
                try
                    _write_figure_script_from_ui!(ui_state; overwrite=true)
                catch err
                    if err isa FigureScriptValidationError || err isa FigureScriptIOError
                        _set_figure_script_error!(ui_state, err)
                    else
                        rethrow()
                    end
                end
            end
            ig.SameLine()
            if ig.Button("Cancel Overwrite")
                ui_state[:figure_script_overwrite_confirm] = ""
            end
            job_running && ig.EndDisabled()
        end

        error_message = get(ui_state, :figure_script_error, "")
        !isempty(error_message) && ig.TextColored((1.0, 0.45, 0.45, 1.0), error_message)
        status_message = get(ui_state, :figure_script_status, "")
        !isempty(status_message) && ig.TextColored((0.45, 0.85, 0.55, 1.0), status_message)
        job_running && ig.TextDisabled("Figure-script worker is running...")
    end
    open_ref[] || (ui_state[:show_figure_script_window] = false)
    ig.End()
end

# ------------------------------------------------------------------
# Modal for missing device metadata (shown each scan when missing)
# ------------------------------------------------------------------
function render_device_info_modal(ui_state)
    # Reset dismissal when root path changes
    current_root = get(ui_state, :root_path, "")
    if get(ui_state, :_modal_last_root_path, "") != current_root
        ui_state[:_modal_last_root_path] = current_root
        ui_state[:dev_info_modal] = true
    end
    # always center
    center = ig.ImVec2(0.5, 0.5)
    @c ig.ImGuiViewport_GetCenter(&center, ig.GetMainViewport())
    ig.SetNextWindowPos(center, ig.ImGuiCond_Always, (0.5, 0.5))

    # Show modal if: missing metadata and user hasn't dismissed it this scan
    if get(ui_state, :dev_info_modal, true) && !get(ui_state, :has_device_metadata, true)
        ig.OpenPopup("Device Metadata Missing")
    end

    opened = get(ui_state, :dev_info_modal, true)

    if @c ig.BeginPopupModal("Device Metadata Missing", &opened, ig.ImGuiWindowFlags_AlwaysAutoResize)
        ig.Text("No device metadata file (device_info.txt) was found.")
        ig.Separator()
        ig.TextWrapped("Create a simple text file named device_info.txt in the TOP folder you opened to add extra info (area, thickness, notes, etc.) for each device.")
        ig.Spacing()
        ig.TextWrapped("How to do it:")
        ig.BulletText("Create a new text file: device_info.txt")
        ig.BulletText("First line (header): device_path, area_um2, t_HZO_nm, ...")
        ig.BulletText("Add a column for each property you want to track.")
        ig.BulletText("Add one line per scope. device_path is matched as an exact slash-separated unit sequence.")
        ig.BulletText("Single units match only exact units, never substrings inside another unit.")
        ig.BulletText("Examples: D1, A9, VI/25um, RuO2test/A9/VI/D1, RuO2test_A10/VI/25um/D1")
        ig.BulletText("Longer exact unit sequences override shorter matches when they set the same field.")
        ig.Spacing()
        ig.TextDisabled("Example:")
        ig.TextDisabled("device_path,   area_um2,   t_HZO_nm,   notes,   active")
        ig.TextDisabled("D1,    12.5,   7.0,   all exact D1 units,   true")
        ig.TextDisabled("RuO2test_A10/VI/25um/D1,    12.4,   7.0,   exact device,  true")
        ig.Spacing()
        ig.TextWrapped("Save the file, then rescan or reopen the folder to load these values.")
        ig.Spacing()
        if ig.Button("Got it")
            opened = false
            ig.CloseCurrentPopup()
        end
        ig.EndPopup()
    end
    ui_state[:dev_info_modal] = opened
end

# Render any additional plot windows opened via right-click context menu.
function render_additional_plot_windows(ui_state)
    proj = get(ui_state, :project, RUO2_PROJECT)
    open_plots = get(ui_state, :open_plot_windows, nothing)
    open_plots === nothing && return
    isempty(open_plots) && return
    to_keep = Vector{Dict{Symbol,Any}}()
    for entry in open_plots
        filepath = get(entry, :filepath, "")
        if isempty(filepath) && haskey(entry, :combined_kind)
            title = get(entry, :title, "Combined Plot")
            id = get(entry, :id, "combined_plot")
            target_id = string(get(entry, :target_id, id))
            open_ref = Ref(true)
            if ig.Begin("Plot: $title###plot_window_$id", open_ref)
                if haskey(entry, :figure)
                    f = entry[:figure]
                    _time!(ui_state, :makie_fig) do
                        MakieFigure("measurement_plot_$id", f; auto_resize_x=true, auto_resize_y=true)
                    end
                elseif _plot_target_loading(ui_state, :extra; target_id=target_id)
                    ig.TextDisabled("Loading plot...")
                else
                    ig.Text("No plot available")
                end
                ig.Separator()
                ig.TextDisabled(title)
            end
            ig.End()
            open_ref[] && push!(to_keep, entry)
            continue
        end
        isempty(filepath) && continue
        title = get(entry, :title, basename(filepath))
        target_id = string(get(entry, :target_id, filepath))
        entry[:target_id] = target_id
        # Refresh / create figure (per-window; no global shared Figure)
        cache_version = _cache_plot_version(ui_state)
        existing_cache_version = get(entry, :cache_version, nothing)
        refresh = !haskey(entry, :figure) || existing_cache_version != cache_version
        if refresh
            request = _extra_plot_window_request(ui_state, proj, entry)
            _queue_plot_job!(ui_state, request)
        end
        # Window (allow user to close)
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
            else
                ig.Text("No plot available")
            end
            ig.Separator()
            ig.TextDisabled(basename(filepath))
        end
        ig.End()
        open_ref[] && push!(to_keep, entry)
    end
    ui_state[:open_plot_windows] = to_keep
end

function _setup_docking_layout!(dockspace_id)
    vp = ig.GetMainViewport()
    sz = unsafe_load(vp.Size)

    ig.DockBuilderRemoveNode(dockspace_id)
    ig.DockBuilderAddNode(dockspace_id, Int(ig.ImGuiDockNodeFlags_DockSpace))
    ig.DockBuilderSetNodeSize(dockspace_id, sz)

    # Vertical split: left column 2/5, right column 3/5
    left  = Ref{UInt32}(0)
    right = Ref{UInt32}(0)
    ig.DockBuilderSplitNode(dockspace_id, ig.ImGuiDir_Left, 2/5, left, right)

    # Left column: top 3/4 hierarchy, bottom 1/4 info
    top_left    = Ref{UInt32}(0)
    bottom_left = Ref{UInt32}(0)
    ig.DockBuilderSplitNode(left[], ig.ImGuiDir_Up, 3/4, top_left, bottom_left)

    # Right column: top 3/4 plot, bottom 1/4 combined plots
    top_right    = Ref{UInt32}(0)
    bottom_right = Ref{UInt32}(0)
    ig.DockBuilderSplitNode(right[], ig.ImGuiDir_Up, 3/4, top_right, bottom_right)

    ig.DockBuilderDockWindow("Hierarchy",         top_left[])
    ig.DockBuilderDockWindow("Plot Area",          top_right[])
    ig.DockBuilderDockWindow("Information Panel",  bottom_left[])
    ig.DockBuilderDockWindow("Combined Plots",     bottom_right[])

    ig.DockBuilderFinish(dockspace_id)
end

function create_window_and_run_loop(root_path::Union{Nothing,String}=nothing; engine=nothing, spawn=1)
    ig.set_backend(:GlfwOpenGL3)
    ui_state = Dict{Symbol,Any}()
    _init_scan_state!(ui_state)
    _init_cache_state!(ui_state)
    _init_bad_state!(ui_state)
    _init_figure_script_state!(ui_state)
    _init_plot_state!(ui_state)
    ui_state[:_frame] = 0
    ctx = ig.CreateContext()
    io = ig.GetIO()
    io.ConfigFlags = unsafe_load(io.ConfigFlags) | ig.ImGuiConfigFlags_DockingEnable
    io.ConfigFlags = unsafe_load(io.ConfigFlags) | ig.ImGuiConfigFlags_ViewportsEnable
    io.ConfigFlags = unsafe_load(io.ConfigFlags) | ig.ImGuiConfigFlags_NavEnableKeyboard
    io.IniFilename = Ptr{Cchar}(C_NULL)   # disable layout persistence
    ig.StyleColorsDark()

    # Load persisted preferences
    prefs = _load_prefs()
    ui_state[:project_preference] = get(prefs, "project", "auto")
    ui_state[:recent_projects] = _parse_recent_projects(prefs)

    if root_path !== nothing && root_path != ""
        _open_project_path!(ui_state, root_path)
    end
    first_frame   = Ref(true)
    setup_layout  = Ref(true)
    ig.render(
        ctx;
        engine,
        window_size=(1920, 1080),
        window_title="Measurement Browser",
        opengl_version=v"3.3",
        spawn,
        wait_events=false,
        on_exit=() -> _print_perf_summary(ui_state),
    ) do
        ui_state[:_frame] += 1
        _poll_scan_events!(ui_state)
        _poll_source_scan_events!(ui_state)
        _poll_plot_events!(ui_state)
        _poll_figure_script_job_events!(ui_state)
        if first_frame[] && !haskey(ui_state, :_gl_info)
            ui_state[:_gl_info] = _collect_gl_info!()
            first_frame[] = false
        end
        dockspace_id = ig.DockSpaceOverViewport(0, ig.GetMainViewport())
        if setup_layout[]
            setup_layout[] = false
            _setup_docking_layout!(dockspace_id)
        end
        if !_scan_running(ui_state)
            _time!(ui_state, :plot_warmup) do
                _ensure_plot_runtime_warmed!(ui_state)
            end
        end
        render_selection_window(ui_state)
        render_project_window(ui_state)
        _time!(ui_state, :info) do
            render_info_window(ui_state)
        end
        _time!(ui_state, :figure_scripts) do
            render_figure_script_window(ui_state)
        end
        _time!(ui_state, :plot) do
            render_plot_window(ui_state)
        end
        _time!(ui_state, :extra_plots) do
            render_additional_plot_windows(ui_state)
        end
        _time!(ui_state, :perf_window) do
            render_perf_window(ui_state)
        end
        _time!(ui_state, :combined_plots) do
            render_combined_plots_window(ui_state)
        end
        # Show metadata guidance modal if needed
        render_device_info_modal(ui_state)
    end
end

function start_browser(root_path::Union{Nothing,String}=nothing)
    create_window_and_run_loop(root_path)
end
