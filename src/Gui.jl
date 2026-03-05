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

using DataPlotter
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
    try
        return TOML.parsefile(path)
    catch
        return Dict{String,Any}()
    end
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
        push!(recents, Dict{String,String}(
            "path" => _normalize_project_path(path),
            "project_preference" => _sanitize_project_preference(pref),
        ))
    end

    return recents
end

function _update_recent_projects(recents::Vector{Dict{String,String}}, path::String, pref::String)
    norm_path = _normalize_project_path(path)
    filter!(entry -> get(entry, "path", "") != norm_path, recents)
    pushfirst!(recents, Dict{String,String}(
        "path" => norm_path,
        "project_preference" => pref,
    ))
    length(recents) > _MAX_RECENT_PROJECTS && resize!(recents, _MAX_RECENT_PROJECTS)
    return recents
end

function _persist_preferences!(ui_state; path::Union{Nothing,String}=nothing)
    prefs = _load_prefs()
    pref = _sanitize_project_preference(string(get(ui_state, :project_preference, "auto")))
    ui_state[:project_preference] = pref
    prefs["project"] = pref

    recents = _parse_recent_projects(prefs)
    if path !== nothing && !isempty(path)
        _update_recent_projects(recents, path, pref)
        prefs["recent_projects"] = recents
    end

    _save_prefs(prefs)
    ui_state[:recent_projects] = recents
end

function _project_preference_for_path(ui_state, path::String)
    recents = get(ui_state, :recent_projects, Dict{String,String}[])
    for entry in recents
        if get(entry, "path", "") == path
            pref = get(entry, "project_preference", "auto")
            return _sanitize_project_preference(pref)
        end
    end
    pref = string(get(ui_state, :project_preference, "auto"))
    return _sanitize_project_preference(pref)
end

function _open_project_path!(ui_state, path::String; persist=true)
    norm_path = _normalize_project_path(path)
    ui_state[:project_preference] = _project_preference_for_path(ui_state, norm_path)
    _scan!(ui_state, norm_path; persist_on_success=persist)
end

function _project_status_text(ui_state)
    if haskey(ui_state, :project)
        pname = project_name(ui_state[:project])
        skipped = get(ui_state, :skipped_count, 0)
        return skipped > 0 ? "Active: $pname ($skipped skipped)" : "Active: $pname"
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
    ui_state[:scan_persist_on_success] = false
    ui_state[:pending_scan_path] = nothing
    ui_state[:pending_scan_persist_on_success] = false
end

function _init_plot_state!(ui_state)
    ui_state[:plot_state] = :idle
    ui_state[:plot_error] = ""
    ui_state[:plot_seq] = 0
    ui_state[:active_plot_id] = 0
    ui_state[:plot_events] = nothing
    ui_state[:active_plot_request] = nothing
    ui_state[:pending_plot_request] = nothing
end

function _plot_running(ui_state)
    return get(ui_state, :plot_state, :idle) in (:loading, :building, :canceling)
end

function _clear_plot_jobs!(ui_state)
    ui_state[:plot_state] = :idle
    ui_state[:plot_error] = ""
    ui_state[:plot_events] = nothing
    ui_state[:active_plot_request] = nothing
    ui_state[:pending_plot_request] = nothing
    ui_state[:active_plot_id] = get(ui_state, :active_plot_id, 0) + 1
end

function _request_plot_cancel!(ui_state)
    _plot_running(ui_state) || return
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
    mtime,
    measurement_kind::Symbol,
    device_params::Dict{Symbol,Any},
    target::Symbol;
    target_id::String,
    plot_key=nothing,
)
    debug_plot_mode = get(ui_state, :debug_plot_mode, false)
    job_key = (
        project_name(proj),
        target,
        target_id,
        filepath,
        mtime,
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
        :mtime => mtime,
        :measurement_kind => measurement_kind,
        :device_params => device_params,
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

    Base.Threads.@spawn begin
        try
            loaded = load_plot_input_for_file(
                request[:project],
                request[:filepath],
                request[:measurement_kind];
                device_params=request[:device_params],
            )
            put!(events, (kind=:loaded, plot_id=plot_id, loaded=loaded))
        catch err
            put!(events, (kind=:error, plot_id=plot_id, error=err, bt=catch_backtrace()))
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
        entry[:mtime] = get(request, :mtime, nothing)
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
            request = get(ui_state, :active_plot_request, nothing)
            request === nothing && continue
            ui_state[:plot_state] = :building
            fig = draw_plot_from_input(
                request[:project],
                request[:measurement_kind],
                msg.loaded;
                DEBUG=request[:debug_plot_mode],
                device_params=request[:device_params],
            )
            _apply_plot_result!(ui_state, request, fig)
            ui_state[:plot_state] = :done
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

function _scan_status_summary(ui_state)
    state = get(ui_state, :scan_state, :idle)
    progress = get(ui_state, :scan_progress, _new_scan_progress())
    total = get(progress, :total_csv, 0)
    processed = get(progress, :processed_csv, 0)
    loaded = get(progress, :loaded_measurements, 0)
    skipped = get(progress, :skipped_csv, 0)

    if state == :counting
        return "Counting CSV files... ($processed found)"
    elseif state == :scanning
        if total > 0
            pct = 100 * processed / total
            return @sprintf("Scanning... %d/%d (%.1f%%), loaded %d, skipped %d", processed, total, pct, loaded, skipped)
        end
        return @sprintf("Scanning... %d processed, loaded %d, skipped %d", processed, loaded, skipped)
    elseif state == :canceling
        return "Canceling scan..."
    elseif state == :canceled
        return "Scan canceled"
    elseif state == :error
        return "Scan failed"
    elseif state == :done
        return "Scan complete"
    end
    return "Idle"
end

function _scan_progress_fraction(ui_state)
    progress = get(ui_state, :scan_progress, _new_scan_progress())
    total = get(progress, :total_csv, 0)
    processed = get(progress, :processed_csv, 0)
    total <= 0 && return 0.0f0
    return Float32(clamp(processed / total, 0, 1))
end

function _render_scan_indicator!(ui_state)
    _scan_running(ui_state) || return

    state = get(ui_state, :scan_state, :idle)
    progress = get(ui_state, :scan_progress, _new_scan_progress())
    total = get(progress, :total_csv, 0)
    processed = get(progress, :processed_csv, 0)

    if state == :counting
        ig.TextDisabled("Scan: counting files...")
        return
    elseif state == :canceling
        ig.TextDisabled("Scan: canceling...")
        return
    end

    if total > 0
        ig.TextDisabled(@sprintf("Scan: %d/%d", processed, total))
        ig.SameLine()
        ig.ProgressBar(_scan_progress_fraction(ui_state), (80, 0))
        return
    end

    ig.TextDisabled(@sprintf("Scan: %d", processed))
end

function _render_plot_indicator!(ui_state)
    _plot_running(ui_state) || return

    state = get(ui_state, :plot_state, :idle)
    if state == :loading
        ig.TextDisabled("Plot: loading data...")
    elseif state == :building
        ig.TextDisabled("Plot: building figure...")
    elseif state == :canceling
        ig.TextDisabled("Plot: canceling...")
    end
end

function _scan_running(ui_state)
    return get(ui_state, :scan_state, :idle) in (:counting, :scanning, :canceling)
end

function _request_scan_cancel!(ui_state)
    token = get(ui_state, :scan_cancel_token, nothing)
    token === nothing && return
    Base.Threads.atomic_xchg!(token, true)
    ui_state[:scan_state] = :canceling
end

function _launch_scan_job!(ui_state, path::String; persist_on_success=false)
    scan_id = get(ui_state, :scan_seq, 0) + 1
    ui_state[:scan_seq] = scan_id
    ui_state[:active_scan_id] = scan_id
    ui_state[:scan_path] = path
    ui_state[:scan_persist_on_success] = persist_on_success
    ui_state[:scan_state] = :counting
    ui_state[:scan_progress] = _new_scan_progress()
    ui_state[:scan_error] = ""
    ui_state[:pending_scan_path] = nothing
    ui_state[:pending_scan_persist_on_success] = false

    events = Channel{NamedTuple}(512)
    cancel_token = Base.Threads.Atomic{Bool}(false)
    ui_state[:scan_events] = events
    ui_state[:scan_cancel_token] = cancel_token
    preferred_project = _preferred_project(ui_state)

    Base.Threads.@spawn begin
        try
            hierarchy = scan_directory(
                path;
                project=preferred_project,
                on_progress=(p) -> put!(events, (kind=:progress, scan_id=scan_id, progress=p)),
                should_cancel=() -> cancel_token[],
                count_first=true,
            )
            put!(events, (kind=:result, scan_id=scan_id, path=path, hierarchy=hierarchy))
        catch err
            if err isa ScanCancelled
                put!(events, (kind=:canceled, scan_id=scan_id))
            else
                put!(events, (kind=:error, scan_id=scan_id, error=err, bt=catch_backtrace()))
            end
        finally
            close(events)
        end
    end
end

function _queue_scan!(ui_state, path::String; persist_on_success=false, force_restart=false)
    norm_path = _normalize_project_path(path)
    if _scan_running(ui_state)
        active_path = get(ui_state, :scan_path, "")
        pending_path = get(ui_state, :pending_scan_path, nothing)

        if !force_restart && !isempty(active_path) && norm_path == active_path
            ui_state[:scan_persist_on_success] = get(ui_state, :scan_persist_on_success, false) || persist_on_success
            ui_state[:pending_scan_path] = nothing
            ui_state[:pending_scan_persist_on_success] = false
            return
        end
        if !force_restart && pending_path !== nothing && norm_path == pending_path
            ui_state[:pending_scan_persist_on_success] = get(ui_state, :pending_scan_persist_on_success, false) || persist_on_success
            return
        end

        ui_state[:pending_scan_path] = norm_path
        ui_state[:pending_scan_persist_on_success] = persist_on_success
        _request_scan_cancel!(ui_state)
        return
    end
    _launch_scan_job!(ui_state, norm_path; persist_on_success)
end

function _finalize_scan!(ui_state)
    ui_state[:scan_events] = nothing
    ui_state[:scan_cancel_token] = nothing
    pending_path = get(ui_state, :pending_scan_path, nothing)
    pending_persist = get(ui_state, :pending_scan_persist_on_success, false)
    ui_state[:pending_scan_path] = nothing
    ui_state[:pending_scan_persist_on_success] = false
    if pending_path !== nothing && !isempty(pending_path)
        _launch_scan_job!(ui_state, pending_path; persist_on_success=pending_persist)
    end
end

function _poll_scan_events!(ui_state)
    events = get(ui_state, :scan_events, nothing)
    events === nothing && return

    while isready(events)
        msg = take!(events)
        msg.scan_id == get(ui_state, :active_scan_id, 0) || continue

        if msg.kind == :progress
            p = msg.progress
            ui_state[:scan_progress] = Dict{Symbol,Any}(
                :phase => p.phase,
                :total_csv => p.total_csv,
                :processed_csv => p.processed_csv,
                :loaded_measurements => p.loaded_measurements,
                :skipped_csv => p.skipped_csv,
                :current_path => p.current_path,
            )
            ui_state[:scan_state] = p.phase
        elseif msg.kind == :result
            _apply_scan_result!(ui_state, msg.path, msg.hierarchy)
            if get(ui_state, :scan_persist_on_success, false)
                _persist_preferences!(ui_state; path=msg.path)
            end
            ui_state[:scan_state] = :done
            _finalize_scan!(ui_state)
        elseif msg.kind == :canceled
            ui_state[:scan_state] = :canceled
            _finalize_scan!(ui_state)
        elseif msg.kind == :error
            ui_state[:scan_state] = :error
            ui_state[:scan_error] = sprint(showerror, msg.error, msg.bt)
            @error "Scan job failed" exception = (msg.error, msg.bt)
            _finalize_scan!(ui_state)
        end
    end
end

# Return the preferred AbstractProject (nothing = auto-detect)
function _preferred_project(ui_state)
    pref = get(ui_state, :project_preference, "auto")
    pref == "auto" && return nothing
    for p in KNOWN_PROJECTS
        project_name(p) == pref && return p
    end
    return nothing
end

# ---------------------------------------------------------------------------
# Centralised scan helper (used by menu bar, project window, initial load)
# ---------------------------------------------------------------------------

function _apply_scan_result!(ui_state, path::String, hierarchy)
    _clear_plot_jobs!(ui_state)
    ui_state[:selected_devices] = HierarchyNode[]
    ui_state[:selected_measurements] = MeasurementInfo[]
    delete!(ui_state, :plot_figure)
    delete!(ui_state, :_last_plot_key)

    ui_state[:hierarchy_root] = hierarchy.root
    ui_state[:all_measurements] = hierarchy.all_measurements
    ui_state[:root_path] = path
    ui_state[:has_device_metadata] = hierarchy.has_device_metadata
    ui_state[:project] = hierarchy.project
    ui_state[:skipped_count] = hierarchy.skipped_count
    all_params = Set{Symbol}()
    for m in hierarchy.all_measurements
        for k in keys(m.device_info.parameters)
            push!(all_params, k)
        end
    end
    ui_state[:device_metadata_keys] = sort!(collect(all_params); by=String)
end

function _scan!(ui_state, path::String; persist_on_success=false, force_restart=false)
    _queue_scan!(ui_state, path; persist_on_success, force_restart)
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

        if ig.Button("Clear timings")
            empty!(get!(() -> Dict{Symbol,Vector{Float64}}(), ui_state, :_timings))
            empty!(get!(() -> Dict{Symbol,Vector{Int}}(), ui_state, :_allocs))
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

            if ig.MenuItem("Reload")
                if haskey(ui_state, :root_path) && !isempty(ui_state[:root_path])
                    @info "Reloading path: $(ui_state[:root_path])"
                    _scan!(ui_state, ui_state[:root_path])
                end
            end

            ig.Separator()
            ig.TextDisabled(_scan_status_summary(ui_state))
            if _scan_running(ui_state)
                ig.ProgressBar(_scan_progress_fraction(ui_state), (-1, 0))
                if ig.MenuItem("Cancel Scan")
                    _request_scan_cancel!(ui_state)
                end
            end
            ig.Separator()
            if ig.MenuItem("Project Settings", C_NULL, get(ui_state, :show_project_window, false))
                ui_state[:show_project_window] = !get(ui_state, :show_project_window, false)
            end
            ig.EndMenu()
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
function _render_hierarchy_tree_panel(ui_state, filter_tree)
    ig.BeginChild("Tree", (0, 0), true)
    ig.SeparatorText("Device Selection")
    _render_scan_indicator!(ui_state)

    device_filter = (device, filter_obj) -> ig.ImGuiTextFilter_PassFilter(filter_obj, device.name, C_NULL)
    _render_selection_toolbar!(
        ui_state, filter_tree, :selected_devices, _all_devices(ui_state), device_filter;
        item_label="devices", filter_id="##tree_filter"
    )

    meta_keys = get(ui_state, :device_metadata_keys, Symbol[])

    # Local helpers tied to filter object
    node_matches(node::HierarchyNode) = ig.ImGuiTextFilter_PassFilter(filter_tree, node.name, C_NULL)
    subtree_match(node::HierarchyNode) = node_matches(node) || any(subtree_match(c) for c in children(node))

    function render_node(node::HierarchyNode, path::Vector{String}=String[], force_show::Bool=false)
        # return if neither the node nor any of its descendants match the filter
        force_show || subtree_match(node) || return

        ui_state[:_node_count] += 1
        ig.TableNextRow()
        ig.TableSetColumnIndex(0)

        full_path = vcat(path, node.name)
        direct_match = force_show || node_matches(node)
        unique_id = join(full_path, "/")
        ig.PushID(unique_id)

        # appearance flags
        is_leaf = isempty(children(node))
        selected_devices = get(ui_state, :selected_devices, HierarchyNode[])
        selected = is_leaf && node in selected_devices
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

        opened = ig.TreeNodeEx(is_leaf ? "" : node.name, flags, node.name)
        # Handle device selection with multi-select support
        if ig.IsItemClicked()
            ui_state[:selected_path] = full_path

            if is_leaf
                io = ig.GetIO()
                shift_held = unsafe_load(io.KeyShift)
                ctrl_held = unsafe_load(io.KeyCtrl)
                selected_devices = get!(ui_state, :selected_devices, HierarchyNode[])

                all_devices = _all_devices(ui_state)
                _update_multi_selection!(selected_devices, node, all_devices, shift_held, ctrl_held)
                ui_state[:selected_devices] = selected_devices

                # Note: measurements are selected independently in the measurements panel
            elseif !isempty(get(ui_state, :selected_devices, HierarchyNode[]))
                ui_state[:selected_devices] = HierarchyNode[]
            end
        end

        # Fill metadata columns
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
            else
                # non-leaf left blank
            end
        end

        # render children
        if opened && !is_leaf
            for c in children(node)
                render_node(c, full_path, direct_match)
            end
            ig.TreePop()
        end
        ig.PopID()
    end

    if haskey(ui_state, :hierarchy_root)
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
            for child in children(ui_state[:hierarchy_root])
                render_node(child, String[], false)
            end
            ig.EndTable()
        end
    else
        ig.Text("No data loaded")
    end
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
    _select_all_visible!(ui_state, state_key, all_items, item_filter, filter_obj)

Handles Ctrl+A "select all visible" functionality for both device and measurement panels.
Only selects items that pass the current filter. Updates ui_state[state_key] with filtered selection.
"""
function _select_all_visible!(ui_state, state_key::Symbol, all_items::Vector{T}, item_filter::Function, filter_obj) where {T}
    if ig.IsKeyPressed(ig.ImGuiKey_A) && ig.IsWindowFocused()
        io = ig.GetIO()
        if unsafe_load(io.KeyCtrl)
            selected_items = get!(ui_state, state_key) do
                T[]
            end
            empty!(selected_items)

            for item in all_items
                if item_filter(item, filter_obj)
                    push!(selected_items, item)
                end
            end

            ui_state[state_key] = selected_items
        end
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

"""
    _count_visible_items(items, filter_obj, item_filter)

Counts how many items pass the current filter. Used for consistent status display
across device and measurement panels. item_filter should return true for items that pass.
"""
function _count_visible_items(items::Vector{T}, filter_obj, item_filter::Function) where {T}
    return count(item -> item_filter(item, filter_obj), items)
end

function _render_selection_toolbar!(
    ui_state,
    filter_obj,
    state_key::Symbol,
    all_items::Vector{T},
    item_filter::Function;
    item_label::String,
    filter_id::String,
) where {T}
    selected_items = get!(ui_state, state_key) do
        T[]
    end
    visible_count = _count_visible_items(all_items, filter_obj, item_filter)
    _render_selection_status!(length(selected_items), visible_count, length(all_items), item_label)

    ig.Text("Filter")
    ig.SameLine()
    _helpmarker("incl,-excl")
    ig.SameLine()
    ig.SetNextItemShortcut(
        ig.ImGuiMod_Ctrl | ig.ImGuiKey_F,
        ig.ImGuiInputFlags_Tooltip
    )
    ig.ImGuiTextFilter_Draw(filter_obj, filter_id, -1)
    _select_all_visible!(ui_state, state_key, all_items, item_filter, filter_obj)
end

# Helper function to collect all leaf nodes (devices) from hierarchy
function _collect_leaf_nodes!(devices::Vector{HierarchyNode}, node::HierarchyNode)
    if isempty(children(node))
        # This is a leaf node (device)
        push!(devices, node)
    else
        # Recursively collect from children
        for child in children(node)
            _collect_leaf_nodes!(devices, child)
        end
    end
end

# Helper function to get all devices for range selection
function _all_devices(ui_state)
    if !haskey(ui_state, :hierarchy_root)
        return HierarchyNode[]
    end

    all_devices = HierarchyNode[]
    _collect_leaf_nodes!(all_devices, ui_state[:hierarchy_root])
    return all_devices
end

function _selected_measurements(ui_state)
    selected_devices = get(ui_state, :selected_devices, HierarchyNode[])
    measurements = MeasurementInfo[]
    for device in selected_devices
        append!(measurements, device.measurements)
    end
    return measurements
end

# Right panel (measurements list) rendering
function _render_measurements_panel(ui_state, filter_meas)
    proj = get(ui_state, :project, RUO2_PROJECT)
    ig.BeginChild("Measurements", (0, 0), true)
    ig.SeparatorText("Measurement Selection")

    selected_devices = get(ui_state, :selected_devices, HierarchyNode[])
    all_measurements = _selected_measurements(ui_state)

    measurement_filter_func = (m, filter_obj) -> begin
        if !ig.ImGuiTextFilter_IsActive(filter_obj)
            return true
        end
        # Evaluate the filter against a single combined string so negative tokens work reliably
        filter_text = string(
            display_label(proj, m), "\n",
            m.clean_title, "\n",
            kind_label(proj, m.measurement_kind)
        )
        return ig.ImGuiTextFilter_PassFilter(filter_obj, filter_text, C_NULL)
    end
    _render_selection_toolbar!(
        ui_state, filter_meas, :selected_measurements, all_measurements, measurement_filter_func;
        item_label="measurements", filter_id="##measurements_filter"
    )

    meas_vec = all_measurements

    if !isempty(selected_devices)
        if length(selected_devices) == 1
            sel_name = join(get(ui_state, :selected_path, [""]), "/")
            ig.Text("Measurements for $sel_name")
            ig.Separator()
        else
            device_names = [d.name for d in selected_devices]
            ig.Text("Measurements from $(length(selected_devices)) devices: $(join(device_names[1:min(3, end)], ", "))$(length(device_names) > 3 ? "..." : "")")
            ig.Separator()
        end
    end

    if !isempty(meas_vec)
        any_shown = false
        for m in meas_vec
            passes = measurement_filter_func(m, filter_meas)
            passes || continue
            any_shown = true
            selected_measurements = get!(ui_state, :selected_measurements, MeasurementInfo[])
            is_selected = m in selected_measurements

            if ig.Selectable(display_label(proj, m), is_selected)
                io = ig.GetIO()
                shift_held = unsafe_load(io.KeyShift)
                ctrl_held = unsafe_load(io.KeyCtrl)

                # Use common multi-select logic
                _update_multi_selection!(selected_measurements, m, meas_vec, shift_held, ctrl_held)
                ui_state[:selected_measurements] = selected_measurements
            end
            # Right-click context menu per measurement entry
            if ig.BeginPopupContextItem()
                if ig.MenuItem("Open Plot in New Window")
                    open_plots = get!(ui_state, :open_plot_windows) do
                        Vector{Dict{Symbol,Any}}()
                    end
                    push!(open_plots, Dict(
                        :target_id => m.filepath,
                        :filepath => m.filepath,
                        :title => m.clean_title,
                        :params => deepcopy(m.device_info.parameters),
                    ))
                end
                ig.EndPopup()
            end
        end
        !any_shown && ig.TextDisabled("No measurements match filter")
    else
        ig.Text("Select one or more devices to view their measurements")
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
        _render_hierarchy_tree_panel(ui_state, filter_tree)
        ig.NextColumn()
        _render_measurements_panel(ui_state, filter_meas)
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
    isfile(filepath) || return
    mtime = Dates.unix2datetime(stat(filepath).mtime)
    plot_key = (filepath, mtime, measurement.parameters)
    last_plot_key = get(ui_state, :_last_plot_key, nothing)
    plot_key == last_plot_key && return

    device_params = merge(measurement.device_info.parameters, measurement.parameters)
    request = _single_plot_job_request(
        ui_state,
        proj,
        filepath,
        mtime,
        measurement.measurement_kind,
        device_params,
        :main;
        target_id="main",
        plot_key=plot_key,
    )
    _queue_plot_job!(ui_state, request)
end

function _compatible_measurements(proj, measurements::Vector{MeasurementInfo}, combined_type)
    combined_type === nothing && return MeasurementInfo[]
    return filter(m -> m.measurement_kind in compatible_kinds(proj, combined_type), measurements)
end

function _build_combined_plot_figure(proj, measurements::Vector{MeasurementInfo}, combined_type)
    compatible = _compatible_measurements(proj, measurements, combined_type)
    length(compatible) < 2 && return nothing, compatible

    try
        paths = [m.filepath for m in compatible]
        dev_params = [merge(m.device_info.parameters, m.parameters) for m in compatible]
        fig = figure_for_files(proj, paths, combined_type; device_params_list=dev_params)
        return fig, compatible
    catch err
        @warn "Combined plot generation failed" error = err
        return nothing, compatible
    end
end

function _open_combined_plot_window!(ui_state, fig, combined_type)
    fig === nothing && return
    open_plots = get!(ui_state, :open_plot_windows) do
        Vector{Dict{Symbol,Any}}()
    end
    counter = get!(ui_state, :_combined_plot_counter) do
        0
    end
    ui_state[:_combined_plot_counter] = counter + 1
    push!(open_plots, Dict(
        :figure => fig,
        :title => "Combined: $(string(combined_type))",
        :id => "combined_$(counter + 1)",
    ))
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

        fig, compatible = _build_combined_plot_figure(proj, selected_measurements, combined_type)
        key = (combined_type, sort([m.filepath for m in compatible]))
        _set_main_plot_figure!(ui_state, fig, key)
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
                fig, _ = _build_combined_plot_figure(proj, selected_measurements, current_type)
                _open_combined_plot_window!(ui_state, fig, current_type)
            end
            ig.EndPopup()
        end
        if !can_generate
            ig.EndDisabled()
        end
    end
    ig.End()
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

        if ig.RadioButton("Auto-detect", pref == "auto")
            ui_state[:project_preference] = "auto"
            changed = true
        end
        ig.SameLine()
        _helpmarker("Detect project automatically from the first matching CSV file found in the folder")

        for p in KNOWN_PROJECTS
            pn = project_name(p)
            if ig.RadioButton(pn, pref == pn)
                ui_state[:project_preference] = pn
                changed = true
            end
            ig.SameLine()
            _helpmarker(project_description(p))
        end

        # ── Persist + rescan on change ───────────────────────────────────────
        if changed
            current_root = get(ui_state, :root_path, "")
            _persist_preferences!(ui_state; path=isempty(current_root) ? nothing : current_root)
            if !isempty(current_root)
                @info "Project preference changed to '$(ui_state[:project_preference])' — rescanning"
                _scan!(ui_state, current_root; force_restart=true)
            end
        end

        ig.Separator()
        ig.TextDisabled(_scan_status_summary(ui_state))
        if _scan_running(ui_state)
            ig.ProgressBar(_scan_progress_fraction(ui_state), (-1, 0))
            if ig.Button("Cancel Scan", (-1, 0))
                _request_scan_cancel!(ui_state)
            end
            ig.Spacing()
        elseif get(ui_state, :scan_state, :idle) == :error
            err = get(ui_state, :scan_error, "")
            !isempty(err) && ig.TextWrapped(err)
            ig.Spacing()
        end
        if ig.Button("Rescan current folder", (-1, 0))
            if haskey(ui_state, :root_path) && !isempty(get(ui_state, :root_path, ""))
                _scan!(ui_state, ui_state[:root_path])
            else
                ig.OpenPopup("no_folder_popup")
            end
        end
        if ig.BeginPopup("no_folder_popup")
            ig.Text("No folder loaded. Use Project → Open Folder first.")
            ig.EndPopup()
        end
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

        elseif isempty(selected_measurements)
            ig.TextDisabled("Select a measurement to view details")
        else
            ig.TextDisabled("Select a single measurement to view details")
        end
        ig.EndTable()
    end
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
        ig.BulletText("Add one line per device. device_path can be just a name (A1) or a full path (CHIP1/SITE1/A2)")
        ig.BulletText("A full path entry overrides a simple name entry for the same leaf.")
        ig.Spacing()
        ig.TextDisabled("Example:")
        ig.TextDisabled("device,   area_um2,   t_HZO_nm,   notes,   active")
        ig.TextDisabled("A1,    12.5,   7.0,   baseline,   true")
        ig.TextDisabled("CHIP1/SITE1/A2,    12.4,   7.0,   override,  true")
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
        # Support figure-only entries (no filepath)
        if isempty(filepath) && haskey(entry, :figure)
            title = get(entry, :title, "Combined Plot")
            id = get(entry, :id, "combined_plot")
            open_ref = Ref(true)
            if ig.Begin("Plot: $title###plot_window_$id", open_ref)
                f = entry[:figure]
                _time!(ui_state, :makie_fig) do
                    MakieFigure("measurement_plot_$id", f; auto_resize_x=true, auto_resize_y=true)
                end
                ig.Separator()
                ig.TextDisabled(title)
            end
            ig.End()
            open_ref[] && push!(to_keep, entry)
            continue
        end
        isempty(filepath) && continue
        if !isfile(filepath)
            continue
        end
        title = get(entry, :title, basename(filepath))
        target_id = string(get(entry, :target_id, filepath))
        entry[:target_id] = target_id
        # Refresh / create figure (per-window; no global shared Figure)
        mtime = Dates.unix2datetime(stat(filepath).mtime)
        existing_mtime = get(entry, :mtime, nothing)
        refresh = !haskey(entry, :figure) || existing_mtime != mtime
        if refresh
            k = detect_kind(proj, basename(filepath))
            device_params = get(entry, :params, Dict{Symbol,Any}())
            if !(device_params isa Dict{Symbol,Any})
                device_params = Dict{Symbol,Any}()
            end
            request = _single_plot_job_request(
                ui_state,
                proj,
                filepath,
                mtime,
                k,
                device_params,
                :extra;
                target_id=target_id,
            )
            _queue_plot_job!(ui_state, request)
        end
        # Window (allow user to close)
        open_ref = Ref(true)
        if ig.Begin("Plot: $title###plot_window_$filepath", open_ref)
            if haskey(entry, :figure)
                f = entry[:figure]
                # Sanitize id for ImGui (avoid slashes)
                id_str = replace(filepath, '/' => '_')
                _time!(ui_state, :makie_fig) do
                    MakieFigure("measurement_plot_$id_str", f; auto_resize_x=true, auto_resize_y=true)
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
        _poll_plot_events!(ui_state)
        if first_frame[] && !haskey(ui_state, :_gl_info)
            ui_state[:_gl_info] = _collect_gl_info!()
            first_frame[] = false
        end
        dockspace_id = ig.DockSpaceOverViewport(0, ig.GetMainViewport())
        if setup_layout[]
            setup_layout[] = false
            _setup_docking_layout!(dockspace_id)
        end
        _time!(ui_state, :device_tree) do
            render_selection_window(ui_state)
        end
        render_project_window(ui_state)
        _time!(ui_state, :info) do
            render_info_window(ui_state)
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
