using Printf
using Statistics: mean
import GLFW
import CImGui as ig
import CImGui.CSyntax: @c
using NativeFileDialog: pick_folder

# ---------------------------------------------------------------------------
# ImGui ini file — pointer lifetime
# ---------------------------------------------------------------------------

# dear imgui stores the IniFilename pointer by address and does NOT copy it; the buffer must
# outlive the entire process.  We hold it here as a module-level Ref so it is never GC'd.
const _IMGUI_INI_BYTES = Ref{Vector{UInt8}}(UInt8[])

using ..Projects:
    KindProfileRow,
    DEFAULT_PROJECT,
    PROJECTS,
    project_description,
    project_name,
    SourceProfileRow,
    source_label
using ..ItemIndex: SourceScan
using ..Cache: ProjectCacheIdentity
import ..Workspace
using ..Workspace:
    WorkspaceStatus,
    cancel_scan!,
    poll_workspace!,
    rebuild_cache!,
    scan_source!,
    source_scan_running

"""Brighten or darken an ImGui button color by `delta` per RGB channel, clamped to [0, 1]."""
shifted_color(color::NTuple{4,Float64}, delta::Real)::NTuple{4,Float64} =
    (clamp(color[1] + delta, 0, 1), clamp(color[2] + delta, 0, 1),
     clamp(color[3] + delta, 0, 1), color[4])

"""Button color for a [`WorkspaceStatus`](@ref) level."""
function status_color(level::Symbol)::NTuple{4,Float64}
    level === :busy && return (0.18, 0.42, 0.78, 1.0)
    level === :fresh && return (0.20, 0.58, 0.30, 1.0)
    level === :stale && return (0.78, 0.58, 0.12, 1.0)
    level === :missing && return (0.82, 0.48, 0.12, 1.0)
    level === :error && return (0.72, 0.18, 0.18, 1.0)
    return (0.34, 0.34, 0.34, 1.0)
end

"""The workspace status snapshot, or a neutral placeholder when no source is open."""
current_status(state::BrowserState)::WorkspaceStatus =
    state.workspace isa Workspace.Workspace ? state.workspace.status :
    WorkspaceStatus(:none, "No Project", "Open a project folder to build a cache.",
        false, nothing, Pair{String,String}[])

"""Render the cache status button and its control popup, colored by the workspace status."""
function _render_cache_toolbar_button!(state::BrowserState)::Nothing
    status = current_status(state)
    color = status_color(status.level)
    ig.PushStyleColor(ig.ImGuiCol_Button, color)
    ig.PushStyleColor(ig.ImGuiCol_ButtonHovered, shifted_color(color, 0.10))
    ig.PushStyleColor(ig.ImGuiCol_ButtonActive, shifted_color(color, -0.08))
    if ig.Button("Cache: $(status.label)")
        ig.OpenPopup("cache_toolbar_popup")
    end
    ig.PopStyleColor()
    ig.PopStyleColor()
    ig.PopStyleColor()
    if ig.BeginItemTooltip()
        ig.TextUnformatted(status.detail)
        ig.EndTooltip()
    end
    ig.SetNextWindowSize((960, 560), ig.ImGuiCond_Always)
    if ig.BeginPopup("cache_toolbar_popup")
        _render_cache_controls!(state)
        ig.EndPopup()
    end
    return nothing
end

const SCAN_KIND_PROFILE_COLUMNS = String[
    "Kind", "Sources", "Items", "Total", "Detect", "Read", "Entries", "Process", "Stats",
]
const SCAN_SOURCE_PROFILE_COLUMNS = String[
    "Source item", "Kind", "Items", "Total", "Detect", "Read", "Entries", "Process", "Stats",
    "Threads",
]
"""Render the per-kind summary for the latest source scan."""
function _render_scan_kind_table(
    rows::Vector{KindProfileRow},
    state::DataGridState,
)::Nothing
    function cell(row_index::Int, column::Int)::String
        row = rows[row_index]
        column == 1 && return String(row.kind)
        column == 2 && return string(row.source_items)
        column == 3 && return string(row.items)
        seconds = column == 4 ? row.total_seconds :
            column == 5 ? row.detect_seconds :
            column == 6 ? row.read_seconds :
            column == 7 ? row.entries_seconds :
            column == 8 ? row.process_seconds : row.stats_seconds
        return @sprintf("%.1f ms", 1000 * seconds)
    end
    render_data_grid!(
        "scan_kind_profile", state;
        n_rows=length(rows), columns=SCAN_KIND_PROFILE_COLUMNS, cell,
        selection_mode=:cells, height=180.0f0)
    return nothing
end

"""Render source-item timings, sorted by total elapsed time."""
function _render_scan_source_table(
    rows::Vector{SourceProfileRow},
    state::DataGridState,
)::Nothing
    function cell(row_index::Int, column::Int)::String
        row = rows[row_index]
        column == 1 && return row.source_item_label
        column == 2 && return String(row.kind)
        column == 3 && return string(row.items)
        column == 10 && return join(row.thread_ids, ", ")
        seconds = column == 4 ? row.total_seconds :
            column == 5 ? row.detect_seconds :
            column == 6 ? row.read_seconds :
            column == 7 ? row.entries_seconds :
            column == 8 ? row.process_seconds : row.stats_seconds
        return @sprintf("%.2f ms", 1000 * seconds)
    end
    cell_link(row::Int, column::Int)::Union{Nothing,String} =
        column == 1 ? rows[row].source_item_path : nothing
    render_data_grid!(
        "scan_source_profile", state;
        n_rows=length(rows), columns=SCAN_SOURCE_PROFILE_COLUMNS, cell, cell_link,
        selection_mode=:cells, height=260.0f0)
    return nothing
end


"""Render project, cache, annotation, workflow, and debug controls."""
function render_menu_bar(state::BrowserState)::Nothing
    workspace = state.workspace
    if ig.BeginMenuBar()
        if ig.BeginMenu("Project")
            ig.TextDisabled(
                workspace isa Workspace.Workspace ?
                "Active: $(source_label(workspace.source))" :
                "No project loaded",
            )
            ig.Separator()

            if ig.MenuItem("Open Folder...")
                path = pick_folder()
                if !isnothing(path) && !isempty(path)
                    @info "Selected path: $path"
                    _open_project_path!(state, path)
                end
            end

            recents = state.recent_projects
            if ig.BeginMenu("Recent Projects")
                if isempty(recents)
                    ig.TextDisabled("No recent projects")
                else
                    for (idx, entry) in enumerate(recents)
                        path = entry.path
                        isempty(path) && continue
                        pref = entry.project_preference
                        label = "$(basename(path)) [$pref]###recent_project_$idx"
                        if ig.MenuItem(label)
                            _open_project_path!(state, path)
                        end
                        if ig.BeginItemTooltip()
                            ig.TextUnformatted(path)
                            ig.EndTooltip()
                        end
                    end
                end
                ig.EndMenu()
            end

            has_root = workspace isa Workspace.Workspace
            source_rescan_running =
                has_root && source_scan_running(workspace)
            rescan_label = source_rescan_running ? "Cancel Rescan" : "Rescan"
            can_rescan = has_root
            !can_rescan && ig.BeginDisabled()
            if ig.MenuItem(rescan_label)
                if source_rescan_running
                    cancel_scan!(workspace)
                else
                    @info "Rescanning source: $(source_label(workspace.source))"
                    scan_source!(workspace)
                end
            end
            !can_rescan && ig.EndDisabled()

            ig.Separator()
            if ig.MenuItem("Project Settings", C_NULL, state.show_project_window)
                state.show_project_window = !state.show_project_window
            end
            ig.EndMenu()
        end
        _render_cache_toolbar_button!(state)
        show_bad_effective = _show_bad_effective(state)
        bad_visibility_toggle_enabled = _tag_state_ready(state) ||
                                        isempty(state.tag_state_error)
        !bad_visibility_toggle_enabled && ig.BeginDisabled()
        if ig.MenuItem("Show Bad", C_NULL, show_bad_effective)
            state.show_bad = !state.show_bad
        end
        if !bad_visibility_toggle_enabled
            ig.EndDisabled()
            if ig.BeginItemTooltip()
                ig.TextUnformatted("Fix tags.txt and rescan before hiding bad items")
                ig.EndTooltip()
            end
        end
        _render_table_inspector_menu!(state)
        if ig.BeginMenu("Debug")
            if ig.MenuItem(
                "Performance Window",
                C_NULL,
                state.show_performance_window,
            )
                state.show_performance_window = !state.show_performance_window
            end
            ig.EndMenu()
        end
        ig.EndMenuBar()
    end
    return nothing
end

"""Draw the workspace status line, progress, source identity, errors, and scan controls."""
function _render_cache_controls!(state::BrowserState)::Nothing
    workspace = state.workspace
    status = current_status(state)
    ig.TextColored(status_color(status.level), status.label)
    ig.TextWrapped(status.detail)
    status.progress === nothing || ig.ProgressBar(status.progress, (-1, 0))

    workspace isa Workspace.Workspace || return nothing
    identity = workspace.cache.identity
    if identity isa ProjectCacheIdentity
        ig.Separator()
        ig.Text("Source: $(identity.source_label)")
        ig.TextWrapped("Source ID: $(identity.source_id)")
        ig.TextWrapped("File: $(identity.cache_path)")
    end

    if !isempty(status.errors)
        ig.Separator()
        ig.TextColored((1.0, 0.35, 0.35, 1.0), "Source Item Errors")
        ig.TextDisabled("Select a source item to show its items")
        shown = min(length(status.errors), 20)
        for index in 1:shown
            source_item_id, message = status.errors[index]
            ig.PushID(source_item_id)
            if ig.TextLink(basename(source_item_id)) && select_source_item!(state, source_item_id)
                state.tree_filter = ""
                state.item_filter = ""
                state.reset_project_filters = true
            end
            if ig.BeginItemTooltip()
                ig.TextUnformatted(source_item_id)
                ig.EndTooltip()
            end
            ig.TextWrapped(first(split(message, '\n'; limit=2)))
            ig.PopID()
            index < shown && ig.Separator()
        end
        length(status.errors) > shown &&
            ig.TextDisabled("$(length(status.errors) - shown) more file errors")
    end

    ig.Separator()
    scan_running = source_scan_running(workspace)
    scan_label = scan_running ? "Cancel Scan" : "Scan Source"
    if ig.Button(scan_label, (-1, 0))
        scan_running ? cancel_scan!(workspace) : scan_source!(workspace)
    end
    rebuild_disabled = scan_running || !(workspace.index.source isa SourceScan)
    rebuild_disabled && ig.BeginDisabled()
    if ig.Button(status.level === :missing ? "Build Cache" : "Rebuild Cache", (-1, 0))
        rebuild_cache!(workspace)
    end
    rebuild_disabled && ig.EndDisabled()
    return nothing
end

"""Render project selection and project-specific settings."""
function render_project_window(state::BrowserState)::Nothing
    state.show_project_window || return nothing

    open_ref = Ref(true)
    if ig.Begin("Project Settings###project_window", open_ref, ig.ImGuiWindowFlags_AlwaysAutoResize)
        workspace = state.workspace
        if workspace isa Workspace.Workspace
            ig.Text("Active: $(source_label(workspace.source))")
        else
            ig.TextDisabled("No folder loaded yet")
        end

        ig.Separator()

        if state.project_locked
            ig.TextDisabled("Project supplied by the Julia caller")
        else
            pref = state.project_preference
            changed = false

            default_project = something(DEFAULT_PROJECT[])
            default_label = "Default ($(project_name(default_project)))"

            if ig.RadioButton(default_label, pref == "auto")
                state.project_preference = "auto"
                changed = true
            end
            ig.SameLine()
            _helpmarker("Use the default project without trying to infer the project from the files in the folder.")

            for p in PROJECTS
                pn = project_name(p)
                if ig.RadioButton(pn, pref == pn)
                    state.project_preference = pn
                    changed = true
                end
                ig.SameLine()
                _helpmarker(project_description(p))
            end

            if changed
                current_root = (
                    workspace isa Workspace.Workspace &&
                    hasproperty(workspace.source, :root_path)
                ) ? workspace.source.root_path : ""
                _persist_preferences!(
                    state;
                    path=isempty(current_root) ? nothing : current_root,
                )
                if workspace isa Workspace.Workspace
                    @info(
                        "Project preference changed to '$(state.project_preference)' - " *
                        "reloading cache",
                    )
                    _open_project_path!(state, current_root)
                end
            end
        end
    end
    open_ref[] || (state.show_project_window = false)
    ig.End()
    return nothing
end

"""Create the initial hierarchy, information, and plot docking layout."""
function _setup_docking_layout!(dockspace_id::UInt32)::Nothing
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

    ig.DockBuilderDockWindow("Hierarchy",         top_left[])
    ig.DockBuilderDockWindow("Plot Area",          right[])
    ig.DockBuilderDockWindow("Information Panel",  bottom_left[])

    ig.DockBuilderFinish(dockspace_id)
    return nothing
end

"""Return whether the native GLFW window requested shutdown."""
function _window_close_requested()::Bool
    window = ig.current_window()
    window === nothing && return false
    return GLFW.WindowShouldClose(window)
end

"""Create the ImGui context with the docking/viewport configuration the browser needs."""
function _init_browser_context!()
    ig.set_backend(:GlfwOpenGL3)
    ctx = ig.CreateContext()
    io = ig.GetIO()
    io.ConfigFlags = unsafe_load(io.ConfigFlags) | ig.ImGuiConfigFlags_DockingEnable
    io.ConfigFlags = unsafe_load(io.ConfigFlags) | ig.ImGuiConfigFlags_ViewportsEnable
    io.ConfigFlags = unsafe_load(io.ConfigFlags) | ig.ImGuiConfigFlags_NavEnableKeyboard

    # Point imgui at a per-machine ini file so [Table] column widths persist across restarts.
    # The pointer must live for the whole process — _IMGUI_INI_BYTES holds the bytes.
    depot = isempty(DEPOT_PATH) ? homedir() : first(DEPOT_PATH)
    ini_dir = joinpath(depot, "measurementbrowser")
    mkpath(ini_dir)
    ini_path = joinpath(ini_dir, "imgui.ini")
    _IMGUI_INI_BYTES[] = vcat(Vector{UInt8}(codeunits(ini_path)), 0x00)
    io.IniFilename = Ptr{Cchar}(pointer(_IMGUI_INI_BYTES[]))

    ig.StyleColorsDark()
    return ctx
end

"""
On macOS, promote the process to a regular foreground app so its window gets a Dock icon and shows
up in Cmd-Tab. When the renderloop is driven from the REPL the process can otherwise stay an
accessory app with no Dock presence. No-op off macOS; never throws, since a failure here must not
take down the render loop.
"""
function _promote_to_foreground_app()
    Sys.isapple() || return nothing
    try
        objc_class(name) = @ccall objc_getClass(name::Cstring)::Ptr{Cvoid}
        objc_sel(name) = @ccall sel_registerName(name::Cstring)::Ptr{Cvoid}
        nsapplication = objc_class("NSApplication")
        nsapplication == C_NULL && return nothing
        app = @ccall objc_msgSend(
            nsapplication::Ptr{Cvoid}, objc_sel("sharedApplication")::Ptr{Cvoid},
        )::Ptr{Cvoid}
        app == C_NULL && return nothing
        # NSApplicationActivationPolicyRegular = 0
        @ccall objc_msgSend(
            app::Ptr{Cvoid}, objc_sel("setActivationPolicy:")::Ptr{Cvoid}, 0::Clong,
        )::Bool
        @ccall objc_msgSend(
            app::Ptr{Cvoid}, objc_sel("activateIgnoringOtherApps:")::Ptr{Cvoid}, true::Bool,
        )::Cvoid
    catch error
        @warn "Could not promote to a macOS foreground app" exception=error
    end
    return nothing
end

function _set_browser_window_hints(window_start::Symbol)::Nothing
    if window_start == :normal
        GLFW.WindowHint(GLFW.FOCUSED, true)
        GLFW.WindowHint(GLFW.FOCUS_ON_SHOW, true)
    elseif window_start == :background
        GLFW.WindowHint(GLFW.FOCUSED, false)
        GLFW.WindowHint(GLFW.FOCUS_ON_SHOW, false)
    else
        error("Unsupported browser window_start '$window_start'; use :normal or :background")
    end
    return nothing
end

"""Render a small indeterminate loading spinner next to a short label."""
function _render_loading_spinner!(label::AbstractString)::Nothing
    radius = 8.0f0
    thickness = 2.4f0
    size = 24.0f0
    cursor = ig.GetCursorScreenPos()
    draw_list = ig.GetWindowDrawList()
    center = ig.ImVec2(cursor.x + size / 2, cursor.y + size / 2)
    phase = Float32(ig.GetTime() * 8.0)
    spokes = 12
    for index in 0:(spokes - 1)
        angle = phase + Float32(2π * index / spokes)
        alpha = 0.18f0 + 0.82f0 * Float32(index + 1) / Float32(spokes)
        inner = radius * 0.45f0
        outer = radius
        p1 = ig.ImVec2(center.x + cos(angle) * inner, center.y + sin(angle) * inner)
        p2 = ig.ImVec2(center.x + cos(angle) * outer, center.y + sin(angle) * outer)
        ig.AddLine(draw_list, p1, p2, ig.GetColorU32(ig.ImGuiCol_Text, alpha), thickness)
    end
    ig.Dummy((size, size))
    ig.SameLine(0.0f0, 10.0f0)
    ig.TextDisabled(label)
    return nothing
end

"""Render the minimal startup surface shown before expensive first-use GUI work finishes."""
function _render_startup_preparation!()::Nothing
    center = ig.ImVec2(0, 0)
    @c ig.ImGuiViewport_GetCenter(&center, ig.GetMainViewport())
    flags = ig.ImGuiWindowFlags_NoDecoration | ig.ImGuiWindowFlags_NoMove |
            ig.ImGuiWindowFlags_NoSavedSettings | ig.ImGuiWindowFlags_AlwaysAutoResize |
            ig.ImGuiWindowFlags_NoInputs
    ig.SetNextWindowPos(center, ig.ImGuiCond_Always, (0.5, 0.5))
    if ig.Begin("###browser_startup_preparation", C_NULL, flags)
        _render_loading_spinner!("Loading")
    end
    ig.End()
    return nothing
end

"""
Run the browser render loop for a prepared state. With `wait=false` the loop runs as a background
task pinned to thread 1 (required for GLFW) and the call returns that task, leaving the REPL live.
"""
function _run_browser(
    state::BrowserState,
    ctx;
    engine,
    spawn,
    wait::Bool,
    exit_after_frames::Union{Nothing,Int}=nothing,
    window_start::Symbol=:normal,
)
    exit_after_frames === nothing || exit_after_frames > 0 ||
        error("exit_after_frames must be positive when provided; got $exit_after_frames")
    _set_browser_window_hints(window_start)
    first_frame   = Ref(true)
    setup_layout  = Ref(true)
    warmup_started = Ref(false)
    return ig.render(
        ctx;
        engine,
        spawn,
        wait,
        window_size=(1920, 1080),
        window_title="Measurement Browser",
        opengl_version=v"3.3",
        wait_events=false,
        on_exit=() -> begin
            _shutdown_background_jobs!(state)
            _save_project_view_if_changed!(state)
            _print_perf_summary(state)
        end,
    ) do
        if _window_close_requested()
            _shutdown_background_jobs!(state)
            return :imgui_exit_loop
        end
        state.performance.frame += 1
        workspace = state.workspace
        workspace isa Workspace.Workspace && poll_workspace!(workspace)
        if exit_after_frames !== nothing && state.performance.frame >= exit_after_frames
            _shutdown_background_jobs!(state)
            return :imgui_exit_loop
        end
        if first_frame[]
            state.performance.gl_info = _gl_info()
            window_start == :normal && _promote_to_foreground_app()
            first_frame[] = false
        end
        if !state.plots.runtime_warmed
            _render_startup_preparation!()
            if warmup_started[]
                _time!(state, :plot_warmup) do
                    ensure_plot_runtime_warmed!(state)
                end
            else
                warmup_started[] = true
            end
            return nothing
        end
        dockspace_id = ig.DockSpaceOverViewport(0, ig.GetMainViewport())
        if setup_layout[]
            setup_layout[] = false
            _setup_docking_layout!(dockspace_id)
        end
        render_selection_window(state)
        render_project_window(state)
        _time!(state, :info) do
            render_info_window(state)
        end
        _time!(state, :table_inspector) do
            render_table_inspector_window(state)
        end
        _time!(state, :plot) do
            render_plot_window(state)
        end
        _time!(state, :extra_plots) do
            render_additional_plot_windows(state)
        end
        _time!(state, :perf_window) do
            render_perf_window(state)
        end
        _save_project_view_if_changed!(state)
        render_cache_rebuild_modal(state)
        render_collection_parameters_modal(state)
    end
end

"""
Open the browser on an already-opened workspace (see `open_workspace`). With `wait=false` the REPL
stays interactive while the window runs; the returned task completes when the window closes.
"""
function open_browser(
    workspace::Workspace.Workspace;
    engine::Any=nothing,
    spawn::Int=1,
    wait::Bool=true,
    window_start::Symbol=:normal,
)
    prefs = _load_prefs()
    state = BrowserState(
        project_locked=true,
        project_preference=project_name(workspace.project),
        recent_projects=_parse_recent_projects(prefs),
    )
    ctx = _init_browser_context!()
    _attach_workspace!(state, workspace)
    return _run_browser(state, ctx; engine, spawn, wait, window_start=window_start)
end
