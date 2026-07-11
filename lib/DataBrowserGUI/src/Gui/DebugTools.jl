import CImGui.CSyntax: @c

"""Toggle one boolean Dear ImGui runtime option and return its new value."""
function _debug_io_option!(label::String, option::Ptr{Bool})::Bool
    value = unsafe_load(option)
    if ig.MenuItem(label, C_NULL, value)
        value = !value
        unsafe_store!(option, value)
    end
    return value
end

"""Create the browser's ImPlot context the first time an ImPlot tool is used."""
function _init_implot_context!(state::BrowserState)::Nothing
    state.implot_context == C_NULL || return nothing
    context = ig.lib.ImPlot_CreateContext()
    context == C_NULL && error("ImPlot failed to create its context")
    state.implot_context = context
    return nothing
end

"""Destroy the browser-owned ImPlot context, if it was created."""
function _shutdown_implot_context!(state::BrowserState)::Nothing
    state.implot_context == C_NULL && return nothing
    ig.lib.ImPlot_DestroyContext(state.implot_context)
    state.implot_context = C_NULL
    return nothing
end

"""Render Dear ImGui's runtime diagnostic switches inside the open Debug menu."""
function _render_debug_options_menu!()::Nothing
    ig.BeginMenu("Debug Options") || return nothing
    io = ig.GetIO()

    _debug_io_option!("Native Debugger Attached", io.ConfigDebugIsDebuggerPresent)
    if ig.IsItemHovered(ig.ImGuiHoveredFlags_DelayNormal) && ig.BeginTooltip()
        ig.TextUnformatted("Enable only while LLDB or GDB is attached.")
        ig.EndTooltip()
    end
    ig.Separator()
    _debug_io_option!("Highlight ID Conflicts", io.ConfigDebugHighlightIdConflicts)
    _debug_io_option!(
        "Offer Item Picker for ID Conflicts",
        io.ConfigDebugHighlightIdConflictsShowItemPicker,
    )
    _debug_io_option!("Test Begin Return Once", io.ConfigDebugBeginReturnValueOnce)
    _debug_io_option!("Test Begin Return Loop", io.ConfigDebugBeginReturnValueLoop)
    _debug_io_option!("Ignore Focus Loss", io.ConfigDebugIgnoreFocusLoss)
    _debug_io_option!("Comment .ini Settings", io.ConfigDebugIniSettings)
    ig.Separator()
    _debug_io_option!("Error Recovery", io.ConfigErrorRecovery)
    _debug_io_option!("Recovery Assertions", io.ConfigErrorRecoveryEnableAssert)
    _debug_io_option!("Recovery Debug Log", io.ConfigErrorRecoveryEnableDebugLog)
    _debug_io_option!("Recovery Tooltips", io.ConfigErrorRecoveryEnableTooltip)
    ig.EndMenu()
    return nothing
end

"""Render the built-in Dear ImGui and ImPlot development tools inside the open Debug menu."""
function render_debug_tools_menu!(state::BrowserState)::Nothing
    if ig.MenuItem("Metrics / Debugger", C_NULL, state.show_imgui_metrics)
        state.show_imgui_metrics = !state.show_imgui_metrics
    end
    if ig.MenuItem("Debug Log", C_NULL, state.show_imgui_debug_log)
        state.show_imgui_debug_log = !state.show_imgui_debug_log
    end
    if ig.MenuItem("ID Stack Tool", C_NULL, state.show_imgui_id_stack)
        state.show_imgui_id_stack = !state.show_imgui_id_stack
    end

    debugger_attached = unsafe_load(ig.GetIO().ConfigDebugIsDebuggerPresent)
    if ig.MenuItem("Item Picker", C_NULL, false, debugger_attached)
        ig.DebugStartItemPicker()
    end
    if !debugger_attached &&
       ig.IsItemHovered(ig.ImGuiHoveredFlags_AllowWhenDisabled |
                        ig.ImGuiHoveredFlags_DelayNormal)
        if ig.BeginTooltip()
            ig.TextUnformatted("Attach a native debugger, then enable")
            ig.TextUnformatted("ConfigDebugIsDebuggerPresent in the Demo window.")
            ig.EndTooltip()
        end
    end
    _render_debug_options_menu!()

    ig.Separator()
    if ig.MenuItem("Style / Font Editor", C_NULL, state.show_imgui_style_editor)
        state.show_imgui_style_editor = !state.show_imgui_style_editor
    end
    if ig.MenuItem("User Guide", C_NULL, state.show_imgui_user_guide)
        state.show_imgui_user_guide = !state.show_imgui_user_guide
    end
    if ig.MenuItem("About Dear ImGui", C_NULL, state.show_imgui_about)
        state.show_imgui_about = !state.show_imgui_about
    end

    if ig.BeginMenu("ImPlot")
        _init_implot_context!(state)
        if ig.MenuItem("Metrics", C_NULL, state.show_implot_metrics)
            state.show_implot_metrics = !state.show_implot_metrics
        end
        if ig.MenuItem("Style Editor", C_NULL, state.show_implot_style_editor)
            state.show_implot_style_editor = !state.show_implot_style_editor
        end
        if ig.MenuItem("User Guide", C_NULL, state.show_implot_user_guide)
            state.show_implot_user_guide = !state.show_implot_user_guide
        end
        ig.Separator()
        if ig.MenuItem("Demo Window", C_NULL, state.show_implot_demo)
            state.show_implot_demo = !state.show_implot_demo
        end
        ig.EndMenu()
    end

    ig.Separator()
    if ig.MenuItem("Dear ImGui Demo", C_NULL, state.show_imgui_demo)
        state.show_imgui_demo = !state.show_imgui_demo
    end
    return nothing
end

"""Render Dear ImGui's style and font tools in one closable window."""
function _render_imgui_style_editor!(state::BrowserState)::Nothing
    state.show_imgui_style_editor || return nothing
    visible = @c ig.Begin("Dear ImGui Style / Font Editor", &state.show_imgui_style_editor)
    if visible
        ig.ShowStyleSelector("Style")
        ig.ShowFontSelector("Font")
        ig.Separator()
        ig.ShowStyleEditor()
    end
    ig.End()
    return nothing
end

"""Render Dear ImGui's interaction guide in a closable window."""
function _render_imgui_user_guide!(state::BrowserState)::Nothing
    state.show_imgui_user_guide || return nothing
    visible = @c ig.Begin("Dear ImGui User Guide", &state.show_imgui_user_guide)
    visible && ig.ShowUserGuide()
    ig.End()
    return nothing
end

"""Render ImPlot's configuration and style tools in one closable window."""
function _render_implot_style_editor!(state::BrowserState)::Nothing
    state.show_implot_style_editor || return nothing
    visible = @c ig.Begin("ImPlot Style Editor", &state.show_implot_style_editor)
    if visible
        ig.lib.ImPlot_ShowStyleSelector("Style")
        ig.lib.ImPlot_ShowColormapSelector("Colormap")
        ig.lib.ImPlot_ShowInputMapSelector("Input map")
        ig.Separator()
        ig.lib.ImPlot_ShowStyleEditor(C_NULL)
    end
    ig.End()
    return nothing
end

"""Render ImPlot's interaction guide in a closable window."""
function _render_implot_user_guide!(state::BrowserState)::Nothing
    state.show_implot_user_guide || return nothing
    visible = @c ig.Begin("ImPlot User Guide", &state.show_implot_user_guide)
    visible && ig.lib.ImPlot_ShowUserGuide()
    ig.End()
    return nothing
end

"""Render every open development tool and the optional Dear ImGui demo."""
function render_debug_tools!(state::BrowserState)::Nothing
    state.show_imgui_metrics && @c ig.ShowMetricsWindow(&state.show_imgui_metrics)
    state.show_imgui_debug_log && @c ig.ShowDebugLogWindow(&state.show_imgui_debug_log)
    state.show_imgui_id_stack && @c ig.ShowIDStackToolWindow(&state.show_imgui_id_stack)
    state.show_imgui_about && @c ig.ShowAboutWindow(&state.show_imgui_about)
    _render_imgui_style_editor!(state)
    _render_imgui_user_guide!(state)

    if state.show_implot_metrics || state.show_implot_style_editor ||
       state.show_implot_user_guide || state.show_implot_demo
        _init_implot_context!(state)
    end
    state.show_implot_metrics &&
        @c ig.lib.ImPlot_ShowMetricsWindow(&state.show_implot_metrics)
    _render_implot_style_editor!(state)
    _render_implot_user_guide!(state)
    state.show_implot_demo &&
        @c ig.lib.ImPlot_ShowDemoWindow(&state.show_implot_demo)

    state.show_imgui_demo && @c ig.ShowDemoWindow(&state.show_imgui_demo)
    return nothing
end
