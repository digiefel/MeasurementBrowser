module MakieImguiIntegration

import CImGui as ig
import ModernGL as gl
import GLFW
import GLMakie
import GLMakie.Makie as Makie

"""GLFW window adapter required by a GLMakie screen embedded in CImGui."""
struct ImMakieWindow
    glfw_window::GLFW.Window
end

"""One embedded Makie figure, its reusable screen, and pointer-gesture state."""
mutable struct ImMakieFigure
    figure::GLMakie.Figure
    screen::GLMakie.Screen{ImMakieWindow}
    was_dragging::Bool
    was_ctrl_left_click::Bool
end

const MAKIE_CONTEXT = Dict{String,ImMakieFigure}()

"""Destroy every Makie screen owned by the browser render context."""
function destroy_context!()::Nothing
    for imfigure in values(MAKIE_CONTEXT)
        empty!(imfigure.figure)
        GLMakie.destroy!(imfigure.screen)
    end

    empty!(MAKIE_CONTEXT)
    return nothing
end

"""Destroy the embedded Makie screen and figure owned by one plot id."""
function destroy_figure!(title_id::String)::Nothing
    imfigure = pop!(MAKIE_CONTEXT, title_id, nothing)
    imfigure === nothing && return nothing
    empty!(imfigure.figure)
    GLMakie.destroy!(imfigure.screen)
    return nothing
end

Base.isopen(window::ImMakieWindow) = isopen(window.glfw_window)
GLFW.SwapBuffers(::ImMakieWindow) = nothing
GLMakie.destroy!(::ImMakieWindow) = nothing
GLMakie.reopen!(screen::GLMakie.Screen{ImMakieWindow}) = screen
GLMakie.set_screen_visibility!(::GLMakie.Screen{ImMakieWindow}, ::Bool) = nothing
GLMakie.framebuffer_size(window::ImMakieWindow) = GLMakie.framebuffer_size(window.glfw_window)
GLMakie.scale_factor(window::ImMakieWindow) = GLMakie.scale_factor(window.glfw_window)
GLMakie.was_destroyed(window::ImMakieWindow) = GLMakie.was_destroyed(window.glfw_window)

GLMakie.ShaderAbstractions.native_switch_context!(x::ImMakieWindow) = GLFW.MakeContextCurrent(x.glfw_window)
GLMakie.ShaderAbstractions.native_context_alive(x::ImMakieWindow) = GLFW.is_initialized() && x.glfw_window != C_NULL

# CImGui forwards input directly below, so GLMakie's normal GLFW event connection must stay disabled.
GLMakie.connect_screen(::GLMakie.Scene, ::GLMakie.Screen{ImMakieWindow}) = nothing

for f in (Makie.window_area, Makie.window_open,
    Makie.mouse_buttons, Makie.mouse_position,
    Makie.scroll,
    Makie.keyboard_buttons, Makie.unicode_input,
    Makie.dropped_files,
    Makie.hasfocus, Makie.entered_window,
    Makie.frame_tick)
    GLMakie.disconnect!(::GLMakie.Screen{ImMakieWindow}, ::typeof(f)) = nothing
end

"""Convert framebuffer pixels into the logical size expected by CImGui."""
function _texture_display_size(
    texture_size::Tuple{Int,Int},
    px_per_unit::Real,
)::Tuple{Float64,Float64}
    scale = Float64(px_per_unit)
    scale > 0 || return (Float64(texture_size[1]), Float64(texture_size[2]))
    return (Float64(texture_size[1]) / scale, Float64(texture_size[2]) / scale)
end

"""Convert CImGui's top-left mouse coordinates to Makie's bottom-left coordinates."""
function _imgui_mouse_to_makie(
    pos::ig.ImVec2,
    origin::ig.ImVec2,
    image_size::Tuple{Float64,Float64},
)::Tuple{Float64,Float64}
    return (Float64(pos.x - origin.x), Float64(image_size[2] - (pos.y - origin.y)))
end

_ctrl_left_as_right_click(io) = unsafe_load(io.MouseCtrlLeftAsRightClick)
_makie_button_for_imgui_right(ctrl_left_as_right, was_ctrl_left_click=false) =
    (ctrl_left_as_right || was_ctrl_left_click) ? Makie.Mouse.left : Makie.Mouse.right
_should_open_axis_popup(was_dragging, was_ctrl_left_click, ctrl_left_as_right) =
    !(was_dragging || was_ctrl_left_click || ctrl_left_as_right)

"""Render the scale choices for one Makie axis direction."""
function _draw_axis_scale_buttons!(scale, index::Int)::Nothing
    if ig.RadioButton("linear##$(index)", scale[] === identity)
        scale[] = identity
    end
    ig.SameLine()
    if ig.RadioButton("pseudolog10##$(index)", scale[] === Makie.pseudolog10)
        scale[] = Makie.pseudolog10
    end
    ig.SameLine()
    if ig.RadioButton("log10##$(index)", scale[] === log10)
        scale[] = log10
    end
    return nothing
end

"""Render the context menu for one hovered Makie axis."""
function _draw_axis_popup!(axis::Makie.Axis)::Nothing
    if ig.BeginPopupContextItem()
        ig.Text("Axis settings")
        ig.Separator()

        ig.Text("X scale:")
        ig.SameLine()
        _draw_axis_scale_buttons!(axis.xscale, 1)

        ig.Text("Y scale:")
        ig.SameLine()
        _draw_axis_scale_buttons!(axis.yscale, 2)

        ig.EndPopup()
    end
    return nothing
end

"""
Render one Makie figure as an interactive CImGui image.

Screens are reused by `title_id`; replacing the figure preserves the OpenGL context.
"""
function MakieFigure(
    title_id::String,
    figure::GLMakie.Figure;
    auto_resize_x::Bool=true,
    auto_resize_y::Bool=false,
)::Bool
    ig.PushID(title_id)
    id = title_id

    if haskey(MAKIE_CONTEXT, id)
        imf = MAKIE_CONTEXT[id]
        if imf.figure !== figure
            # Reuse the existing screen, but only clear the previous figure
            # after GLMakie has detached it from the screen. Clearing it before
            # display(screen, f) prevents GLMakie from walking the old scene and
            # releasing its render objects.
            old_figure = imf.figure
            scene = Makie.get_scene(figure)
            scene.events.window_open[] = true
            display(imf.screen, figure)
            empty!(old_figure)
            imf.figure = figure
            @debug "replaced figure for " id
        end
    else
        window = ig.current_window()
        makie_window = ImMakieWindow(window)
        screen = GLMakie.Screen(; window=makie_window, start_renderloop=false)
        MAKIE_CONTEXT[id] = ImMakieFigure(figure, screen, false, false)
        scene = Makie.get_scene(figure)
        scene.events.window_open[] = true
        display(screen, figure)
        @debug "created context for " id
    end

    imfigure = MAKIE_CONTEXT[id]
    scene = Makie.get_scene(imfigure.figure)

    region_avail = ig.GetContentRegionAvail()
    region_size = (Int(region_avail.x), Int(region_avail.y))
    scene_size = size(scene)
    new_size = (auto_resize_x ? region_size[1] : scene_size[1],
        auto_resize_y ? region_size[2] : scene_size[2])

    if scene_size != new_size && all(new_size .> 0)
        scene.events.window_area[] = GLMakie.Rect2i(0, 0, Int(new_size[1]), Int(new_size[2]))
        resize!(figure, new_size[1], new_size[2])
    end

    GLMakie.poll_updates(imfigure.screen)
    do_render = GLMakie.requires_update(imfigure.screen)
    do_render && GLMakie.render_frame(imfigure.screen; resize_buffers=false)

    gl.glBindFramebuffer(gl.GL_FRAMEBUFFER, 0) # Unbind the framebuffer to avoid issues with ImGui rendering

    # The color texture is what we need to render as an image. We add it to the
    # drawlist and then create an InvisibleButton of the same size to create a
    # space in the layout that can respond to key presses and clicks etc (which
    # a regular Image() can't do).
    color_buffer = imfigure.screen.framebuffer.buffers[:color]
    drawlist = ig.GetWindowDrawList()
    cursor_pos = ig.GetCursorScreenPos()
    texture_size = size(color_buffer)
    image_size = _texture_display_size(texture_size, imfigure.screen.px_per_unit[])
    ig.AddImage(drawlist,
        ig.ImTextureRef(ig.ImTextureID(color_buffer.id)),
        cursor_pos,
        (cursor_pos.x + image_size[1], cursor_pos.y + image_size[2]),
        (0, 1), (1, 0))

    ig.InvisibleButton("figure_image", image_size)

    # Update the scene events
    if scene.events.hasfocus[] != ig.IsItemHovered()
        scene.events.hasfocus[] = ig.IsItemHovered()
    end
    if scene.events.entered_window[] != ig.IsItemHovered()
        scene.events.entered_window[] = ig.IsItemHovered()
    end

    io = ig.GetIO()
    if ig.IsItemHovered()
        pos = ig.GetMousePos()
        new_pos = _imgui_mouse_to_makie(pos, cursor_pos, image_size)
        if new_pos != scene.events.mouseposition[]
            scene.events.mouseposition[] = new_pos
            @debug "mouse position updated", scene.events.mouseposition
        end

        ctrl_left_as_right = _ctrl_left_as_right_click(io)

        # Makie gestures read modifier state when the mouse event arrives.
        for (igkey, glfwkey) in ((ig.lib.ImGuiKey_X, Int(GLFW.KEY_X)),
            (ig.lib.ImGuiKey_Y, Int(GLFW.KEY_Y)),
            (ig.lib.ImGuiKey_LeftCtrl, Int(GLFW.KEY_LEFT_CONTROL)))
            if ig.IsKeyPressed(igkey) ||
                    (ctrl_left_as_right && igkey == ig.lib.ImGuiKey_LeftCtrl)
                event = Makie.KeyEvent(Makie.Keyboard.Button(glfwkey), Makie.Keyboard.press)
                scene.events.keyboardbutton[] = event
                @debug "keyboard button pressed event", scene.events.keyboardbutton
            elseif ig.IsKeyReleased(igkey)
                event = Makie.KeyEvent(Makie.Keyboard.Button(glfwkey), Makie.Keyboard.release)
                scene.events.keyboardbutton[] = event
                @debug "keyboard button released event", scene.events.keyboardbutton
            end
        end

        # macOS may report Ctrl+left-click as right-click; translate it back for Makie.
        if ig.IsKeyPressed(ig.ImGuiKey_MouseLeft)
            scene.events.mousebutton[] = Makie.MouseButtonEvent(Makie.Mouse.left, Makie.Mouse.press)
            @debug "mouse button pressed event", scene.events.mousebutton
        elseif ig.IsKeyReleased(ig.ImGuiKey_MouseLeft)
            scene.events.mousebutton[] = Makie.MouseButtonEvent(Makie.Mouse.left, Makie.Mouse.release)
            @debug "mouse button released event", scene.events.mousebutton
        end

        if ig.IsKeyPressed(ig.ImGuiKey_MouseRight)
            makie_button = _makie_button_for_imgui_right(ctrl_left_as_right)
            imfigure.was_ctrl_left_click = ctrl_left_as_right
            scene.events.mousebutton[] = Makie.MouseButtonEvent(makie_button, Makie.Mouse.press)
            @debug "mouse button pressed event", scene.events.mousebutton
        elseif ig.IsKeyReleased(ig.ImGuiKey_MouseRight)
            makie_button = _makie_button_for_imgui_right(false, imfigure.was_ctrl_left_click)
            scene.events.mousebutton[] = Makie.MouseButtonEvent(makie_button, Makie.Mouse.release)
            imfigure.was_ctrl_left_click = false
            @debug "mouse button released event", scene.events.mousebutton
        end

        if !imfigure.was_ctrl_left_click && ig.IsMouseDragging(ig.lib.ImGuiMouseButton_Right)
            imfigure.was_dragging = true
        end

        wheel_y = unsafe_load(io.MouseWheel)
        wheel_x = unsafe_load(io.MouseWheelH)
        if (wheel_x, wheel_y) != scene.events.scroll[]
            scene.events.scroll[] = (wheel_x, wheel_y)
            @debug "scrolling" scene.events.scroll
        end
    elseif imfigure.was_ctrl_left_click && ig.IsKeyReleased(ig.ImGuiKey_MouseRight)
        scene.events.mousebutton[] = Makie.MouseButtonEvent(Makie.Mouse.left, Makie.Mouse.release)
        imfigure.was_ctrl_left_click = false
        @debug "mouse button released event", scene.events.mousebutton
    end

    ctrl_left_as_right = _ctrl_left_as_right_click(ig.GetIO())
    if !ctrl_left_as_right && ig.IsMouseClicked(ig.lib.ImGuiMouseButton_Right)
        imfigure.was_dragging = false
    end

    if _should_open_axis_popup(imfigure.was_dragging, imfigure.was_ctrl_left_click, ctrl_left_as_right)
        for content in figure.content
            if content isa Makie.Axis && GLMakie.is_mouseinside(content)
                _draw_axis_popup!(content)
            end
        end
    end

    ig.PopID()

    return do_render
end

"""Return the Makie theme used by figures embedded in CImGui."""
function _imgui_theme()::Makie.Attributes
    theme = Makie.theme_light()
    theme.Legend.framevisible = true
    theme.Legend.padding = (5, 5, 5, 5)
    theme.GLMakie = Makie.Attributes(fxaa=false)

    return theme
end

function __init__()::Nothing
    ig.atrenderexit(destroy_context!)
    Makie.set_theme!(_imgui_theme())
    return nothing
end

end
