import CImGui as ig

@testset "debug tools state" begin
    state = DataBrowserGUI.Browser.BrowserState()

    @test !state.show_imgui_metrics
    @test !state.show_imgui_debug_log
    @test !state.show_imgui_id_stack
    @test !state.show_imgui_style_editor
    @test !state.show_imgui_user_guide
    @test !state.show_imgui_about
    @test !state.show_imgui_demo
    @test !state.show_implot_metrics
    @test !state.show_implot_style_editor
    @test !state.show_implot_user_guide
    @test !state.show_implot_demo
    @test state.implot_context == C_NULL

    imgui_context = ig.CreateContext()
    try
        DataBrowserGUI.Browser._init_implot_context!(state)
        @test state.implot_context != C_NULL
        @test ig.lib.ImPlot_GetCurrentContext() == state.implot_context

        context = state.implot_context
        DataBrowserGUI.Browser._init_implot_context!(state)
        @test state.implot_context == context

        DataBrowserGUI.Browser._shutdown_implot_context!(state)
        @test state.implot_context == C_NULL
    finally
        DataBrowserGUI.Browser._shutdown_implot_context!(state)
        ig.DestroyContext(imgui_context)
    end
end
