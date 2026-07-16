import CImGui as ig

@testset "debug tools state" begin
    state = DataBrowserGUI.Browser.BrowserState()

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
