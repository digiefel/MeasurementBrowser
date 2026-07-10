using DataBrowserGUI
using DataBrowserPlots
using Test

const Browser = DataBrowserGUI.Browser

@testset "plots extension boundary" begin
    registry = getfield(Browser, :_GUI_EXTENSION_TYPES)
    original = copy(registry)
    try
        empty!(registry)
        @test isempty(Browser._instantiate_extensions())

        Browser.register_gui_extension!(DataBrowserPlots.PlotsExtension)
        extensions = Browser._instantiate_extensions()
        @test length(extensions) == 1
        @test extensions[1] isa DataBrowserPlots.PlotsExtension
    finally
        empty!(registry)
        append!(registry, original)
    end
end
