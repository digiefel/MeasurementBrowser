using Test, DataBrowserPlots, DataBrowserGUI

@testset "plot views preserve fixed and live selections" begin
    browser = DataBrowserGUI.Browser
    state = browser.BrowserState()
    original = Dict{String,Any}(
        "main_plot" => Dict("id" => "main", "title" => "Live", "plot_kind" => "",
            "live" => true, "items" => ["one"]),
        "plot_windows" => [Dict("id" => "plot_1", "title" => "Comparison", "plot_kind" => "",
            "live" => false, "items" => ["one", "two"])])
    first = PlotsExtension()
    browser.load_view!(first, state, original)
    saved = browser.save_view(first, state)
    restored = PlotsExtension()
    browser.load_view!(restored, state, saved)
    @test browser.save_view(restored, state) == saved
    @test saved["main_plot"] == original["main_plot"]
    @test saved["plot_windows"] == original["plot_windows"]
end
