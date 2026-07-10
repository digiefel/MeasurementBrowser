using DataBrowserGUI
using Test

const Browser = DataBrowserGUI.Browser

mutable struct TestGuiExtension <: Browser.GuiExtension
    calls::Dict{Symbol,Int}
    ready::Bool
    persisted::Int
end

TestGuiExtension() = TestGuiExtension(Dict{Symbol,Int}(), false, 0)

function Browser.init!(ext::TestGuiExtension, ::Browser.BrowserState)
    ext.calls[:init] = get(ext.calls, :init, 0) + 1
    return nothing
end

function Browser.menu!(ext::TestGuiExtension, ::Browser.BrowserState)
    ext.calls[:menu] = get(ext.calls, :menu, 0) + 1
    return nothing
end

function Browser.draw!(ext::TestGuiExtension, ::Browser.BrowserState)
    ext.calls[:draw] = get(ext.calls, :draw, 0) + 1
    return nothing
end

function Browser.reset!(ext::TestGuiExtension, ::Browser.BrowserState)
    ext.calls[:reset] = get(ext.calls, :reset, 0) + 1
    ext.ready = false
    return nothing
end

function Browser.shutdown!(ext::TestGuiExtension, ::Browser.BrowserState)
    ext.calls[:shutdown] = get(ext.calls, :shutdown, 0) + 1
    return nothing
end

Browser.is_ready(ext::TestGuiExtension, ::Browser.BrowserState) = ext.ready

function Browser.save_view(ext::TestGuiExtension, ::Browser.BrowserState)
    return Dict{String,Any}("persisted" => ext.persisted)
end

function Browser.load_view!(ext::TestGuiExtension, ::Browser.BrowserState, view::Dict{String,Any})
    ext.persisted = Int(get(view, "persisted", 0))
    return nothing
end

@testset "gui extensions" begin
    registry = getfield(Browser, :_GUI_EXTENSION_TYPES)
    original = copy(registry)
    try
        empty!(registry)
        Browser.register_gui_extension!(TestGuiExtension)

        state = Browser.BrowserState()
        state.extensions = Browser._instantiate_extensions()
        @test length(state.extensions) == 1
        ext = state.extensions[1]
        @test ext isa TestGuiExtension

        Browser.init!(ext, state)
        @test ext.calls[:init] == 1

        Browser.menu!(ext, state)
        Browser.draw!(ext, state)
        @test ext.calls[:menu] == 1
        @test ext.calls[:draw] == 1

        ext.ready = true
        @test Browser._extensions_ready(state)

        ext.persisted = 7
        view = Browser.PersistedProjectView(
            project="ext-test",
            extensions=Dict(
                "UnknownExtension" => Dict{String,Any}("kept" => true),
            ),
        )
        state.saved_project_view = view
        merged = Browser._persisted_extensions(state)
        @test merged["UnknownExtension"] == Dict{String,Any}("kept" => true)
        @test merged["TestGuiExtension"] == Dict{String,Any}("persisted" => 7)

        view.extensions["TestGuiExtension"] = Dict{String,Any}("persisted" => 42)
        Browser.load_view!(ext, state, view.extensions["TestGuiExtension"])
        @test ext.persisted == 42

        Browser.reset!(ext, state)
        @test ext.calls[:reset] == 1
        @test !ext.ready

        Browser.shutdown!(ext, state)
        @test ext.calls[:shutdown] == 1
    finally
        empty!(registry)
        append!(registry, original)
    end
end
