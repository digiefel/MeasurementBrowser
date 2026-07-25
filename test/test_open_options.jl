using DataBrowser
using Test

const OO_CACHE = DataBrowserCache

@testset "workspace open options capture and replay" begin
    mktempdir() do dir
        project = DataBrowser.define_project("OpenOptions_$(basename(dir))")
        workspace = DataBrowser.open_workspace(
            project, dir;
            recursive=false,
            metadata_file="custom_meta.txt",
            cache=false,
            background_processing=false,
        )
        reopened = nothing
        try
            options = workspace.open_options
            @test options.recursive == false
            @test options.metadata_file == "custom_meta.txt"
            @test options.cache == false
            @test options.background_processing == false
            # `rebuild` is a one-shot action, never replayed on reopen.
            @test !haskey(options, :rebuild)

            reopened = DataBrowser.open_workspace(project, dir; options...)
            @test reopened.source.recursive == false
            @test reopened.source.metadata_file == "custom_meta.txt"
            @test reopened.cache.db isa OO_CACHE.MemoryCacheDB
            @test reopened.open_options == options
        finally
            reopened === nothing || DataBrowser.close_workspace!(reopened)
            DataBrowser.close_workspace!(workspace)
        end
    end

    mktempdir() do dir
        project = DataBrowser.define_project("OpenOptionsSource_$(basename(dir))")
        source = DataBrowser.DirectorySource(dir; recursive=false, metadata_file=nothing)
        workspace = DataBrowser.open_workspace(project, source; cache=false)
        try
            # The low-level source entry records the same complete set as the path entry.
            @test workspace.open_options.recursive == false
            @test workspace.open_options.metadata_file === nothing
            @test workspace.open_options.cache == false
        finally
            DataBrowser.close_workspace!(workspace)
        end
    end
end
