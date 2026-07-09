using DataBrowser
using Test

const OO_CACHE = DataBrowserCore.Cache

@testset "workspace open options capture and replay" begin
    mktempdir() do dir
        write(joinpath(dir, "sample.csv"), "a,b\n1,2\n")
        project = DataBrowser.define_project("OpenOptions_$(basename(dir))")
        workspace = DataBrowser.open_workspace(
            project, dir;
            recursive=false,
            metadata_file="custom_meta.txt",
            profile_internal=false,
            profile_cpu=false,
            profile_output=nothing,
            crash_trace=nothing,
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
            @test options.profile_internal == false
            @test options.profile_cpu == false
            @test options.profile_output === nothing
            @test options.crash_trace === nothing
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
