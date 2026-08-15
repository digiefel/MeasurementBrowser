using DataBrowser
using Test

const OO_CACHE = DataBrowserCache

@testset "source copy and workspace reopen" begin
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
            @test workspace.source.recursive == false
            @test workspace.source.metadata_file == "custom_meta.txt"
            @test workspace.disk_cache == false
            @test workspace.background_processing == false

            cloned = copy(workspace.source)
            @test cloned isa DataBrowser.DirectorySource
            @test cloned !== workspace.source
            @test cloned.root_path == workspace.source.root_path
            @test cloned.recursive == false
            @test cloned.metadata_file == "custom_meta.txt"
            @test cloned.watcher_task === nothing

            reopened = DataBrowser.open_workspace(
                project,
                cloned;
                cache=workspace.disk_cache,
                background_processing=workspace.background_processing,
            )
            @test reopened.source.recursive == false
            @test reopened.source.metadata_file == "custom_meta.txt"
            @test reopened.cache.db isa OO_CACHE.MemoryCacheDB
            @test reopened.disk_cache == false
            @test reopened.background_processing == false
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
            @test workspace.source.recursive == false
            @test workspace.source.metadata_file === nothing
            @test workspace.disk_cache == false
        finally
            DataBrowser.close_workspace!(workspace)
        end
    end
end
