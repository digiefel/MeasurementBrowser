using DataBrowser
using Test

const OO_CACHE = DataBrowserCache

@testset "source copy and modify_workspace!" begin
    mktempdir() do dir
        project = DataBrowser.define_project("OpenOptions_$(basename(dir))")
        workspace = DataBrowser.open_workspace(
            project, dir;
            recursive=false,
            metadata_file="custom_meta.txt",
            cache=false,
            background_processing=false,
        )
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

            previous_source = workspace.source
            modified = DataBrowser.modify_workspace!(workspace)
            @test modified === workspace
            @test workspace.source !== previous_source
            @test workspace.source.recursive == false
            @test workspace.source.metadata_file == "custom_meta.txt"
            @test workspace.cache.db isa OO_CACHE.MemoryCacheDB
            @test workspace.disk_cache == false
            @test workspace.background_processing == false
        finally
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

            DataBrowser.modify_workspace!(
                workspace;
                source=DataBrowser.DirectorySource(
                    dir; recursive=true, metadata_file=nothing),
            )
            @test workspace.source.recursive == true
            @test workspace.source.metadata_file === nothing
            @test workspace.disk_cache == false
        finally
            DataBrowser.close_workspace!(workspace)
        end
    end
end
