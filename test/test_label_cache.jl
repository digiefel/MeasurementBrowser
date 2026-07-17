using DataBrowser
using Test
using DataBrowserCore.Workspace: wait_workspace_idle!

# Labels are resolved once from live values during interpretation and persisted on the records.
# A warm reopen must deliver them without running any project label code; a changed source item
# must re-resolve them on the next scan.

const LABEL_CALLS = Ref(0)

function _label_cache_project(name::AbstractString)
    project = DataBrowser.define_project(name)
    DataBrowser.register_item!(
        project,
        :reading;
        detect=file -> endswith(file.filename, ".csv"),
        read=file -> begin
            value = parse(Float64, split(readchomp(file.filepath), '\n')[2])
            (x=[value],)
        end,
        label=(data, _metadata) -> begin
            LABEL_CALLS[] += 1
            "reading-$(only(data.x))"
        end,
    )
    return project
end

@testset "record labels are cached and refreshed when the source item changes" begin
    mktempdir() do dir
        filepath = joinpath(dir, "a.csv")
        write(filepath, "x\n1\n")
        project = _label_cache_project("LabelCache_$(basename(dir))")

        LABEL_CALLS[] = 0
        workspace = DataBrowser.open_workspace(
            project, test_source(project, dir); background_processing=true)
        try
            wait_workspace_idle!(workspace)
            record = only(collect(values(workspace.index.items)))
            @test DataBrowser.label(record) == "reading-1.0"
            @test LABEL_CALLS[] == 1
            source_row = only(values(read(workspace.cache.db.source_items)))
            @test DataBrowser.label(source_row) == "a.csv"
        finally
            DataBrowser.close_workspace!(workspace)
        end

        # Warm reopen: the persisted label is delivered from the record without rerunning
        # read or label project code.
        LABEL_CALLS[] = 0
        reopened = DataBrowser.open_workspace(
            project, test_source(project, dir); background_processing=true)
        try
            wait_workspace_idle!(reopened)
            record = only(collect(values(reopened.index.items)))
            @test DataBrowser.label(record) == "reading-1.0"
            @test LABEL_CALLS[] == 0
        finally
            DataBrowser.close_workspace!(reopened)
        end

        # The file changed while the workspace was closed: the stale label is re-resolved
        # from the live value on the next scan.
        write(filepath, "x\n2\n")
        LABEL_CALLS[] = 0
        refreshed = DataBrowser.open_workspace(
            project, test_source(project, dir); background_processing=true)
        try
            wait_workspace_idle!(refreshed)
            record = only(collect(values(refreshed.index.items)))
            @test DataBrowser.label(record) == "reading-2.0"
            @test LABEL_CALLS[] == 1
            source_row = only(values(read(refreshed.cache.db.source_items)))
            @test DataBrowser.label(source_row) == "a.csv"
        finally
            DataBrowser.close_workspace!(refreshed)
        end
    end
end
