using Test

@testset "threaded workspace smoke test" begin
    mktempdir() do dir
        script_path = joinpath(dir, "threaded_workspace_smoke.jl")
        write(script_path, raw"""
using MeasurementBrowser

@assert Base.Threads.nthreads() >= 4

root = mktempdir()
try
    write(
        joinpath(root, "metadata.txt"),
        "collection_path,wafer,owner\nroot,W,collection\nroot/g0,W0,group\nroot/g1,W1,group\nroot/g2,W2,group\nroot/g3,W3,group\n",
    )
    for index in 1:48
        write(joinpath(root, "item_$index.txt"), "value=$index\n")
    end

    project = MeasurementBrowser.define_project("Threaded workspace smoke")
    MeasurementBrowser.register_item!(
        project,
        :text;
        detect=file -> endswith(file.filename, ".txt"),
        read=function (file)
            sleep(0.002)
            return read(file.filepath, String)
        end,
        entries=function (file, data)
            index = parse(Int, replace(splitext(file.filename)[1], "item_" => ""))
            group = "g$(index % 4)"
            return [MeasurementBrowser.DataItem(
                kind=:text,
                collection=["root", group],
                label=file.filename,
                parameters=Dict{Symbol,Any}(:index => index, :owner => "item"),
                data=data,
                id=file.filepath * "#text",
            )]
        end,
        stats=item -> Dict{Symbol,Any}(
            :bytes => length(item.data),
            :wafer => item.parameters[:wafer],
            :owner => item.parameters[:owner],
        ),
    )

    workspace = MeasurementBrowser.open_workspace(
        project,
        MeasurementBrowser.DirectorySource(root),
    )
    try
        deadline = time() + 20
        while time() < deadline
            MeasurementBrowser.Workspace.poll_workspace!(workspace)
            if workspace.scan.state in (:done, :unchanged) &&
               workspace.analysis.state == :done
                break
            end
            sleep(0.01)
        end

        @assert workspace.scan.state in (:done, :unchanged)
        @assert workspace.analysis.state == :done
        @assert length(workspace.index.items) == 48
        for record in values(workspace.index.items)
            @assert startswith(String(record.parameters[:wafer]), "W")
            @assert record.parameters[:owner] == "item"
            @assert record.stats[:bytes] > 0
            @assert record.stats[:owner] == "item"
        end
    finally
        MeasurementBrowser.close_workspace!(workspace)
    end
finally
    rm(root; recursive=true, force=true)
end
""")

        cmd = `$(Base.julia_cmd()) --project=$(Base.active_project()) --threads=4 --startup-file=no $script_path`
        @test success(cmd)
    end
end
