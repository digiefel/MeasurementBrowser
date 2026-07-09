using CSV
using DataFrames: DataFrame, nrow
using DataBrowser
using Test

const MB = DataBrowser

function write_stress_csvs(root::AbstractString; files::Int=96, rows::Int=2000)::Nothing
    for file_index in 1:files
        path = joinpath(root, "stress_$(file_index).csv")
        open(path, "w") do io
            println(io, "x,y,z")
            for row in 1:rows
                println(io, "$(row),$(file_index + row / 10),$(file_index * row)")
            end
        end
    end
    return nothing
end

function build_csv_stress_project(items_per_file::Int)::MB.Project
    project = MB.define_project("GLFW scan startup stress")
    MB.register_item!(
        project,
        :table;
        detect=file -> endswith(file.filename, ".csv"),
        read=file -> DataFrame(CSV.File(file.filepath)),
        entries=(file, data) -> [
            MB.DataItem(
                kind=:table,
                collection=["stress"],
                label="$(file.filename) item $item_index",
                metadata=Dict{Symbol,Any}(:item_index => item_index),
                data=data,
                id=file.filepath * "#table_$item_index",
            )
            for item_index in 1:items_per_file
        ],
        analyze=item -> Dict{Symbol,Any}(:rows => nrow(item.data)),
    )
    return project
end

function run_browser_frames(
    workspace::DataBrowserCore.Workspace.Workspace;
    frames::Int=12,
    window_start::Symbol=:background,
)::Nothing
    state = DataBrowserGUI.Browser.BrowserState(
        project_locked=true,
        project_preference=workspace.project.name,
    )
    ctx = DataBrowserGUI.Browser._init_browser_context!()
    DataBrowserGUI.Browser._attach_workspace!(state, workspace)
    task = DataBrowserGUI.Browser._run_browser(
        state,
        ctx;
        engine=nothing,
        spawn=1,
        wait=false,
        exit_after_frames=frames,
        window_start=window_start,
    )
    wait(task)
    return nothing
end

if get(ENV, "MB_GLFW_SCAN_STRESS", "") == "1"
    @testset "GLFW startup during threaded scan stress" begin
        @test Threads.nthreads() >= 2
        mktempdir() do root
            files = parse(Int, get(ENV, "MB_GLFW_SCAN_STRESS_FILES", "96"))
            rows = parse(Int, get(ENV, "MB_GLFW_SCAN_STRESS_ROWS", "2000"))
            frames = parse(Int, get(ENV, "MB_GLFW_SCAN_STRESS_FRAMES", "12"))
            items_per_file = parse(Int, get(ENV, "MB_GLFW_SCAN_STRESS_ITEMS_PER_FILE", "1"))
            window_start = Symbol(get(ENV, "MB_GLFW_SCAN_STRESS_WINDOW_START", "background"))
            files > 0 || error("MB_GLFW_SCAN_STRESS_FILES must be positive, got $files")
            rows > 0 || error("MB_GLFW_SCAN_STRESS_ROWS must be positive, got $rows")
            frames > 0 || error("MB_GLFW_SCAN_STRESS_FRAMES must be positive, got $frames")
            items_per_file > 0 ||
                error("MB_GLFW_SCAN_STRESS_ITEMS_PER_FILE must be positive, got $items_per_file")
            window_start in (:normal, :background) ||
                error("MB_GLFW_SCAN_STRESS_WINDOW_START must be normal or background; got $window_start")
            @info "Running GLFW scan startup stress" files rows frames items_per_file window_start threads=Threads.nthreads()
            write_stress_csvs(root; files, rows)
            workspace = MB.open_workspace(
                build_csv_stress_project(items_per_file),
                DirectorySource(root),
            )
            try
                run_browser_frames(workspace; frames=frames, window_start=window_start)
                @test workspace.scan.state in (:done, :unchanged, :canceling, :canceled)
            finally
                MB.close_workspace!(workspace)
            end
        end
    end
else
    @info "Skipping opt-in GLFW stress test. Set MB_GLFW_SCAN_STRESS=1 to run it."
end
