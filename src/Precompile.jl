using GLMakie: Axis, Figure, contents, lines!
using DataFrames: DataFrame, nrow

@setup_workload begin
    fixture_dir = normpath(joinpath(@__DIR__, "..", "test", "fixtures", "TASE"))
    filenames = [
        "TASESNS1c1f_A_2TSNJunction_11_20260224_111623_298K_FourTerminalIV.csv",
        "TASESNS1c1f_A_2TSNJunction_31_20260224_111700_298K_FourTerminalIV.csv",
    ]
    scan_dir = mktempdir()
    foreach(filename -> cp(
        joinpath(fixture_dir, filename),
        joinpath(scan_dir, filename);
        force=true,
    ), filenames)

    @compile_workload begin
        # Warm the registry pipeline (define -> register -> scan -> plot) the way real projects use it.
        project = define_project("Precompile")
        register_item!(
            project,
            :iv;
            detect=file -> endswith(file.filename, ".csv"),
            read=function (file)
                rows = Tuple{Float64,Float64}[]
                for line in Iterators.drop(readlines(file.filepath), 1)
                    isempty(strip(line)) && continue
                    parts = split(line, ',')
                    length(parts) >= 2 || continue
                    a = tryparse(Float64, parts[1])
                    b = tryparse(Float64, parts[2])
                    (a === nothing || b === nothing) && continue
                    push!(rows, (a, b))
                end
                DataFrame(i=first.(rows), v=last.(rows))
            end,
            entries=(file, data) -> [DataItem(
                kind=:iv,
                collection=[splitext(file.filename)[1]],
                label=file.filename,
                data=data,
            )],
            stats=item -> Dict{Symbol,Any}(:rows => nrow(item.data)),
        )
        register_plot!(
            project,
            :iv;
            label="I-V",
            setup=(_workspace, _items) -> (figure = Figure(); Axis(figure[1, 1]); figure),
            draw=function (_workspace, items, figure)
                axis = only(contents(figure[1, 1]))
                for item in items
                    df = item.data
                    nrow(df) == 0 && continue
                    lines!(axis, df.i, df.v)
                end
                nothing
            end,
        )

        records = [
            only(items_for_file(project, joinpath(fixture_dir, filename)))
            for filename in filenames
        ]
        precompile_depot = mktempdir()
        pushfirst!(DEPOT_PATH, precompile_depot)
        workspace = nothing
        try
            workspace = Workspace.Workspace(project, DirectorySource(fixture_dir))
            state = Browser.BrowserState(
                workspace=workspace,
                project_locked=true,
                project_preference=project_name(project),
            )
            Browser.current_status(state)
            Browser._project_visible_selection(state)
            Browser._make_timings_figure(state.performance.live_plots)
            Browser._make_build_figure(state.performance.live_plots)
            Browser._sample_build_progress!(state.performance.live_plots, workspace)
            Workspace.replace_item_index!(workspace, ItemIndex.Hierarchy(records, workspace.source))
            read_item_data(workspace, records)
            plot_kind = RegisteredPlot{:iv,Symbol("I-V")}
            items = Workspace.materialize_items(workspace, records)
            figure = setup_plot(workspace, plot_kind, items)
            plot_data!(workspace, plot_kind, items, figure)
        finally
            workspace === nothing || close_workspace!(workspace)
            popfirst!(DEPOT_PATH)
            rm(precompile_depot; force=true, recursive=true)
        end
        source = DirectorySource(scan_dir)
        scan_source(project, source)
        scan_source(project, source; count_first=true)
    end
end
