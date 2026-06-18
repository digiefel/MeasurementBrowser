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
            entries=(file, _data) -> [DataItem(
                kind=:iv,
                collection=[splitext(file.filename)[1]],
                label=file.filename,
            )],
            stats=(_item, data) -> Dict{Symbol,Any}(:rows => nrow(data)),
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
        workspace = Workspace.Workspace(project, fixture_dir)
        read_item_data(workspace, records)
        plot_kind = RegisteredPlot{:iv,Symbol("I-V")}
        figure = setup_plot(workspace, plot_kind, records)
        plot_data!(workspace, plot_kind, records, figure)
        source = RegisteredProjectSource(project, scan_dir)
        scan_source(project, source)
        scan_source(project, source; count_first=true)
    end
end
