using CSV
using DataFrames: DataFrame, nrow
using GLMakie: Axis, Figure, contents, lines!

"""
Minimal fixture project used only to exercise engine pipelines during precompilation.

User project callbacks are not precompilable for other projects; this registry exists solely to
drive scan, interpret, process, analyze, collection analysis, cache writes, and plotting through the
package machinery using bundled test fixtures.
"""
function _precompile_project()
    project = define_project("Precompile")
    register_item!(project, :iv;
        detect=file -> endswith(file.filename, ".csv"),
        read=file -> DataFrame(CSV.File(file.filepath; ntasks=1)),
        label=(data, metadata) -> metadata[:filename],
        collection=(data, metadata) -> [splitext(metadata[:filename])[1]],
        process=function (input, metadata)
            data = copy(input)
            data.engine_warm = data.Current_A ./ max.(abs.(data.VoltageHigh_V), 1e-12)
            return data
        end,
        analyze=(data, metadata) -> Dict{Symbol,Any}(:rows => nrow(data)),
    )
    register_collection_analysis!(project, :iv;
        analyze=(data, metadata) -> Dict{Symbol,Any}(:items => length(data)),
    )
    register_plot!(project, :iv; label="I-V",
        setup=(_workspace, _items) -> (figure = Figure(); Axis(figure[1, 1]); figure),
        draw=function (_workspace, items, figure)
            axis = only(contents(figure[1, 1]))
            for item in items
                df = item_data(item)
                nrow(df) == 0 && continue
                hasproperty(df, :VoltageHigh_V) && hasproperty(df, :Current_A) || continue
                lines!(axis, df.VoltageHigh_V, df.Current_A)
            end
            nothing
        end)
    return project
end

@setup_workload begin
    fixture_dir = normpath(joinpath(@__DIR__, "..", "test", "fixtures", "TASE"))
    filenames = [
        "TASESNS1c1f_A_2TSNJunction_11_20260224_111623_298K_FourTerminalIV.csv",
        "TASESNS1c1f_A_2TSNJunction_31_20260224_111700_298K_FourTerminalIV.csv",
    ]
    fixture_paths = [joinpath(fixture_dir, filename) for filename in filenames]

    @compile_workload begin
        project = _precompile_project()
        records = [only(items_for_file(project, path)) for path in fixture_paths]

        precompile_depot = mktempdir()
        pushfirst!(DEPOT_PATH, precompile_depot)
        workspace = nothing
        try
            workspace = open_workspace(
                project,
                fixture_dir;
            )
            Workspace.wait_workspace_idle!(workspace; timeout=45)
            Workspace.refresh_status!(workspace)
            records = collect(values(workspace.index.items))
            if !isempty(records)
                select_items!(workspace, records)
                items = Workspace.materialize_items(workspace, records)
                plot_kind = only(registered_plot_kinds(project, :iv))
                figure = setup_plot(workspace, plot_kind, items)
                plot_data!(workspace, plot_kind, items, figure)
                Cache.read_item_data(workspace.cache.db, records; stage=:interpreted)
                Cache.read_item_data(workspace.cache.db, records; stage=:processed)
            end
            Workspace.workspace_memory_snapshot(workspace)
            Cache.cache_built(workspace.cache.db) && Cache.load_cache_index(workspace.cache.db)

            state = Browser.BrowserState(
                workspace=workspace,
                project_locked=true,
                project_preference=project_name(project),
            )
            Browser.current_status(state)
            Browser._project_visible_selection(state)
            Browser._update_live_timings!(state.performance.live_plots, state.performance.timings)
            Browser._sample_build!(state.performance.live_plots, workspace)
        finally
            workspace === nothing || close_workspace!(workspace)
            popfirst!(DEPOT_PATH)
            rm(precompile_depot; force=true, recursive=true)
        end
    end
end
