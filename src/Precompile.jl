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
        measurements = [
            only(measurements_for_file(TASE_PROJECT, joinpath(fixture_dir, filename)))
            for filename in filenames
        ]
        workspace = Workspace.Workspace(TASE_PROJECT, fixture_dir)
        read_measurement_data(workspace, measurements)
        figure = setup_plot(workspace, TASEFourTerminalIVPlot, measurements)
        plot_data!(workspace, TASEFourTerminalIVPlot, measurements, figure)
        infer_measurement_group("precompile_group", measurements[1:1], measurements)
        scan_source(scan_dir; project=TASE_PROJECT)
        scan_source(scan_dir; project=TASE_PROJECT, count_first=true)
    end
end
