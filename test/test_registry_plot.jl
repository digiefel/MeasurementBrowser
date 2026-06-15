using MeasurementBrowser
using CSV
using DataFrames: DataFrame
using GLMakie: Figure, Axis
using Test

const MB = MeasurementBrowser

@testset "registry plot bridge" begin
    dir = mktempdir()
    write(joinpath(dir, "m.csv"), "x,y\n1,2\n3,4\n")

    project = MB.define_project("PlotBridge")
    MB.register_measurement!(
        project,
        :table;
        detect=file -> endswith(file.filename, ".csv"),
        read=file -> DataFrame(CSV.File(file.filepath)),
        measurements=(file, data) -> [MB.MeasurementInfo(
            filepath=file.filepath,
            measurement_kind=:table,
            device_info=MB.DeviceInfo(["dev"]),
            timestamp=file.timestamp,
            clean_title=file.filename,
        )],
    )
    drew = Ref(0)
    MB.register_plot!(
        project,
        :table;
        label="Table Plot",
        setup=(workspace, measurements, processed) -> Figure(),
        draw=(workspace, measurements, processed, figure) -> (drew[] = length(processed); Axis(figure[1, 1]); nothing),
    )
    MB.register_plot!(
        project,
        :table;
        label="Table Summary",
        setup=(workspace, measurements, processed) -> Figure(),
        draw=(workspace, measurements, processed, figure) -> nothing,
    )
    MB.register_plot!(
        project,
        :table;
        label="Table Summary",
        setup=(workspace, measurements, processed) -> Figure(),
        draw=(workspace, measurements, processed, figure) -> (drew[] = -1; nothing),
    )

    # Identity helpers used by the GUI and persistence.
    table_plot = MB.RegistryPlot{:table,Symbol("Table Plot")}
    table_summary = MB.RegistryPlot{:table,Symbol("Table Summary")}
    @test MB.registered_plot_kinds(project, :table) == [table_plot, table_summary]
    @test MB.registered_plot_kinds(project, :missing) == Type{<:MB.PlotKind}[]
    @test MB.plot_kind_label(project, table_plot) == "Table Plot"
    @test MB.plot_kind_name(table_plot) == "table::Table Plot"
    @test MB.plot_kind_from_name("table::Table Plot") === table_plot
    @test MB.plot_kind_from_name("table") === nothing

    # The bridge runs the registered setup/draw callbacks via the engine's plot dispatch.
    measurements = measurements_for_file(project, joinpath(dir, "m.csv"))
    workspace = MB.Workspace.Workspace(project, dir)
    figure = setup_plot(workspace, table_plot, measurements)
    @test figure isa Figure
    @test plot_data!(workspace, table_plot, measurements, figure) === nothing
    @test drew[] == length(measurements)
    figure = setup_plot(workspace, table_summary, measurements)
    @test plot_data!(workspace, table_summary, measurements, figure) === nothing
    @test drew[] == -1

    # An unregistered kind errors clearly rather than dispatching nowhere.
    @test_throws KeyError setup_plot(
        workspace,
        MB.RegistryPlot{:missing,Symbol("Missing")},
        measurements,
    )
end
