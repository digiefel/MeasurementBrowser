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

    # Identity helpers used by the GUI and persistence.
    @test MB.registered_plot_kind(project, :table) === MB.RegistryPlot{:table}
    @test MB.registered_plot_kind(project, :missing) === nothing
    @test MB.plot_kind_label(project, MB.RegistryPlot{:table}) == "Table Plot"
    @test MB.plot_kind_name(MB.RegistryPlot{:table}) == "table"
    @test MB.plot_kind_from_name("table") === MB.RegistryPlot{:table}

    # The bridge runs the registered setup/draw callbacks via the engine's plot dispatch.
    measurements = measurements_for_file(project, joinpath(dir, "m.csv"))
    workspace = MB.Workspace.Workspace(project, dir)
    figure = setup_plot(workspace, MB.RegistryPlot{:table}, measurements)
    @test figure isa Figure
    @test plot_data!(workspace, MB.RegistryPlot{:table}, measurements, figure) === nothing
    @test drew[] == length(measurements)

    # An unregistered kind errors clearly rather than dispatching nowhere.
    @test_throws ErrorException setup_plot(workspace, MB.RegistryPlot{:missing}, measurements)
end
