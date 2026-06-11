using GLMakie
using DataFrames: nrow

"""Create the figure used for one or more TASE four-terminal measurements."""
function setup_plot(
    ::Workspace.Workspace{TASEProject},
    plot_kind::Type{TASEFourTerminalIVPlot},
    measurements::Vector{MeasurementInfo},
)::Figure
    isempty(measurements) && error("TASE plot requires at least one measurement")
    title = length(measurements) == 1 ? only(measurements).clean_title : "TASE IV overlay"
    fig = Figure(size=(700, 500))
    Axis(fig[1, 1], xlabel="Current (µA)", ylabel="Voltage Drop (mV)", title=title)
    return fig
end

"""Draw the selected TASE measurements into an existing figure."""
function plot_data!(
    workspace::Workspace.Workspace{TASEProject},
    plot_kind::Type{TASEFourTerminalIVPlot},
    measurements::Vector{MeasurementInfo},
    figure::Figure,
)::Nothing
    axes = contents(figure[1, 1])
    isempty(axes) && error("TASE plot figure has no axis")
    ax = only(axes)
    data = read_measurement_data(workspace, measurements)
    for (measurement, df) in zip(measurements, data)
        nrow(df) == 0 && continue
        label = measurement.clean_title
        lines!(ax, df.i .* 1e6, df.v .* 1e3; linewidth=2, label)
        scatter!(ax, df.i .* 1e6, df.v .* 1e3; markersize=4)
    end
    length(measurements) > 1 && axislegend(ax)
    return nothing
end
