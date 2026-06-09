using GLMakie
using DataFrames

using DataLoader: read_tlm_4p

function load_source_data(
    ::TASEProject,
    source_file::SourceFile;
    measurement::Union{Nothing,MeasurementInfo}=nothing,
)::DataFrame
    _check_cancel()
    if measurement !== nothing && measurement.measurement_kind !== :four_terminal_iv
        error("TASE cannot load measurement data for $(measurement.measurement_kind)")
    end
    detect_kind(TASE_PROJECT, source_file.filename) === :four_terminal_iv ||
        error("TASE cannot load source data for $(source_file.filepath)")
    return read_tlm_4p(source_file.filename, dirname(source_file.filepath))
end

function setup_plot(
    ::TASEProject,
    plot_kind::Type{TASEFourTerminalIVPlot},
    measurements::Vector{MeasurementInfo},
)::Figure
    isempty(measurements) && error("TASE plot requires at least one measurement")
    title = length(measurements) == 1 ? only(measurements).clean_title : "TASE IV overlay"
    fig = Figure(size=(700, 500))
    Axis(fig[1, 1], xlabel="Current (µA)", ylabel="Voltage Drop (mV)", title=title)
    return fig
end

function plot_data!(
    project::TASEProject,
    plot_kind::Type{TASEFourTerminalIVPlot},
    measurements::Vector{MeasurementInfo},
    figure::Figure,
)::Nothing
    axes = contents(figure[1, 1])
    isempty(axes) && error("TASE plot figure has no axis")
    ax = only(axes)
    data = read_measurement_data(project, measurements)
    for (measurement, df) in zip(measurements, data)
        nrow(df) == 0 && continue
        label = measurement.clean_title
        lines!(ax, df.current_source .* 1e6, df.voltage_drop .* 1e3; linewidth=2, label)
        scatter!(ax, df.current_source .* 1e6, df.voltage_drop .* 1e3; markersize=4)
    end
    length(measurements) > 1 && axislegend(ax)
    return nothing
end
