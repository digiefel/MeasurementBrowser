"""
RuO2Project.jl - Project dispatch methods for RuO2test ferroelectric measurements
"""

using Dates
include("RuO2/Display.jl")
include("RuO2/CVSweepIO.jl")
include("RuO2/Interpretation.jl")
include("RuO2/PUNDIO.jl")
include("RuO2/PUNDAnalysis.jl")
include("RuO2/Stats.jl")
include("RuO2/PlotHelpers.jl")

struct RuO2PUNDPlot <: PlotKind end
struct RuO2IVSweepPlot <: PlotKind end
struct RuO2TLM4PointPlot <: PlotKind end
struct RuO2CVSweepPlot <: PlotKind end
struct RuO2TLMAnalysisPlot <: PlotKind end
struct RuO2TLMTemperaturePlot <: PlotKind end
struct RuO2PUNDFatiguePlot <: PlotKind end

include("RuO2/PlotSingle.jl")
include("RuO2/PlotCombined.jl")

# ---------------------------------------------------------------------------
# Interface implementations
# ---------------------------------------------------------------------------

project_name(::RuO2Project) = "RuO2"
project_description(::RuO2Project) = "Ferroelectric RuO2test measurements"

function _plot_job_data(
    project::RuO2Project,
    plot_kind::Type{<:PlotKind},
    measurements::Vector{MeasurementInfo};
    debug::Bool=false,
)
    if debug
        plot_kind === RuO2PUNDPlot ||
            error("Debug plots are not implemented for $(project_name(project)) $(plot_kind)")
        length(measurements) == 1 ||
            error("RuO2 PUND debug plot requires exactly one measurement")
        measurement = only(measurements)
        df = only(read_measurement_data(project, [measurement]))
        return _ruo2_plot_data(measurement, df; debug)
    end

    return nothing
end

function _plot_job_figure(
    project::RuO2Project,
    plot_kind::Type{<:PlotKind},
    measurements::Vector{MeasurementInfo},
    data;
    debug::Bool=false,
    device_params::Vector{Dict{Symbol,Any}}=Dict{Symbol,Any}[],
)
    if debug
        plot_kind === RuO2PUNDPlot ||
            error("Debug plots are not implemented for $(project_name(project)) $(plot_kind)")
        return debug_plot(project, measurements, data; device_params, plot_kind)
    end

    fig = setup_plot(project, plot_kind, measurements)
    plot_data!(project, plot_kind, measurements, fig)
    return fig
end

# ---------------------------------------------------------------------------
# Registration
# ---------------------------------------------------------------------------

const RUO2_PROJECT = RuO2Project()
push!(KNOWN_PROJECTS, RUO2_PROJECT)
_default_project[] = RUO2_PROJECT
