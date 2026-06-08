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

available_plot_kinds(::RuO2Project)::Vector{Type{<:PlotKind}} = [
    RuO2PUNDPlot,
    RuO2IVSweepPlot,
    RuO2TLM4PointPlot,
    RuO2CVSweepPlot,
    RuO2TLMAnalysisPlot,
    RuO2TLMTemperaturePlot,
    RuO2PUNDFatiguePlot,
]

plot_kind_label(::Type{RuO2PUNDPlot})::String = "PUND"
plot_kind_label(::Type{RuO2IVSweepPlot})::String = "I-V Sweep"
plot_kind_label(::Type{RuO2TLM4PointPlot})::String = "TLM 4-Point"
plot_kind_label(::Type{RuO2CVSweepPlot})::String = "C-V Sweep"
plot_kind_label(::Type{RuO2TLMAnalysisPlot})::String = "TLM Analysis"
plot_kind_label(::Type{RuO2TLMTemperaturePlot})::String = "TLM vs Temperature"
plot_kind_label(::Type{RuO2PUNDFatiguePlot})::String = "PUND Fatigue"

plot_kind_description(::Type{RuO2PUNDPlot})::String = "PUND, PN, and wakeup waveforms"
plot_kind_description(::Type{RuO2IVSweepPlot})::String = "I-V and breakdown sweeps"
plot_kind_description(::Type{RuO2TLM4PointPlot})::String = "Single TLM 4-point measurement"
plot_kind_description(::Type{RuO2CVSweepPlot})::String = "C-V sweep overlay"
plot_kind_description(::Type{RuO2TLMAnalysisPlot})::String = "Width-normalized resistance vs length from multiple TLM measurements"
plot_kind_description(::Type{RuO2TLMTemperaturePlot})::String = "Sheet resistance or resistivity vs temperature"
plot_kind_description(::Type{RuO2PUNDFatiguePlot})::String = "P-E curve evolution and remnant polarization vs fatigue cycles"

plot_kind_measurement_kinds(::Type{RuO2PUNDPlot})::Vector{Symbol} =
    [:pund, :pn, :wakeup_pn, :wakeup_pund]
plot_kind_measurement_kinds(::Type{RuO2IVSweepPlot})::Vector{Symbol} =
    [:iv, :breakdown, :unknown]
plot_kind_measurement_kinds(::Type{RuO2TLM4PointPlot})::Vector{Symbol} = [:tlm4p]
plot_kind_measurement_kinds(::Type{RuO2CVSweepPlot})::Vector{Symbol} = [:cvsweep]
plot_kind_measurement_kinds(::Type{RuO2TLMAnalysisPlot})::Vector{Symbol} = [:tlm4p]
plot_kind_measurement_kinds(::Type{RuO2TLMTemperaturePlot})::Vector{Symbol} = [:tlm4p]
plot_kind_measurement_kinds(::Type{RuO2PUNDFatiguePlot})::Vector{Symbol} = [:pund]

plot_kind_min_measurements(::Type{RuO2TLMAnalysisPlot})::Int = 2
plot_kind_min_measurements(::Type{RuO2TLMTemperaturePlot})::Int = 2
plot_kind_min_measurements(::Type{RuO2PUNDFatiguePlot})::Int = 2

supports_plot_kind(::Type{RuO2PUNDFatiguePlot}, measurement::MeasurementInfo)::Bool =
    measurement.measurement_kind === :pund && haskey(measurement.parameters, :fatigue_idx)

function default_plot_kind(::RuO2Project, measurement::MeasurementInfo)::PlotKindSelection
    measurement.measurement_kind in (:pund, :pn, :wakeup_pn, :wakeup_pund) && return RuO2PUNDPlot
    measurement.measurement_kind in (:iv, :breakdown, :unknown) && return RuO2IVSweepPlot
    measurement.measurement_kind === :tlm4p && return RuO2TLM4PointPlot
    measurement.measurement_kind === :cvsweep && return RuO2CVSweepPlot
    return nothing
end

# ---------------------------------------------------------------------------
# Registration
# ---------------------------------------------------------------------------

const RUO2_PROJECT = RuO2Project()
push!(KNOWN_PROJECTS, RUO2_PROJECT)
_default_project[] = RUO2_PROJECT
