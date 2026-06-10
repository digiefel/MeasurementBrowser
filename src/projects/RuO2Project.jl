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

# ---------------------------------------------------------------------------
# Registration
# ---------------------------------------------------------------------------

const RUO2_PROJECT = RuO2Project()
push!(KNOWN_PROJECTS, RUO2_PROJECT)
_default_project[] = RUO2_PROJECT
