"""
TASEProject.jl - Project dispatch methods for TASE four-terminal IV measurements

Filename pattern:
  {Chip}_{Facet}_{DeviceType}_{DeviceID}_{YYYYMMDD}_{HHMMSS}_{T}K_FourTerminalIV.csv

Example:
  TASESNS1c1f_A_2TSNJunction_11_20260224_111623_298K_FourTerminalIV.csv
"""

using GLMakie
include("TASE/Display.jl")
include("TASE/Interpretation.jl")
include("TASE/Analysis.jl")

struct TASEFourTerminalIVPlot <: PlotKind end

include("TASE/Plotting.jl")

const REGEX_TASE = r"^([^_]+)_([^_]+)_([^_]+)_(\d+)_\d{8}_\d{6}_\d+K_FourTerminalIV\.csv"i

# ---------------------------------------------------------------------------
# Interface implementations
# ---------------------------------------------------------------------------

project_name(::TASEProject) = "TASE"
project_description(::TASEProject) = "GaN TASE four-terminal IV"

available_plot_kinds(::TASEProject)::Vector{Type{<:PlotKind}} = [TASEFourTerminalIVPlot]
plot_kind_label(::Type{TASEFourTerminalIVPlot})::String = "Four-terminal IV"
plot_kind_description(::Type{TASEFourTerminalIVPlot})::String = "Four-terminal IV overlay"
plot_kind_measurement_kinds(::Type{TASEFourTerminalIVPlot})::Vector{Symbol} = [:four_terminal_iv]
default_plot_kind(::TASEProject, measurement::MeasurementInfo)::PlotKindSelection =
    measurement.measurement_kind === :four_terminal_iv ? TASEFourTerminalIVPlot : nothing

# ---------------------------------------------------------------------------
# Registration
# ---------------------------------------------------------------------------

const TASE_PROJECT = TASEProject()
push!(KNOWN_PROJECTS, TASE_PROJECT)
