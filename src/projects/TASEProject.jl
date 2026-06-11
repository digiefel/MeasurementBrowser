"""
TASEProject.jl - Project dispatch methods for TASE four-terminal IV measurements

Filename pattern:
  {Chip}_{Facet}_{DeviceType}_{DeviceID}_{YYYYMMDD}_{HHMMSS}_{T}K_FourTerminalIV.csv

Example:
  TASESNS1c1f_A_2TSNJunction_11_20260224_111623_298K_FourTerminalIV.csv
"""

using GLMakie

struct TASEProject <: AbstractProject end

include("TASE/Display.jl")
include("TASE/Interpretation.jl")

struct TASEFourTerminalIVPlot <: PlotKind end

include("TASE/Plotting.jl")

const REGEX_TASE = r"^([^_]+)_([^_]+)_([^_]+)_(\d+)_\d{8}_\d{6}_\d+K_FourTerminalIV\.csv"i

# ---------------------------------------------------------------------------
# Interface implementations
# ---------------------------------------------------------------------------

project_name(::TASEProject) = "TASE"
project_description(::TASEProject) = "GaN TASE four-terminal IV"

# ---------------------------------------------------------------------------
# Registration
# ---------------------------------------------------------------------------

const TASE_PROJECT = TASEProject()
push!(PROJECTS, TASE_PROJECT)
