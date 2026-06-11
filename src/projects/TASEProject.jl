"""
TASEProject.jl - Project dispatch methods for TASE four-terminal IV measurements

Filename pattern:
  {Chip}_{Facet}_{DeviceType}_{DeviceID}_{YYYYMMDD}_{HHMMSS}_{T}K_FourTerminalIV.csv

Example:
  TASESNS1c1f_A_2TSNJunction_11_20260224_111623_298K_FourTerminalIV.csv
"""

using GLMakie
using DataFrames: DataFrame
using DataLoader: read_tlm_4p

struct TASEProject <: AbstractProject end

include("TASE/Display.jl")
include("TASE/Interpretation.jl")

"""Load the direct four-terminal IV table for one TASE source file."""
function load_source_data(
    ::TASEProject,
    source_file::SourceFile;
    measurement::Union{Nothing,MeasurementInfo}=nothing,
)::DataFrame
    check_cancel()
    if measurement !== nothing && measurement.measurement_kind !== :four_terminal_iv
        error("TASE cannot load measurement data for $(measurement.measurement_kind)")
    end
    detect_kind(TASE_PROJECT, source_file.filename) === :four_terminal_iv ||
        error("TASE cannot load source data for $(source_file.filepath)")
    return read_tlm_4p(source_file.filename, dirname(source_file.filepath))
end

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
