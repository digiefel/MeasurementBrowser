"""
RuO2Project.jl - Project dispatch methods for RuO2test ferroelectric measurements
"""

using Dates
using DataFrames: DataFrame
using DataLoader: read_iv_sweep, read_tlm_4p

struct RuO2Project <: AbstractProject end

include("RuO2/Display.jl")
include("RuO2/CVSweepIO.jl")
include("RuO2/Interpretation.jl")
include("RuO2/PUNDIO.jl")
include("RuO2/PUNDAnalysis.jl")
include("RuO2/Stats.jl")

"""Load the direct table for one RuO2 source file or logical measurement."""
function load_source_data(
    ::RuO2Project,
    source_file::SourceFile;
    measurement::Union{Nothing,MeasurementInfo}=nothing,
)::DataFrame
    check_cancel()
    kind = measurement === nothing ?
        detect_kind(RUO2_PROJECT, source_file.filename) :
        measurement.measurement_kind
    if kind === :pund || kind === :pn
        if is_pund_fatigue_file(source_file.filepath)
            measurement === nothing && error("PUND fatigue data requires a logical measurement")
            fatigue_count = Int(measurement.parameters[:fatigue_idx])
            fatigue_df = read_pund_fatigue_file(source_file.filepath)
            return _select_pund_fatigue_cycle(fatigue_df, fatigue_count)
        end
        return read_pund_file(source_file.filename, dirname(source_file.filepath))
    elseif kind === :wakeup_pn || kind === :wakeup_pund
        measurement === nothing && error("PUND wakeup data requires a logical measurement")
        wakeup_V = Float64(measurement.parameters[:wakeup_V])
        return _select_pund_wakeup_readout(source_file.filepath, wakeup_V)
    elseif kind === :iv || kind === :breakdown
        return read_iv_sweep(source_file.filename, dirname(source_file.filepath))
    elseif kind === :tlm4p
        return read_tlm_4p(source_file.filename, dirname(source_file.filepath))
    elseif kind === :cvsweep
        return read_cv_sweep(source_file.filename, dirname(source_file.filepath))
    end
    error("RuO2 source data API is not implemented for $kind")
end

project_name(::RuO2Project) = "RuO2"
project_description(::RuO2Project) = "Ferroelectric RuO2test measurements"

const RUO2_PROJECT = RuO2Project()
push!(PROJECTS, RUO2_PROJECT)
DEFAULT_PROJECT[] = RUO2_PROJECT
