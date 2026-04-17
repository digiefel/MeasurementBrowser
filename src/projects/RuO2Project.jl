"""
RuO2Project.jl - Project dispatch methods for RuO2test ferroelectric measurements
"""

using Dates
include("RuO2/Display.jl")
include("RuO2/Interpretation.jl")
include("RuO2/PlotHelpers.jl")
include("RuO2/PlotSingle.jl")
include("RuO2/PlotCombined.jl")

# ---------------------------------------------------------------------------
# Regex patterns (RuO2-specific)
# ---------------------------------------------------------------------------
const REGEX_RUO2_CHIP_DIR = r"^RuO2[^/]*_[^/]*$"
const REGEX_RUO2_IDENTIFIER = r"^((?:RuO2)[^()\[\];\s]+)$"
const REGEX_RUO2_BRACKET_IDENTIFIER = r"\[((?:RuO2)[^()\[\];\s]+)"

# ---------------------------------------------------------------------------
# Interface implementations
# ---------------------------------------------------------------------------

project_name(::RuO2Project) = "RuO2"
project_description(::RuO2Project) = "Ferroelectric RuO2test measurements"

function expand_measurement(::RuO2Project, meas::MeasurementInfo)::Vector{MeasurementInfo}
    expanded = _ruo2_expand_multi_device(meas)
    return vcat([_ruo2_expand_pund_fatigue(m) for m in expanded]...)
end

# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------

function _ruo2_expand_multi_device(meas::MeasurementInfo)::Vector{MeasurementInfo}
    meas.measurement_kind == :breakdown || return [meas]
    dev = last(meas.device_info.location)
    if (m = match(r"^([A-Z][0-9]+)([A-Z][0-9]+)$", dev)) === nothing
        return [meas]
    end
    parts = m.captures
    loc = copy(meas.device_info.location)
    return [MeasurementInfo(
        item_id(file_id(meas.filepath); split=p),
        meas.filename,
        meas.filepath,
        replace(meas.clean_title, dev => p),
        meas.measurement_kind,
        meas.timestamp,
        DeviceInfo(vcat(loc[1:end-1], [p]), copy(meas.device_info.parameters)),
        copy(meas.parameters),
        meas.wakeup_pulse_count,
    ) for p in parts]
end

"""
Scan a PUND fatigue CSV to extract unique cycle numbers and peak voltage amplitude.
Returns (sorted_cycles::Vector{Int}, voltage_V::Union{Float64,Nothing}).
"""
function _ruo2_scan_fatigue_file(filepath::AbstractString; should_cancel::Union{Nothing,Function}=nothing)
    cycles = Int[]
    peak_V = 0.0
    last_cycle = nothing
    saw_header = false

    open(filepath, "r") do io
        for raw_line in eachline(io)
            _check_cancel(should_cancel)

            line = strip(raw_line)
            isempty(line) && continue

            if !saw_header
                line == "Cycle,Time_s,Voltage_V,Current_A" || continue
                saw_header = true
                continue
            end

            cycle_end = findfirst(==(','), line)
            cycle_end === nothing && error("Malformed fatigue row in '$filepath': missing cycle delimiter")

            time_end = findnext(==(','), line, nextind(line, cycle_end))
            time_end === nothing && error("Malformed fatigue row in '$filepath': missing time delimiter")

            voltage_end = findnext(==(','), line, nextind(line, time_end))
            voltage_end === nothing && error("Malformed fatigue row in '$filepath': missing voltage delimiter")

            cycle = parse(Int, SubString(line, firstindex(line), prevind(line, cycle_end)))
            voltage = abs(parse(Float64, SubString(line, nextind(line, time_end), prevind(line, voltage_end))))

            if last_cycle === nothing
                push!(cycles, cycle)
                last_cycle = cycle
            elseif cycle != last_cycle
                cycle > last_cycle || error(
                    "Fatigue file '$filepath' contains non-monotonic cycle blocks: saw cycle $cycle after $last_cycle",
                )
                push!(cycles, cycle)
                last_cycle = cycle
            end

            voltage > peak_V && (peak_V = voltage)
        end
    end

    saw_header || error("Fatigue file '$filepath' is missing the fatigue CSV data header")

    voltage = peak_V > 0 ? round(peak_V; digits=1) : nothing
    return cycles, voltage
end

function _ruo2_expand_pund_fatigue(meas::MeasurementInfo)::Vector{MeasurementInfo}
    meas.measurement_kind == :pund_fatigue || return [meas]
    cycles, voltage_V = try
        _ruo2_scan_fatigue_file(meas.filepath)
    catch e
        @warn "Could not read fatigue cycles from $(meas.filepath)" error = e
        return MeasurementInfo[]
    end
    isempty(cycles) && return MeasurementInfo[]
    expanded = MeasurementInfo[]
    sizehint!(expanded, length(cycles))
    for cycle in cycles
        parameters = copy(meas.parameters)
        parameters[:fatigue_cycle] = cycle
        voltage_V !== nothing && (parameters[:voltage_V] = voltage_V)
        push!(expanded, MeasurementInfo(
            item_id(file_id(meas.filepath); cycle=cycle),
            meas.filename,
            meas.filepath,
            meas.clean_title * " cycle $cycle",
            :pund,
            meas.timestamp,
            DeviceInfo(copy(meas.device_info.location), copy(meas.device_info.parameters)),
            parameters,
            nothing,
        ))
    end
    return expanded
end

function combined_plot_types(::RuO2Project)
    return [
        (nothing, "None", "No combined plot selected"),
        (:tlm_analysis, "TLM Analysis", "Width-normalized resistance vs length from multiple TLM 4-point measurements"),
        (:tlm_temperature, "TLM vs Temperature", "Sheet resistance vs temperature (groups by site/chip)"),
        (:pund_fatigue, "PUND Fatigue", "P-E curve evolution or remnant polarization vs fatigue cycles from PUND measurements"),
    ]
end

function compatible_kinds(::RuO2Project, combined_kind::Symbol)
    if combined_kind === :tlm_analysis || combined_kind === :tlm_temperature
        return [:tlm4p]
    elseif combined_kind === :pund_fatigue
        return [:pund, :wakeup]
    end
    return Symbol[]
end

# ---------------------------------------------------------------------------
# Registration
# ---------------------------------------------------------------------------

const RUO2_PROJECT = RuO2Project()
push!(KNOWN_PROJECTS, RUO2_PROJECT)
_default_project[] = RUO2_PROJECT
