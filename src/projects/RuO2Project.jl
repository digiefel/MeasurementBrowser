"""
RuO2Project.jl - Project dispatch methods for RuO2test ferroelectric measurements
"""

using Dates
include("RuO2/Display.jl")
include("RuO2/Interpretation.jl")
include("RuO2/PUNDFatigueIO.jl")
include("RuO2/Stats.jl")
include("RuO2/PlotHelpers.jl")
include("RuO2/PlotSingle.jl")
include("RuO2/PlotCombined.jl")

# ---------------------------------------------------------------------------
# Regex patterns (RuO2-specific)
# ---------------------------------------------------------------------------
const REGEX_RUO2_CHIP_DIR = r"^RuO2[^/]*_[^/]*$"
const REGEX_RUO2_IDENTIFIER = r"^((?:RuO2)[^()\[\];\s]+)$"
const REGEX_RUO2_BRACKET_IDENTIFIER = r"\[((?:RuO2)[^()\[\];\s]+)"
const REGEX_RUO2_TIMESTAMP_FILENAME = r"^(RuO2.+?)_\d{8}_\d{6}_"

# ---------------------------------------------------------------------------
# Interface implementations
# ---------------------------------------------------------------------------

project_name(::RuO2Project) = "RuO2"
project_description(::RuO2Project) = "Ferroelectric RuO2test measurements"

"""
Attach global pulse-count metadata after RuO2 measurements have been expanded.
Counts are scoped per physical device so plots and information-panel stats read the same history.
"""
function annotate_measurements!(::RuO2Project, measurements::Vector{MeasurementInfo})
    by_device = Dict{String,Vector{MeasurementInfo}}()
    for measurement in measurements
        delete!(measurement.parameters, :global_pund_pulse_count)
        device_key = device_path_key(measurement.device_info)
        push!(get!(by_device, device_key, MeasurementInfo[]), measurement)
    end

    for device_measurements in values(by_device)
        sort!(device_measurements, by=measurement_timestamp_key)
        pulse_count = 0
        for measurement in device_measurements
            measurement.measurement_kind === :pund || continue
            fatigue_cycle = get(measurement.parameters, :fatigue_cycle, nothing)
            pulse_count = fatigue_cycle === nothing ? pulse_count + 1 : Int(fatigue_cycle)
            measurement.parameters[:global_pund_pulse_count] = pulse_count
        end
    end
    return nothing
end

"""
Add RuO2 summary values that depend on the selected device's full measurement set.
This covers pulse-history totals that are not meaningful for a single measurement row.
"""
function augment_measurements_stats!(
    ::RuO2Project,
    stats::Dict{Symbol,Any},
    measurements::Vector{MeasurementInfo},
)
    global_counts = Int[
        Int(measurement.parameters[:global_pund_pulse_count])
        for measurement in measurements
        if haskey(measurement.parameters, :global_pund_pulse_count)
    ]
    if !isempty(global_counts)
        stats[:global_pund_pulse_count] = maximum(global_counts)
    end

    return stats
end

# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------

"""
Split paired-device breakdown files into per-device measurement entries.

Some RuO2 breakdown runs encode two adjacent devices in one filename; expanding them
preserves device-level navigation and stats without duplicating the source file.
"""
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
        copy(meas.stats),
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
    columns = nothing

    open(filepath, "r") do io
        for raw_line in eachline(io)
            _check_cancel(should_cancel)

            line = strip(raw_line)
            isempty(line) && continue

            if columns === nothing
                if startswith(line, "Cycle,")
                    columns = _ruo2_fatigue_columns(filepath, line)
                end
                continue
            end

            cycle, voltage = _ruo2_parse_fatigue_scan_row(filepath, line, columns)

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

    columns !== nothing || error(
        "Fatigue file '$filepath' is missing a data header with columns Cycle, Time_s, Voltage_V, Current_A",
    )

    voltage = peak_V > 0 ? round(peak_V; digits=1) : nothing
    return cycles, voltage
end

function _ruo2_fatigue_columns(filepath::AbstractString, header::AbstractString)
    columns = strip.(split(header, ','))
    indices = Dict(column => index for (index, column) in pairs(columns))
    required = ("Cycle", "Time_s", "Voltage_V", "Current_A")
    missing = [column for column in required if !haskey(indices, column)]
    isempty(missing) || error(
        "Fatigue file '$filepath' data header is missing required columns: $(join(missing, ", "))",
    )
    return (
        cycle=indices["Cycle"],
        time=indices["Time_s"],
        voltage=indices["Voltage_V"],
        current=indices["Current_A"],
        count=length(columns),
    )
end

function _ruo2_parse_fatigue_scan_row(
    filepath::AbstractString,
    line::AbstractString,
    columns::NamedTuple,
)
    field = 1
    field_start = firstindex(line)
    index = field_start
    cycle = nothing
    voltage = nothing

    while index <= lastindex(line)
        comma = findnext(',', line, index)
        field_stop = comma === nothing ? lastindex(line) : prevind(line, comma)
        if field == columns.cycle
            cycle = parse(Int, SubString(line, field_start, field_stop))
        elseif field == columns.voltage
            voltage = abs(parse(Float64, SubString(line, field_start, field_stop)))
        end

        comma === nothing && break
        field += 1
        index = nextind(line, comma)
        field_start = index
    end

    field == columns.count || error(
        "Malformed fatigue row in '$filepath': expected $(columns.count) columns, got $field",
    )
    cycle !== nothing || error("Fatigue row in '$filepath' is missing Cycle")
    voltage !== nothing || error("Fatigue row in '$filepath' is missing Voltage_V")
    return cycle, voltage
end

"""
List RuO2 combined-plot modes exposed by the project selector.

The returned metadata drives UI labels and constrains downstream combined plotting to
analyses that make sense for RuO2 measurement families.
"""
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
        return [:pund]
    end
    return Symbol[]
end

# ---------------------------------------------------------------------------
# Registration
# ---------------------------------------------------------------------------

const RUO2_PROJECT = RuO2Project()
push!(KNOWN_PROJECTS, RUO2_PROJECT)
_default_project[] = RUO2_PROJECT
