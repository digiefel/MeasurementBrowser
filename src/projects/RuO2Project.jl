"""
RuO2Project.jl - Project dispatch methods for RuO2test ferroelectric measurements
"""

using Dates
using DataFrames
using DataPlotter: load_fe_pund_single, load_iv_sweep_single, load_tlm_4p_single, load_wakeup_single,
                   plot_fe_pund, plot_iv_sweep_single, plot_tlm_4p, plot_wakeup,
                   plot_tlm_combined, plot_tlm_temperature, plot_pund_fatigue

# ---------------------------------------------------------------------------
# Regex patterns (RuO2-specific)
# ---------------------------------------------------------------------------
const REGEX_DEVICE = r"RuO2test_([A-Z0-9]+)_([A-Z0-9]+)_([A-Z0-9]+(?:W[0-9]+)?)"
const REGEX_DEVICE_NEW = r"^RuO2test_([^_]+)_([^_]+)_([^_]+)_([^_]+)_"

# ---------------------------------------------------------------------------
# Interface implementations
# ---------------------------------------------------------------------------

project_name(::RuO2Project) = "RuO2"
project_description(::RuO2Project) = "Ferroelectric RuO2test measurements"

accepts_file(::RuO2Project, filename::String) = occursin("ruo2test_", lowercase(filename))

function parse_device_info(::RuO2Project, filename::String)
    if (m = match(REGEX_DEVICE_NEW, filename)) !== nothing
        caps = filter(!isnothing, collect(m.captures))
        return DeviceInfo(String.(caps))
    elseif (m = match(REGEX_DEVICE, filename)) !== nothing
        caps = filter(!isnothing, collect(m.captures))
        return DeviceInfo(String.(caps))
    end
    error("Unrecognized RuO2 device filename format: $filename")
end

function detect_kind(::RuO2Project, filename::String)::Symbol
    lower = lowercase(filename)
    if occursin("pund_fatigue", lower) || occursin("pund fatigue", lower)
        return :pund_fatigue
    elseif occursin("fe pund", lower) || occursin("fepund", lower)
        return :pund
    elseif occursin("i_v sweep", lower) || occursin("iv sweep", lower)
        return :iv
    elseif occursin("tlm_4p", lower) || occursin("tlm", lower)
        return :tlm4p
    elseif occursin("break", lower) || occursin("breakdown", lower)
        return :breakdown
    elseif occursin("wakeup", lower)
        return :wakeup
    else
        return :unknown
    end
end

function kind_label(::RuO2Project, kind::Symbol)::String
    kind === :pund && return "FE PUND"
    kind === :pund_fatigue && return "PUND Fatigue"
    kind === :iv && return "I-V Sweep"
    kind === :tlm4p && return "TLM 4-Point"
    kind === :breakdown && return "Breakdown"
    kind === :wakeup && return "Wakeup"
    return "Unknown"
end

function display_label(proj::RuO2Project, meas::MeasurementInfo)
    label = kind_label(proj, meas.measurement_kind)
    temp = get(meas.parameters, :temperature_K, nothing)
    parts = Any[meas.timestamp, label]

    if meas.measurement_kind == :wakeup
        meas.wakeup_pulse_count !== nothing && push!(parts, "$(meas.wakeup_pulse_count)×")
    elseif meas.measurement_kind == :pund
        m = match(r"(\d+(?:\.\d+)?)V", meas.filename)
        voltage = m !== nothing ? tryparse(Float64, m.captures[1]) : get(meas.parameters, :voltage_V, nothing)
        voltage !== nothing && push!(parts, "$(voltage)V")
        haskey(meas.parameters, :fatigue_cycle) && push!(parts, "cycle $(meas.parameters[:fatigue_cycle]) (fatigue)")
    end

    temp !== nothing && push!(parts, "$(temp)K")
    return join(parts, " ")
end

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
        meas.filename,
        meas.filepath,
        replace(meas.clean_title, dev => p),
        meas.measurement_kind,
        meas.timestamp,
        DeviceInfo(vcat(loc[1:end-1], [p]), deepcopy(meas.device_info.parameters)),
        deepcopy(meas.parameters),
        meas.wakeup_pulse_count,
    ) for p in parts]
end

"""
Scan a PUND fatigue CSV to extract unique cycle numbers and peak voltage amplitude.
Returns (sorted_cycles::Vector{Int}, voltage_V::Union{Float64,Nothing}).
"""
function _ruo2_scan_fatigue_file(filepath::AbstractString)
    cycles = Set{Int}()
    peak_V = 0.0
    open(filepath, "r") do io
        readline(io)  # skip header
        for line in eachline(io)
            isempty(line) && continue
            idx1 = findfirst(',', line)
            idx1 === nothing && continue
            idx2 = findnext(',', line, idx1 + 1)
            idx2 === nothing && continue
            idx3 = findnext(',', line, idx2 + 1)
            v_end = idx3 === nothing ? lastindex(line) : idx3 - 1
            try
                push!(cycles, parse(Int, @view line[1:idx1-1]))
                v = abs(parse(Float64, @view line[idx2+1:v_end]))
                v > peak_V && (peak_V = v)
            catch
                continue
            end
        end
    end
    voltage = peak_V > 0 ? round(peak_V; digits=1) : nothing
    return sort!(collect(cycles)), voltage
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
    extra = Dict{Symbol,Any}()
    voltage_V !== nothing && (extra[:voltage_V] = voltage_V)
    return [MeasurementInfo(
        meas.filename,
        meas.filepath,
        meas.clean_title * " cycle $c",
        :pund,
        meas.timestamp,
        DeviceInfo(copy(meas.device_info.location), deepcopy(meas.device_info.parameters)),
        merge(deepcopy(meas.parameters), extra, Dict{Symbol,Any}(:fatigue_cycle => c)),
        nothing,
    ) for c in cycles]
end

# ---------------------------------------------------------------------------
# Plot dispatch
# ---------------------------------------------------------------------------

function load_plot_input_for_file(::RuO2Project, path::AbstractString, kind::Union{Symbol,Nothing}; kwargs...)
    device_params = get(kwargs, :device_params, nothing)
    if kind === :pund
        return load_fe_pund_single(path; device_params)
    elseif kind === :iv
        return load_iv_sweep_single(path)
    elseif kind === :tlm4p
        return load_tlm_4p_single(path)
    elseif kind === :breakdown
        loaded = load_iv_sweep_single(path)
        loaded === nothing && return nothing
        return (df=loaded.df, title=loaded.title * " (Breakdown)")
    elseif kind === :wakeup
        return load_wakeup_single(path)
    else
        # Keep historical fallback behavior for unknown kinds.
        @warn "figure_for_file: Invalid measurement type. Attempting I-V sweep reader."
        return load_iv_sweep_single(path)
    end
end

function draw_plot_from_input(::RuO2Project, kind::Union{Symbol,Nothing}, loaded; kwargs...)
    loaded === nothing && return nothing

    if kind === :pund
        area_um2 = hasproperty(loaded, :area_um2) ? loaded.area_um2 : nothing
        return plot_fe_pund(loaded.df, loaded.title; area_um2, kwargs...)
    elseif kind === :iv || kind === :breakdown || kind === :unknown || kind === nothing
        return plot_iv_sweep_single(loaded.df, loaded.title; kwargs...)
    elseif kind === :tlm4p
        return plot_tlm_4p(loaded.df, loaded.title; kwargs...)
    elseif kind === :wakeup
        return plot_wakeup(loaded.df, loaded.title; kwargs...)
    end

    return nothing
end

function figure_for_file(proj::RuO2Project, path::AbstractString, kind::Union{Symbol,Nothing}; kwargs...)
    try
        loaded = load_plot_input_for_file(proj, path, kind; kwargs...)
        return draw_plot_from_input(proj, kind, loaded; kwargs...)
    catch err
        @warn "figure_for_file failed" path kind error = err
        return nothing
    end
end

function figure_for_files(::RuO2Project, paths::Vector{String}, combined_kind::Symbol;
                          device_params_list::Vector{Dict{Symbol,Any}}=Dict{Symbol,Any}[], kwargs...)
    isempty(paths) && return nothing
    try
        if combined_kind === :tlm_analysis
            return plot_tlm_combined(paths; device_params_list=device_params_list, kwargs...)
        elseif combined_kind === :tlm_temperature
            return plot_tlm_temperature(paths; device_params_list=device_params_list, kwargs...)
        elseif combined_kind === :pund_fatigue
            return plot_pund_fatigue(paths; device_params_list=device_params_list, kwargs...)
        else
            @warn "figure_for_files: Unknown combined plot kind: $combined_kind"
            return nothing
        end
    catch err
        @warn "figure_for_files failed" paths combined_kind error = err
        return nothing
    end
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
