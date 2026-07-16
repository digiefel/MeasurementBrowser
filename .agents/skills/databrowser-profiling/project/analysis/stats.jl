using DataFrames: AbstractDataFrame, DataFrame, nrow
using Statistics: mean

# ---------------------------------------------------------------------------
# PUND waveform stats
# ---------------------------------------------------------------------------

"""Per-item PUND stats: voltage envelope, pulse frequency, and Pr when an area is known."""
function pund_stats(data::AbstractDataFrame, metadata)::Dict{Symbol,Any}
    nrow(data) == 0 && return Dict{Symbol,Any}()
    return pund_stats_from_waveform(data, metadata)
end

function pund_stats_from_waveform(
    df::AbstractDataFrame,
    device_params::AbstractDict,
)::Dict{Symbol,Any}
    nrow(df) > 0 || error("Cannot compute PUND stats from an empty dataframe")
    V = df.voltage
    n_pre = min(9, length(V))
    round_one(x) = (y = round(x; digits=1); iszero(y) ? 0.0 : y)
    V_base = round_one(mean(V[1:n_pre]))
    V_min = round_one(minimum(V))
    V_max = round_one(maximum(V))
    stats = Dict{Symbol,Any}(
        :V_base => V_base,
        :V_min => V_min,
        :V_max => V_max,
        :V_amp => round_one(max(abs(V_max - V_base), abs(V_min - V_base))),
    )

    pulses = pund_pulses(df)
    frequency_kHz = pund_frequency_kHz(df, pulses)
    isfinite(frequency_kHz) && (stats[:frequency_kHz] = frequency_kHz)

    area_um2 = get(device_params, :area_um2, nothing)
    if area_um2 !== nothing && isfinite(Float64(area_um2)) && Float64(area_um2) > 0
        analyzed = (hasproperty(df, :Q_FE) && hasproperty(df, :pulse_idx)) ? df : analyze_pund(df; pulses)
        stats[:Pr_uCcm2] = round(pund_pr_value(analyzed, Float64(area_um2)); digits=3)
    end
    return stats
end

function pund_pulses(df::AbstractDataFrame)::Vector{UnitRange{Int}}
    return hasproperty(df, :pulse_idx) ? _pulse_ranges_from_indices(df.pulse_idx) :
        detect_pund_pulses(df.time, df.voltage, df.current).pulses
end

function pund_frequency_kHz(
    df::AbstractDataFrame,
    pulses::Vector{UnitRange{Int}}=pund_pulses(df),
)::Float64
    isempty(pulses) && return NaN
    span = df.time[last(pulses[end])] - df.time[first(pulses[1])]
    span > 0 || return NaN
    value = round(1.0 / ((span / length(pulses)) * 1e3); digits=1)
    return iszero(value) ? 0.0 : value
end

function _pulse_ranges_from_indices(pulse_idx)::Vector{UnitRange{Int}}
    ranges = UnitRange{Int}[]
    start = nothing
    current = 0
    for i in eachindex(pulse_idx)
        pulse = pulse_idx[i]
        if pulse > 0 && pulse != current
            start !== nothing && push!(ranges, start:(i - 1))
            start = i
            current = pulse
        elseif pulse == 0 && start !== nothing
            push!(ranges, start:(i - 1))
            start = nothing
            current = 0
        end
    end
    start !== nothing && push!(ranges, start:lastindex(pulse_idx))
    return ranges
end

"""Mean remnant polarization (µC/cm²) over the valid PUND loop(s)."""
function pund_pr_value(df::DataFrame, area_um2::Float64)::Float64
    pid = df.pulse_idx
    maxpid = maximum(pid)
    group_size = maxpid == 1 ? 1 : maxpid % 5 == 0 ? 5 : maxpid % 4 == 0 ? 4 : 2
    n_groups = maxpid ÷ group_size
    area_cm2 = area_um2 / 1e8
    pr_vals = Float64[]
    for rep in 1:n_groups
        base = (rep - 1) * group_size
        mask = if group_size == 1
            pid .== base + 1
        else
            P_code = group_size == 5 ? base + 2 : base + 1
            N_code = group_size == 5 ? base + 4 : group_size == 4 ? base + 3 : base + 2
            (pid .== P_code) .| (pid .== N_code)
        end
        Q = df.Q_FE
        valid = mask .& isfinite.(Q)
        if any(valid)
            P = ((Q[valid] .- mean(Q[valid])) ./ area_cm2) .* 1e6
            push!(pr_vals, 0.5 * (maximum(P) - minimum(P)))
        end
    end
    isempty(pr_vals) && error("Expected at least one PUND loop")
    return mean(pr_vals)
end

# ---------------------------------------------------------------------------
# Current-voltage processing
# ---------------------------------------------------------------------------

"""Add resistance columns and a linear R fit to a current-voltage table."""
function iv_resistance(data::DataFrame)::DataFrame
    mask = isfinite.(data.i) .& isfinite.(data.v)
    i = Float64.(data.i[mask])
    v = Float64.(data.v[mask])
    slope, offset = _ols(i, v)
    resistance = v ./ i
    return DataFrame(
        i=i,
        v=v,
        resistance_ohm=resistance,
        valid_resistance=isfinite.(resistance) .& (abs.(resistance) .< 1e9),
        fit_resistance_ohm=fill(slope, length(i)),
        fit_offset_v=fill(offset, length(i)),
    )
end

iv_stats(data::DataFrame, _metadata)::Dict{Symbol,Any} =
    nrow(data) == 0 ? Dict{Symbol,Any}() :
    Dict{Symbol,Any}(:resistance_ohm => first(data.fit_resistance_ohm))

function _ols(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})::Tuple{Float64,Float64}
    n = length(x)
    n == 0 && return (NaN, NaN)
    sx = sum(x); sy = sum(y); sxx = sum(x .^ 2); sxy = sum(x .* y)
    denom = n * sxx - sx^2
    abs(denom) < 1e-20 && return (NaN, NaN)
    slope = (n * sxy - sx * sy) / denom
    return (slope, (sy - slope * sx) / n)
end

# ---------------------------------------------------------------------------
# Expansion helpers (one file -> several measurements)
# ---------------------------------------------------------------------------

"""Read a wakeup file and attach header values needed by every expanded item."""
function read_wakeup(file)
    data = read_pund_wakeup_file(file.filepath)
    header = read_header_comments(file.filepath)
    pulse_type = lowercase(get(header, "read_pulse_type", "both"))
    metadata = merge(
        ruo2_file_metadata(file),
        Dict{Symbol,Any}(
            :read_pulse_type => pulse_type,
            :wakeup_count => _header_float(header, "fatigue_count"),
            :wakeup_f => _header_float(header, "fatigue_freq"),
        ),
    )
    return (data=data, metadata=metadata)
end

"""A wakeup file holds one readout per configured voltage and PN/PUND segment."""
function wakeup_measurements(data::DataFrame, metadata)::Vector
    pulse_type = metadata[:read_pulse_type]
    segments = pulse_type == "pund" ? [:pund] : pulse_type == "pn" ? [:pn] : [:pn, :pund]
    return [
        ruo2_entry(data, metadata; voltage=v, segment=seg)
        for v in sort(unique(data.wakeup_V)) for seg in segments
    ]
end

"""Select one wakeup voltage's readout and return the requested PN or PUND segment."""
function process_wakeup(data::DataFrame, metadata)::DataFrame
    rows = data.wakeup_V .== metadata[:voltage]
    readout = DataFrame(time=data.time[rows], voltage=data.voltage[rows], current=data.current[rows])
    result = analyze_pund_and_pn(readout)
    segment = metadata[:segment] === :pn ? result.pn : result.pund
    return segment === nothing ? DataFrame() : DataFrame(segment)
end

"""Paired-device breakdown files split into one measurement per device."""
function breakdown_measurements(data::DataFrame, metadata)::Vector
    device = String(metadata[:device])
    paired = match(r"^([A-Z][0-9]+)([A-Z][0-9]+)$", device)
    paired === nothing && return [ruo2_entry(data, metadata)]
    return [
        ruo2_entry(
            data,
            metadata;
            device=String(p),
        )
        for p in paired.captures
    ]
end

_header_float(header::AbstractDict, key)::Float64 = parse(Float64, get(header, key, "NaN"))

"""Parse `#key: value` comment lines from the head of a CSV export."""
function read_header_comments(filepath::AbstractString; max_lines::Int=50)::Dict{String,String}
    fields = Dict{String,String}()
    open(filepath, "r") do io
        for (n, line) in enumerate(eachline(io))
            n > max_lines && break
            m = match(r"^#\s*([^:]+):\s*(.*)$", strip(line))
            m === nothing || (fields[strip(m.captures[1])] = strip(m.captures[2]))
        end
    end
    return fields
end
