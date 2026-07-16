using DataFrames: DataFrame, nrow
using GLMakie
using Statistics: mean

const SI_PREFIXES = Dict(
    -24 => "y",
    -21 => "z",
    -18 => "a",
    -15 => "f",
    -12 => "p",
    -9 => "n",
    -6 => "μ",
    -3 => "m",
    0 => "",
    3 => "k",
    6 => "M",
    9 => "G",
    12 => "T",
    15 => "P",
    18 => "E",
    21 => "Z",
    24 => "Y",
)

function engineering_label(value::Real, unit::AbstractString; sigdigits::Int=4)::String
    isfinite(value) || return "$(value) $unit"
    value == 0 && return "0 $unit"

    exponent = clamp(3 * floor(Int, log10(abs(value)) / 3), -24, 24)
    scaled = value / 10.0^exponent
    prefix = get(SI_PREFIXES, exponent, "")
    return "$(round(scaled; sigdigits)) $(prefix)$unit"
end

function with_plot_units(df::DataFrame, measurement)::DataFrame
    out = DataFrame(df)
    nrow(out) == 0 && return out

    hasproperty(out, :time) && (out.time_us = out.time .* 1e6)
    hasproperty(out, :current) && (out.current_uA = out.current .* 1e6)
    hasproperty(out, :I_FE) && (out.i_fe_uA = out.I_FE .* 1e6)
    if hasproperty(out, :Q_FE)
        finite_q = filter(isfinite, out.Q_FE)
        q_center = isempty(finite_q) ? 0.0 : mean(finite_q)
        q_centered = out.Q_FE .- q_center
        out.q_centered_pC = q_centered .* 1e12

        area_um2 = get(metadata(measurement), :area_um2, nothing)
        if area_um2 !== nothing && Float64(area_um2) > 0
            out.p_fe_uC_cm2 = (q_centered ./ (Float64(area_um2) / 1e8)) .* 1e6
        end
    end
    return out
end

function pulse_group_size(df::DataFrame)::Int
    hasproperty(df, :pulse_idx) || return 1
    max_pulse = maximum(df.pulse_idx)
    max_pulse <= 1 && return 1
    max_pulse % 5 == 0 && return 5
    max_pulse % 4 == 0 && return 4
    max_pulse % 2 == 0 && return 2
    return max_pulse
end

function pulse_groups(df::DataFrame)::Vector{Tuple{Int,BitVector}}
    hasproperty(df, :pulse_idx) || return Tuple{Int,BitVector}[]
    isempty(df.pulse_idx) && return Tuple{Int,BitVector}[]

    group_size = pulse_group_size(df)
    max_pulse = maximum(df.pulse_idx)
    groups = Tuple{Int,BitVector}[]
    for rep in 1:(max_pulse ÷ group_size)
        pulse_range = ((rep - 1) * group_size + 1):(rep * group_size)
        mask = BitVector([pulse in pulse_range for pulse in df.pulse_idx])
        any(mask) && push!(groups, (rep, mask))
    end
    return groups
end
