using DataLoader: read_fe_pund, read_pund_wakeup_amplitude
using DataAnalysis: detect_pund_pulses, analyze_pund
using DataFrames: nrow
using Statistics: mean

"""
Stats attached to one PUND logical measurement. These describe the selected waveform,
not the voltage/frequency settings requested in the file header.
"""
function pund_stats_from_waveform(df, device_params::Dict{Symbol,Any}=Dict{Symbol,Any}())::Dict{Symbol,Any}
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
        :frequency_kHz => NaN,
        :Pr_max_uCcm2 => NaN,
    )

    det = detect_pund_pulses(df.time, df.voltage, df.current)
    if !isempty(det.pulses)
        t_start = df.time[first(det.pulses[1])]
        t_end = df.time[last(det.pulses[end])]
        total_span = t_end - t_start
        if total_span > 0
            pulse_period = total_span / length(det.pulses)
            stats[:frequency_kHz] = round_one(1.0 / (pulse_period * 1e3))
        end
    end

    if haskey(device_params, :area_um2)
        area_um2 = Float64(device_params[:area_um2])
        if isfinite(area_um2) && area_um2 > 0
            df_an = hasproperty(df, :Q_FE) && hasproperty(df, :pulse_idx) ?
                df : analyze_pund(df; pulses=det.pulses)
            pr_vals = pund_pr_values(df_an, area_um2)
            isempty(pr_vals) || (stats[:Pr_max_uCcm2] = round(maximum(pr_vals); digits=3))
        end
    end

    return stats
end

"""
File-level entry point for PUND stats. `measurement_params` selects the standalone,
wakeup, or fatigue waveform before stats are computed from that selected waveform.
"""
function compute_pund_stats(
    filepath::AbstractString,
    measurement_params::Dict{Symbol,Any},
    device_params::Dict{Symbol,Any}=Dict{Symbol,Any}(),
)::Dict{Symbol,Any}
    fname = basename(filepath)
    dir = dirname(filepath)
    fatigue_count = measurement_params[:fatigue_count]
    wakeup_V = measurement_params[:wakeup_V]

    if fatigue_count > 0
        fatigue_df = _load_ruo2_pund_fatigue_file(filepath)
        return pund_stats_from_waveform(
            _select_pund_fatigue_cycle(fatigue_df, Int(fatigue_count)),
            device_params,
        )
    elseif isfinite(Float64(wakeup_V))
        df = read_pund_wakeup_amplitude(fname, dir, Float64(wakeup_V), 1)
        return pund_stats_from_waveform(df, device_params)
    end

    return pund_stats_from_waveform(read_fe_pund(fname, dir), device_params)
end

"""
Convert analyzed switching charge into remnant polarization per repetition.
"""
function pund_pr_values(df, area_um2::Float64)
    pid = df.pulse_idx
    maxpid = maximum(pid)
    group_size = maxpid == 1 ? 1 : maxpid % 5 == 0 ? 5 : maxpid % 4 == 0 ? 4 : 2
    n_groups = maxpid ÷ group_size
    area_cm2 = area_um2 / 1e8

    pr_vals = Float64[]
    for rep in 1:n_groups
        base = (rep - 1) * group_size
        rep_mask = if group_size == 1
            pid .== base + 1
        else
            P_code = group_size == 5 ? base + 2 : base + 1
            N_code = group_size == 5 ? base + 4 : group_size == 4 ? base + 3 : base + 2
            (pid .== P_code) .| (pid .== N_code)
        end

        Q_FE = df.Q_FE
        valid = rep_mask .& isfinite.(Q_FE) .& .!isnan.(Q_FE)
        if any(valid)
            Qc = Q_FE[valid] .- mean(Q_FE[valid])
            P_FE = (Qc ./ area_cm2) .* 1e6
            ymin = minimum(P_FE)
            ymax = maximum(P_FE)
            if isfinite(ymin) && isfinite(ymax)
                push!(pr_vals, 0.5 * (ymax - ymin))
            end
        end
    end
    return pr_vals
end
