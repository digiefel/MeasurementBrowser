using DataLoader: read_fe_pund, read_pund_wakeup_amplitude, read_pund_wakeup_reps
using DataAnalysis: detect_pund_pulses, analyze_pund, analyze_pund_and_pn
using DataFrames: nrow
using Statistics: mean

"""
    compute_pund_stats(filepath, measurement_params, device_params) -> Dict{Symbol,Any}
    compute_pund_stats_from_analyzed_plot(analyzed, device_params) -> Dict{Symbol,Any}

Compute data-derived statistics for a PUND measurement by reading the CSV.
Returns a Dict with keys:
  :V_base        — DC voltage level, from the first 9 points
  :V_min         — minimum applied voltage
  :V_max         — maximum applied voltage
  :V_amp         — largest excursion away from V_base
  :frequency_kHz — inverse period of a single triangular pulse
  :Pr_max_uCcm2  — maximum remnant polarization across repetitions (µC/cm², if area_um2 available)
"""
function pund_voltage_stats(df)::Dict{Symbol,Any}
    nrow(df) > 0 || error("Cannot compute PUND voltage stats from an empty dataframe")
    V = df.voltage
    n_pre = min(9, length(V))
    one_digit(x) = (y = round(x; digits=1); iszero(y) ? 0.0 : y)
    V_base = one_digit(mean(V[1:n_pre]))
    V_min = one_digit(minimum(V))
    V_max = one_digit(maximum(V))
    return Dict{Symbol,Any}(
        :V_base => V_base,
        :V_min => V_min,
        :V_max => V_max,
        :V_amp => one_digit(max(abs(V_max - V_base), abs(V_min - V_base))),
    )
end

function _pund_stats_from_analysis_df(df, device_params::Dict{Symbol,Any})::Dict{Symbol,Any}
    result = pund_voltage_stats(df)

    t = df.time
    V = df.voltage
    I = df.current

    det = detect_pund_pulses(t, V, I)
    if !isempty(det.pulses)
        t_start = t[first(det.pulses[1])]
        t_end = t[last(det.pulses[end])]
        total_span = t_end - t_start
        if total_span > 0
            pulse_period = total_span / length(det.pulses)
            result[:frequency_kHz] = round(1.0 / (pulse_period * 1e3); digits=1)
        end
    end

    df_an = hasproperty(df, :Q_FE) && hasproperty(df, :pulse_idx) ?
        df : analyze_pund(df; pulses=det.pulses)
    area_um2 = haskey(device_params, :area_um2) ? Float64(device_params[:area_um2]) : NaN

    if isfinite(area_um2) && area_um2 > 0
        pid = df_an.pulse_idx
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

            Q_FE = df_an.Q_FE
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

        if !isempty(pr_vals)
            result[:Pr_max_uCcm2] = round(maximum(pr_vals); digits=3)
        end
    end

    return result
end

function compute_pund_stats_from_analyzed_plot(
    analyzed,
    device_params::Dict{Symbol,Any},
)::Dict{Symbol,Any}
    hasproperty(analyzed, :df) || error("Analyzed PUND plot is missing df")
    return _pund_stats_from_analysis_df(analyzed.df, device_params)
end

function compute_pund_voltage_stats(
    filepath::AbstractString,
    measurement_params::Dict{Symbol,Any},
)::Dict{Symbol,Any}
    fname = basename(filepath)
    dir = dirname(filepath)
    fatigue_count = measurement_params[:fatigue_count]
    wakeup_V = measurement_params[:wakeup_V]

    if fatigue_count > 0
        fatigue_df = _load_ruo2_pund_fatigue_file(filepath)
        return pund_voltage_stats(_select_pund_fatigue_cycle(fatigue_df, Int(fatigue_count)))
    elseif isfinite(Float64(wakeup_V))
        df = read_pund_wakeup_amplitude(fname, dir, Float64(wakeup_V), 1)
        return pund_voltage_stats(df)
    end

    return pund_voltage_stats(read_fe_pund(fname, dir))
end

function compute_pund_stats(
    filepath::AbstractString,
    measurement_params::Dict{Symbol,Any},
    device_params::Dict{Symbol,Any},
)::Dict{Symbol,Any}
    fname = basename(filepath)
    dir = dirname(filepath)
    fatigue_count = measurement_params[:fatigue_count]
    wakeup_V = measurement_params[:wakeup_V]

    if fatigue_count > 0
        fatigue_df = _load_ruo2_pund_fatigue_file(filepath)
        df = _select_pund_fatigue_cycle(fatigue_df, Int(fatigue_count))
        return _pund_stats_from_analysis_df(df, device_params)
    elseif isfinite(Float64(wakeup_V))
        df = read_pund_wakeup_amplitude(fname, dir, Float64(wakeup_V), 1)
        analyzed = analyze_pund_and_pn(df).pund
        return _pund_stats_from_analysis_df(analyzed, device_params)
    end

    return _pund_stats_from_analysis_df(read_fe_pund(fname, dir), device_params)
end
