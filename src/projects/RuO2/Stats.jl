using DataLoader: read_fe_pund, read_pund_fatigue_cycle, read_pund_wakeup_amplitude, read_pund_wakeup_reps
using DataAnalysis: detect_pund_pulses, analyze_pund, analyze_pund_and_pn
using DataFrames: nrow
using Statistics: mean

"""
    compute_pund_stats(filepath, measurement_params, device_params) -> Dict{Symbol,Any}
    compute_pund_stats_from_analyzed_plot(analyzed, device_params) -> Dict{Symbol,Any}

Compute data-derived statistics for a PUND measurement by reading the CSV.
Returns a Dict with keys:
  :voltage_max_V      — peak positive voltage relative to baseline
  :voltage_baseline_V — DC voltage level (mean of first 9 points)
  :voltage_min_V      — peak negative voltage relative to baseline
  :frequency_kHz      — inverse period of a single triangular pulse
  :Pr_max_uCcm2       — maximum remnant polarization across repetitions (µC/cm², if area_um2 available)
Returns an empty Dict on failure (missing file, unreadable data, etc.).
"""
function _pund_stats_from_analysis_df(df, device_params::Dict{Symbol,Any})::Dict{Symbol,Any}
    result = Dict{Symbol,Any}()
    nrow(df) == 0 && return result

    t = df.time
    V = df.voltage
    I = df.current

    n_pre = min(9, length(V))
    V_baseline = mean(V[1:n_pre])

    result[:voltage_baseline_V] = round(V_baseline; digits=1)
    result[:voltage_max_V] = round(maximum(V); digits=1)
    result[:voltage_min_V] = round(minimum(V); digits=1)

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
    try
        hasproperty(analyzed, :df) || return Dict{Symbol,Any}()
        return _pund_stats_from_analysis_df(analyzed.df, device_params)
    catch
        return Dict{Symbol,Any}()
    end
end

function compute_pund_stats(filepath::AbstractString,
    measurement_params::Dict{Symbol,Any},
    device_params::Dict{Symbol,Any})::Dict{Symbol,Any}

    isfile(filepath) || return Dict{Symbol,Any}()
    result = Dict{Symbol,Any}()

    try
        fname = basename(filepath)
        dir = dirname(filepath)

        if haskey(measurement_params, :fatigue_cycle)
            cycle = Int(measurement_params[:fatigue_cycle])
            df = read_pund_fatigue_cycle(fname, dir, cycle)
            nrow(df) == 0 && return result
            merge!(result, _pund_stats_from_analysis_df(df, device_params))
        elseif haskey(measurement_params, :amplitude_V)
            amp = Float64(measurement_params[:amplitude_V])
            rep = Int(get(measurement_params, :wakeup_rep, 1))
            df = read_pund_wakeup_amplitude(fname, dir, amp, rep)
            nrow(df) == 0 && return result
            analyzed = analyze_pund_and_pn(df).pund
            merge!(result, _pund_stats_from_analysis_df(analyzed, device_params))
        else
            df = read_fe_pund(fname, dir)
            nrow(df) == 0 && return result
            merge!(result, _pund_stats_from_analysis_df(df, device_params))
        end
    catch
        # Return whatever we managed to compute (or empty on total failure)
    end

    return result
end
