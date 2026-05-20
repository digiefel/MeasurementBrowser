using DataLoader: read_pund_file, read_pund_wakeup_file
using DataAnalysis: detect_pund_pulses, analyze_pund
using DataFrames: nrow
using Statistics: mean

"""Compute PUND history and waveform stats after all measurements for a device are known."""
function compute_and_add_measurement_stats!(
    ::RuO2Project,
    measurements::Vector{MeasurementInfo},
    files::Vector{SourceFile},
)
    by_device = Dict{String,Vector{MeasurementInfo}}()
    headers = Dict(file.filepath => file.header_summary for file in files)
    fatigue_files = Dict{String,Any}()

    for measurement in measurements
        device_key = device_path_key(measurement.device_info)
        push!(get!(by_device, device_key, MeasurementInfo[]), measurement)
    end

    for device_measurements in values(by_device)
        sort!(device_measurements, by=measurement_timestamp_key)
        wakeup_counts = Dict{Tuple{Float64,Float64},Float64}()
        seen_wakeup_events = Set{Tuple{String,Float64,Float64}}()
        wakeup_count_so_far = 0.0
        wakeup_f_so_far = NaN
        wakeup_V_so_far = NaN
        fatigue_count_so_far = 0
        fatigue_f_so_far = NaN
        fatigue_V_so_far = NaN

        for measurement in device_measurements
            kind = measurement.measurement_kind
            (kind === :pund || kind === :pn || kind === :wakeup_pn || kind === :wakeup_pund) || continue

            params = measurement.parameters
            stats = measurement.stats
            empty!(stats)
            if kind === :wakeup_pn || kind === :wakeup_pund
                header = headers[measurement.filepath]
                wakeup_V = Float64(params[:wakeup_V])
                wakeup_f = parse(Float64, first(split(header["fatigue_freq"], ',')))
                local_count = parse(Float64, first(split(header["fatigue_count"], ',')))

                condition_key = (wakeup_V, wakeup_f)
                event_key = (measurement.filepath, wakeup_V, wakeup_f)
                if !(event_key in seen_wakeup_events)
                    wakeup_counts[condition_key] = get(wakeup_counts, condition_key, 0.0) + local_count
                    push!(seen_wakeup_events, event_key)
                end

                wakeup_count_so_far = wakeup_counts[condition_key]
                wakeup_f_so_far = wakeup_f
                wakeup_V_so_far = wakeup_V
                stats[:wakeup_count] = wakeup_count_so_far
                stats[:wakeup_f] = wakeup_f_so_far
                stats[:wakeup_V] = wakeup_V_so_far
                stats[:fatigue_count] = fatigue_count_so_far
                stats[:fatigue_f] = fatigue_f_so_far
                stats[:fatigue_V] = fatigue_V_so_far
            elseif kind === :pund
                source_is_fatigue = is_pund_fatigue_file(measurement.filepath)
                if source_is_fatigue
                    header = headers[measurement.filepath]
                    fatigue_count_so_far = Int(params[:fatigue_idx])
                    fatigue_f_so_far = parse(Float64, first(split(header["fatigue_freq"], ',')))
                    fatigue_V_so_far = parse(Float64, first(split(header["vmax"], ',')))
                end

                stats[:wakeup_count] = wakeup_count_so_far
                stats[:wakeup_f] = wakeup_f_so_far
                stats[:wakeup_V] = wakeup_V_so_far
                stats[:fatigue_count] = fatigue_count_so_far
                stats[:fatigue_f] = fatigue_f_so_far
                stats[:fatigue_V] = fatigue_V_so_far
            end

            merge!(stats, compute_pund_stats(
                measurement,
                measurement.device_info.parameters;
                fatigue_files=fatigue_files,
            ))
        end
    end
    return nothing
end

"""
Compute waveform stats for one already-expanded PUND-family logical measurement.
"""
function compute_pund_stats(
    measurement::MeasurementInfo,
    device_params::Dict{Symbol,Any}=Dict{Symbol,Any}(),
    ;
    fatigue_files::Dict{String,Any}=Dict{String,Any}(),
)::Dict{Symbol,Any}
    filepath = measurement.filepath
    fname = basename(filepath)
    dir = dirname(filepath)
    params = measurement.parameters

    if is_pund_fatigue_file(filepath)
        fatigue_count = Int(params[:fatigue_idx])
        fatigue_df = get!(fatigue_files, filepath) do
            _load_ruo2_pund_fatigue_file(filepath)
        end
        return pund_stats_from_waveform(
            _select_pund_fatigue_cycle(fatigue_df, fatigue_count),
            device_params,
        )
    elseif measurement.measurement_kind === :wakeup_pn || measurement.measurement_kind === :wakeup_pund
        df = _select_pund_wakeup_readout(filepath, Float64(params[:wakeup_V]))
        return pund_stats_from_waveform(df, device_params)
    end

    return pund_stats_from_waveform(read_pund_file(fname, dir), device_params)
end

function _select_pund_wakeup_readout(filepath::AbstractString, wakeup_V::Float64)
    wakeup_df = read_pund_wakeup_file(filepath)
    readout_df = wakeup_df[wakeup_df.wakeup_V .== wakeup_V, [:time, :voltage, :current]]
    nrow(readout_df) > 0 || error("No wakeup readout for wakeup_V=$wakeup_V in $filepath")
    return readout_df
end

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
