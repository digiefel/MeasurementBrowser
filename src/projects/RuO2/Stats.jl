using DataFrames: DataFrame, nrow
using DataAnalysis: extract_tlm_geometry_from_params
using Statistics: mean

"""Analyze one PUND-family measurement and add plotting-unit columns."""
function process_pund_measurement_data(
    measurement::MeasurementInfo,
    data::DataFrame,
)::DataFrame
    kind = measurement.measurement_kind
    segment = kind === :pn || kind === :wakeup_pn ? :pn :
        kind === :wakeup_pund ? :pund :
        nothing
    group_size = 5
    analyzed = if segment === :pn
        result = analyze_pund_and_pn(data)
        group_size = result.pn_group_size
        result.pn
    elseif segment === :pund
        result = analyze_pund_and_pn(data)
        group_size = result.pund_group_size
        result.pund
    else
        analyze_pund(data)
    end

    analyzed === nothing && return DataFrame()
    processed = DataFrame(analyzed)
    processed.time_us = processed.time .* 1e6
    processed.current_uA = processed.current .* 1e6
    processed.i_fe_uA = processed.I_FE .* 1e6
    processed.q_fe_pC = processed.Q_FE .* 1e12
    processed.pulse_group_size = fill(group_size, nrow(processed))

    finite_q = filter(isfinite, processed.Q_FE)
    q_mean = isempty(finite_q) ? 0.0 : mean(finite_q)
    q_centered = processed.Q_FE .- q_mean
    processed.q_centered_pC = q_centered .* 1e12

    params = merge(measurement.device_info.parameters, measurement.parameters)
    area_um2 = get(params, :area_um2, nothing)
    if area_um2 !== nothing
        area_cm2 = Float64(area_um2) / 1e8
        processed.p_fe_uC_cm2 = (q_centered ./ area_cm2) .* 1e6
    end
    return processed
end

"""Fit one current-voltage measurement and add resistance columns."""
function process_iv_measurement_data(
    measurement::MeasurementInfo,
    data::DataFrame,
)::DataFrame
    mask = isfinite.(data.i) .& isfinite.(data.v)
    i = Float64.(data.i[mask])
    v = Float64.(data.v[mask])
    n = length(i)
    sx = sum(i)
    sy = sum(v)
    sxx = sum(i .^ 2)
    sxy = sum(i .* v)
    denom = n * sxx - sx^2
    resistance_fit = NaN
    voltage_offset = NaN
    if abs(denom) > 1e-20
        resistance_fit = (n * sxy - sx * sy) / denom
        voltage_offset = (sy - resistance_fit * sx) / n
    end

    params = merge(measurement.device_info.parameters, measurement.parameters)
    length_um, width_um = extract_tlm_geometry_from_params(params)
    sheet_resistance = (
        isfinite(resistance_fit) &&
        isfinite(length_um) &&
        isfinite(width_um) &&
        length_um > 0
    ) ? resistance_fit * width_um / length_um : NaN
    resistance = v ./ i

    return DataFrame(
        i=i,
        v=v,
        resistance_ohm=resistance,
        valid_resistance=isfinite.(resistance) .& (abs.(resistance) .< 1e9),
        fit_resistance_ohm=fill(resistance_fit, n),
        fit_offset_v=fill(voltage_offset, n),
        sheet_resistance_ohm=fill(sheet_resistance, n),
    )
end

"""Load and process one logical RuO2 measurement through the workspace."""
function process_measurement_data(
    workspace::Workspace.Workspace{RuO2Project},
    measurement::MeasurementInfo,
)::DataFrame
    data = only(read_measurement_data(workspace, [measurement]))
    measurement.measurement_kind in (:pund, :pn, :wakeup_pn, :wakeup_pund) &&
        return process_pund_measurement_data(measurement, data)
    measurement.measurement_kind in (:iv, :breakdown, :tlm4p) &&
        return process_iv_measurement_data(measurement, data)
    return data
end

"""Compute PUND history and waveform stats after all measurements for a device are known."""
function compute_and_add_measurement_stats!(
    ::RuO2Project,
    measurements::Vector{MeasurementInfo},
    files::Vector{SourceFile},
    ;
    on_progress::Union{Nothing,Function}=nothing,
)::Vector{MeasurementAnalysisFailure}
    by_device = Dict{String,Vector{MeasurementInfo}}()
    fatigue_files = Dict{String,Any}()
    fatigue_cycle_rows = Dict{String,Dict{Int,Vector{Int}}}()
    fatigue_pulses = Dict{String,Vector{UnitRange{Int}}}()
    wakeup_files = Dict{String,Any}()
    wakeup_voltage_rows = Dict{String,Dict{Float64,Vector{Int}}}()
    wakeup_splits = Dict{Tuple{String,Float64},Any}()
    failures = MeasurementAnalysisFailure[]
    total = count(_ruo2_needs_pund_stats, measurements)
    processed = 0

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
            _ruo2_needs_pund_stats(measurement) || continue
            check_cancel()

            params = measurement.parameters
            stats = measurement.stats
            empty!(stats)
            waveform_stats = Dict{Symbol,Any}()
            try
                waveform_stats = compute_pund_stats(
                    measurement,
                    measurement.device_info.parameters;
                    fatigue_files=fatigue_files,
                    fatigue_cycle_rows=fatigue_cycle_rows,
                    fatigue_pulses=fatigue_pulses,
                    wakeup_files=wakeup_files,
                    wakeup_voltage_rows=wakeup_voltage_rows,
                    wakeup_splits=wakeup_splits,
                )
            catch err
                bt = catch_backtrace()
                @error(
                    "Measurement analysis failed",
                    project=project_name(RUO2_PROJECT),
                    file=measurement.filepath,
                    measurement=measurement.unique_id,
                    exception=(err, bt),
                )
                push!(failures, MeasurementAnalysisFailure(
                    measurement.filepath,
                    measurement.unique_id,
                    sprint(showerror, err),
                ))
            end

            if kind === :wakeup_pn || kind === :wakeup_pund
                wakeup_V = Float64(params[:wakeup_V])
                wakeup_f = Float64(params[:wakeup_f])
                local_count = Float64(params[:wakeup_count])

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
                stats[:fatigue_count] = fatigue_count_so_far
                isfinite(wakeup_f_so_far) && (stats[:wakeup_f] = wakeup_f_so_far)
                isfinite(wakeup_V_so_far) && (stats[:wakeup_V] = wakeup_V_so_far)
                isfinite(fatigue_f_so_far) && (stats[:fatigue_f] = fatigue_f_so_far)
                isfinite(fatigue_V_so_far) && (stats[:fatigue_V] = fatigue_V_so_far)
            elseif kind === :pund
                source_is_fatigue = is_pund_fatigue_file(measurement.filepath)
                if source_is_fatigue
                    fatigue_count_so_far = Int(params[:fatigue_idx])
                    # Current RuO2 fatigue files use the same frequency/amplitude
                    # for fatigue pulses and the selected readout waveform.
                    haskey(waveform_stats, :frequency_kHz) &&
                        (fatigue_f_so_far = waveform_stats[:frequency_kHz] * 1000)
                    haskey(waveform_stats, :V_amp) &&
                        (fatigue_V_so_far = waveform_stats[:V_amp])
                end

                stats[:wakeup_count] = wakeup_count_so_far
                stats[:fatigue_count] = fatigue_count_so_far
                isfinite(wakeup_f_so_far) && (stats[:wakeup_f] = wakeup_f_so_far)
                isfinite(wakeup_V_so_far) && (stats[:wakeup_V] = wakeup_V_so_far)
                isfinite(fatigue_f_so_far) && (stats[:fatigue_f] = fatigue_f_so_far)
                isfinite(fatigue_V_so_far) && (stats[:fatigue_V] = fatigue_V_so_far)
            end

            merge!(stats, waveform_stats)
            processed += 1
            on_progress !== nothing && on_progress((
                total=total,
                processed=processed,
                current_path=measurement.filepath,
            ))
        end
    end
    return failures
end

function _ruo2_needs_pund_stats(measurement::MeasurementInfo)::Bool
    kind = measurement.measurement_kind
    return kind === :pund || kind === :pn || kind === :wakeup_pn || kind === :wakeup_pund
end

"""
Compute waveform stats for one already-expanded PUND-family logical measurement.
"""
function compute_pund_stats(
    measurement::MeasurementInfo,
    device_params::Dict{Symbol,Any}=Dict{Symbol,Any}(),
    ;
    fatigue_files::Dict{String,Any}=Dict{String,Any}(),
    fatigue_cycle_rows::Dict{String,Dict{Int,Vector{Int}}}=Dict{String,Dict{Int,Vector{Int}}}(),
    fatigue_pulses::Dict{String,Vector{UnitRange{Int}}}=Dict{String,Vector{UnitRange{Int}}}(),
    wakeup_files::Dict{String,Any}=Dict{String,Any}(),
    wakeup_voltage_rows::Dict{String,Dict{Float64,Vector{Int}}}=Dict{String,Dict{Float64,Vector{Int}}}(),
    wakeup_splits::Dict{Tuple{String,Float64},Any}=Dict{Tuple{String,Float64},Any}(),
)::Dict{Symbol,Any}
    filepath = measurement.filepath
    fname = basename(filepath)
    dir = dirname(filepath)
    params = measurement.parameters

    if is_pund_fatigue_file(filepath)
        fatigue_count = Int(params[:fatigue_idx])
        fatigue_df = get!(fatigue_files, filepath) do
            read_pund_fatigue_file(filepath)
        end
        rows_by_cycle = get!(fatigue_cycle_rows, filepath) do
            _rows_by_value(Int, fatigue_df.cycle)
        end
        df = _select_pund_fatigue_cycle(fatigue_df, fatigue_count, rows_by_cycle)
        pulses = get!(fatigue_pulses, filepath) do
            detect_pund_pulses(df.time, df.voltage, df.current).pulses
        end
        return pund_stats_from_waveform(df, device_params; pulses)
    elseif measurement.measurement_kind === :wakeup_pn || measurement.measurement_kind === :wakeup_pund
        wakeup_df = get!(wakeup_files, filepath) do
            read_pund_wakeup_file(filepath)
        end
        rows_by_voltage = get!(wakeup_voltage_rows, filepath) do
            _rows_by_value(Float64, wakeup_df.wakeup_V)
        end
        wakeup_V = Float64(params[:wakeup_V])
        split = get!(wakeup_splits, (filepath, wakeup_V)) do
            df = _select_pund_wakeup_readout(wakeup_df, wakeup_V, rows_by_voltage)
            analyze_pund_and_pn(df)
        end
        selected = measurement.measurement_kind === :wakeup_pn ? split.pn : split.pund
        return pund_stats_from_waveform(selected, device_params)
    end

    return pund_stats_from_waveform(read_pund_file(fname, dir), device_params)
end

function _rows_by_value(::Type{T}, values)::Dict{T,Vector{Int}} where {T}
    rows = Dict{T,Vector{Int}}()
    for (index, value) in pairs(values)
        push!(get!(rows, T(value), Int[]), Int(index))
    end
    return rows
end

function _select_pund_fatigue_cycle(
    fatigue_df,
    cycle::Integer,
    rows_by_cycle::Dict{Int,Vector{Int}},
)
    rows = get(rows_by_cycle, Int(cycle), Int[])
    isempty(rows) && error("PUND fatigue data has no rows for cycle $cycle")
    return DataFrame(
        time=fatigue_df.time[rows],
        current=fatigue_df.current[rows],
        voltage=fatigue_df.voltage[rows],
    )
end

function _select_pund_wakeup_readout(filepath::AbstractString, wakeup_V::Float64)
    wakeup_df = read_pund_wakeup_file(filepath)
    rows_by_voltage = _rows_by_value(Float64, wakeup_df.wakeup_V)
    return _select_pund_wakeup_readout(wakeup_df, wakeup_V, rows_by_voltage)
end

function _select_pund_wakeup_readout(
    wakeup_df,
    wakeup_V::Float64,
    rows_by_voltage::Dict{Float64,Vector{Int}},
)
    rows = get(rows_by_voltage, wakeup_V, Int[])
    isempty(rows) && error("No wakeup readout for wakeup_V=$wakeup_V")
    return DataFrame(
        time=wakeup_df.time[rows],
        voltage=wakeup_df.voltage[rows],
        current=wakeup_df.current[rows],
    )
end

"""
Stats attached to one PUND logical measurement. These describe the selected waveform,
not the voltage/frequency settings requested in the file header.
"""
function pund_stats_from_waveform(
    df,
    device_params::Dict{Symbol,Any}=Dict{Symbol,Any}();
    pulses::Union{Nothing,Vector{UnitRange{Int}}}=nothing,
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

    waveform_pulses = hasproperty(df, :pulse_idx) ?
        _pulse_ranges_from_indices(df.pulse_idx) :
        pulses === nothing ? detect_pund_pulses(df.time, df.voltage, df.current).pulses : pulses
    if !isempty(waveform_pulses)
        t_start = df.time[first(waveform_pulses[1])]
        t_end = df.time[last(waveform_pulses[end])]
        total_span = t_end - t_start
        if total_span > 0
            pulse_period = total_span / length(waveform_pulses)
            stats[:frequency_kHz] = round_one(1.0 / (pulse_period * 1e3))
        end
    end

    if haskey(device_params, :area_um2)
        area_um2 = Float64(device_params[:area_um2])
        if isfinite(area_um2) && area_um2 > 0
            df_an = hasproperty(df, :Q_FE) && hasproperty(df, :pulse_idx) ?
                df : analyze_pund(df; pulses=waveform_pulses)
            stats[:Pr_uCcm2] = round(pund_pr_value(df_an, area_um2); digits=3)
        end
    end

    return stats
end

function _pulse_ranges_from_indices(pulse_idx)
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

"""Return mean Pr over the valid loop(s) represented by this logical measurement."""
function pund_pr_value(df, area_um2::Float64)
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
    isempty(pr_vals) && error("Expected at least one PUND loop")
    return mean(pr_vals)
end
