"""
    find_zero_crossing(x, y) -> Union{Float64, Nothing}

Find the x-value where y crosses zero using linear interpolation.
Returns the first zero crossing found, or nothing if no crossing exists.
"""
function find_zero_crossing(x::AbstractVector, y::AbstractVector)
    for i in 1:(length(y)-1)
        y1, y2 = y[i], y[i+1]
        (isnan(y1) || isnan(y2)) && continue
        if y1 * y2 < 0  # sign change
            # Linear interpolation: y = y1 + (y2-y1)/(x2-x1) * (x - x1) = 0
            # x = x1 - y1 * (x2-x1) / (y2-y1)
            x_cross = x[i] - y1 * (x[i+1] - x[i]) / (y2 - y1)
            return x_cross
        elseif y1 == 0.0
            return x[i]
        end
    end
    return nothing
end

"""
    fit_rc_transient(t, I, dV_dt; min_points=5) -> NamedTuple

Fit an RC transient response to current data during a voltage ramp.

Model: I(t) = I_final + (I_initial - I_final) * exp(-(t-t0)/τ)

Fitting strategy:
- τ (time constant): fitted from first 20% of data (steepest transient)
- C (capacitance): calculated from last 20% of data (steady state: I_ss = C * dV/dt)

Returns NamedTuple with:
- C: capacitance (Farads)
- tau: RC time constant (seconds)
- R: series resistance (Ohms)
- I_initial: current at start of edge (Amps)
- I_final: steady-state current (Amps)
- fit_quality: R² of the linearized fit (0-1, higher is better)

Returns nothing if fit fails (insufficient data, bad fit, etc.)
"""
function fit_rc_transient(t::AbstractVector, I::AbstractVector, dV_dt::Float64;
    min_points::Int=5, tau_min::Float64=1e-9, tau_max::Float64=1e-3)

    n = length(t)
    n < min_points && return nothing
    abs(dV_dt) < 1e-12 && return nothing  # No voltage change

    # I_initial: current at the start (first few points)
    I_initial = median(filter(isfinite, I[1:min(3, n)]))
    !isfinite(I_initial) && return nothing

    # I_final (steady state): estimate from last 20% of data
    start_idx = max(1, 4n ÷ 5)
    I_latter = @view I[start_idx:end]
    I_final = median(filter(isfinite, I_latter))
    !isfinite(I_final) && return nothing

    # Need a significant difference between initial and final
    delta_I = I_final - I_initial
    abs(delta_I) < 1e-12 && return nothing

    # Capacitance from steady state (last 20%): C = |I_final / dV_dt|
    C = abs(I_final / dV_dt)
    !isfinite(C) && return nothing
    C < 1e-15 && return nothing  # Unreasonably small capacitance

    # Fit τ using only the first 20% of data (steepest transient)
    # Model: I(t) = I_final + (I_initial - I_final) * exp(-(t-t0)/τ)
    # Rearranging: (I - I_final) / (I_initial - I_final) = exp(-(t-t0)/τ)
    # Let y = (I - I_final) / (I_initial - I_final), then ln(y) = -(t-t0)/τ

    t0 = t[1]
    end_idx_fit = max(min_points, n ÷ 5)  # First 20%

    y_vals = Float64[]
    t_vals = Float64[]

    for i in 1:end_idx_fit
        Ii = I[i]
        !isfinite(Ii) && continue

        # y = (I - I_final) / (I_initial - I_final)
        y = (Ii - I_final) / (I_initial - I_final)

        # For valid exponential decay, y should be positive and ≤ 1
        # y=1 at t=t0, y→0 as t→∞
        0.01 < y < 1.5 || continue

        push!(y_vals, y)
        push!(t_vals, t[i] - t0)
    end

    length(y_vals) < min_points && return nothing

    # Linear fit of ln(y) vs t: ln(y) = -t/τ  =>  slope = -1/τ
    ln_y = log.(y_vals)

    # Simple linear regression: ln(y) = a + b*t
    n_fit = length(t_vals)
    t_mean = mean(t_vals)
    ln_y_mean = mean(ln_y)

    num = sum((t_vals .- t_mean) .* (ln_y .- ln_y_mean))
    den = sum((t_vals .- t_mean) .^ 2)

    abs(den) < 1e-30 && return nothing

    slope = num / den
    # slope = -1/τ, so τ = -1/slope
    slope >= 0 && return nothing  # Should be negative for decay

    tau = -1.0 / slope

    # Sanity check on τ
    (tau < tau_min || tau > tau_max) && return nothing

    # Calculate R from τ = RC
    R = tau / C
    !isfinite(R) && return nothing

    # Compute fit quality (R²) on the first 20%
    intercept = ln_y_mean - slope * t_mean
    ln_y_pred = intercept .+ slope .* t_vals
    ss_res = sum((ln_y .- ln_y_pred) .^ 2)
    ss_tot = sum((ln_y .- ln_y_mean) .^ 2)
    r_squared = ss_tot > 1e-30 ? 1.0 - ss_res / ss_tot : 0.0

    # Require reasonable fit quality
    r_squared < 0.3 && return nothing

    return (C=C, tau=tau, R=R, I_initial=I_initial, I_final=I_final, fit_quality=r_squared)
end

"""
    analyze_ud_pulses(t, V, I, pulses, n_groups) -> NamedTuple

Analyze U (up) and D (down) pulses from a PUND sequence to extract
capacitance and RC time constants.

For each U/D pulse:
1. Take the first half (rising from 0V to peak)
2. Fit τ from first 20% (steepest transient), C from last 20% (steady state)
3. Model: I(t) = I_final + (I_initial - I_final) * exp(-t/τ)

Returns NamedTuple with per-group results:
- groups: Vector of NamedTuples with (group_idx, C_U, tau_U, R_U, fit_quality_U, C_D, ...)
- C_avg: Median capacitance across all successful fits
- tau_avg: Median time constant
- R_avg: Median resistance
"""
function analyze_ud_pulses(t::AbstractVector, V::AbstractVector, I::AbstractVector,
    pulses::Vector{UnitRange{Int}}, n_groups::Int; has_initial_poling::Bool=true)

    group_results = NamedTuple{
        (:group_idx, :C_U, :tau_U, :R_U, :fit_quality_U, :C_D, :tau_D, :R_D, :fit_quality_D),
        Tuple{Int,
              Union{Float64,Nothing}, Union{Float64,Nothing}, Union{Float64,Nothing}, Union{Float64,Nothing},
              Union{Float64,Nothing}, Union{Float64,Nothing}, Union{Float64,Nothing}, Union{Float64,Nothing}}
    }[]

    all_C = Float64[]
    all_tau = Float64[]
    all_R = Float64[]

    for group_idx in 1:n_groups
        group_size = has_initial_poling ? 5 : 4
        U_offset = has_initial_poling ? 3 : 2
        D_offset = has_initial_poling ? 5 : 4
        U_pulse_num = (group_idx - 1) * group_size + U_offset
        D_pulse_num = (group_idx - 1) * group_size + D_offset

        U_pulse_num > length(pulses) && continue
        D_pulse_num > length(pulses) && continue

        U_range = pulses[U_pulse_num]
        D_range = pulses[D_pulse_num]

        C_U, tau_U, R_U, fit_quality_U = nothing, nothing, nothing, nothing
        C_D, tau_D, R_D, fit_quality_D = nothing, nothing, nothing, nothing

        # Analyze U pulse - entire rising edge (first half of triangular pulse)
        # fit_rc_transient uses first 20% for τ, last 20% for C
        if length(U_range) >= 10
            t_U = t[U_range]
            V_U = V[U_range]
            I_U = I[U_range]

            # Find apex (max |V|) - rising edge is 1:apex_idx
            apex_idx = argmax(abs.(V_U))

            if apex_idx > 10
                t_rise = t_U[1:apex_idx]
                V_rise = V_U[1:apex_idx]
                I_rise = I_U[1:apex_idx]

                dV_dt = (V_rise[end] - V_rise[1]) / (t_rise[end] - t_rise[1])
                if abs(dV_dt) > 1e-6
                    result = fit_rc_transient(t_rise, I_rise, dV_dt)
                    if result !== nothing
                        C_U, tau_U, R_U = result.C, result.tau, result.R
                        fit_quality_U = result.fit_quality
                        push!(all_C, result.C)
                        push!(all_tau, result.tau)
                        push!(all_R, result.R)
                    end
                end
            end
        end

        # Analyze D pulse - entire rising edge (going negative)
        if length(D_range) >= 10
            t_D = t[D_range]
            V_D = V[D_range]
            I_D = I[D_range]

            apex_idx = argmax(abs.(V_D))

            if apex_idx > 10
                t_rise = t_D[1:apex_idx]
                V_rise = V_D[1:apex_idx]
                I_rise = I_D[1:apex_idx]

                dV_dt = (V_rise[end] - V_rise[1]) / (t_rise[end] - t_rise[1])
                if abs(dV_dt) > 1e-6
                    result = fit_rc_transient(t_rise, I_rise, dV_dt)
                    if result !== nothing
                        C_D, tau_D, R_D = result.C, result.tau, result.R
                        fit_quality_D = result.fit_quality
                        push!(all_C, result.C)
                        push!(all_tau, result.tau)
                        push!(all_R, result.R)
                    end
                end
            end
        end

        push!(group_results, (
            group_idx=group_idx,
            C_U=C_U, tau_U=tau_U, R_U=R_U, fit_quality_U=fit_quality_U,
            C_D=C_D, tau_D=tau_D, R_D=R_D, fit_quality_D=fit_quality_D
        ))
    end

    C_avg = isempty(all_C) ? nothing : median(all_C)
    tau_avg = isempty(all_tau) ? nothing : median(all_tau)
    R_avg = isempty(all_R) ? nothing : median(all_R)

    return (groups=group_results, C_avg=C_avg, tau_avg=tau_avg, R_avg=R_avg)
end

"""
    detect_pund_pulses(t, V, I; smooth_window=9, dV_thresh_mult=5.0,
        expand_win=5, min_pulse_len=20, min_V_amp_mult=5.0,
        trough_thresh=0.05) -> NamedTuple

Detect triangular voltage pulses in a PUND waveform and check polarity
consistency. Returns a diagnostic NamedTuple with all intermediate results;
never errors on inconsistent polarity.

Detection has two stages (fallback chain):
  1. Derivative-based: finds pulses via smoothed dV/dt threshold + expansion.
  2. Trough-based (fallback): splits at near-zero |V| regions.
"""
function detect_pund_pulses(t, V, I;
    smooth_window::Int=9,
    dV_thresh_mult::Float64=5.0,
    expand_win::Int=5,
    min_pulse_len::Int=20,
    min_V_amp_mult::Float64=5.0,
    trough_thresh::Float64=0.05)

    N = length(t)

    # Baseline current offset
    n0 = min(10, N)
    baseline_I = mean(skipmissing(filter(!isnan, I[1:n0])))
    I_shifted = I .- baseline_I

    # Smoothed derivative and threshold (1/2 of max |dV|)
    dV = smoothdata([0.0; diff(V)], :movmedian, smooth_window)
    dV_max = maximum(abs.(dV))
    dV_threshold = dV_max / 2

    # Voltage baseline: mean of first few points (may be non-zero, e.g. 1 V DC offset)
    n_pre = min(9, length(V))
    V_baseline = mean(V[1:n_pre])
    Vc = V .- V_baseline
    V_noise = mean(abs.(Vc[1:n_pre]))
    min_V_amp = V_noise * min_V_amp_mult

    # Stage 1: derivative-based detection
    pulses = UnitRange{Int}[]
    detection_method = :none
    pulse_mask = abs.(dV) .> dV_threshold
    expanded_mask = copy(pulse_mask)

    # TODO disable this for now. Probably not needed and overengineered.
    # if dV_max > 1e-12
    #     detection_method = :derivative
    #     for i in (expand_win+1):(length(pulse_mask)-expand_win)
    #         if any(pulse_mask[(i-expand_win):(i+expand_win)])
    #             expanded_mask[i] = true
    #         end
    #     end
    #     all_pulses = _true_runs(expanded_mask)

    #     for pulse in all_pulses
    #         if length(pulse) >= min_pulse_len
    #             pulse_V_range = maximum(abs.(Vc[pulse]))
    #             if pulse_V_range >= min_V_amp
    #                 push!(pulses, pulse)
    #             end
    #         end
    #     end
    # end

    # Stage 2: trough-based fallback (only if stage 1 found nothing)
    if isempty(pulses)
        detection_method = :trough_fallback
        absVc = abs.(Vc)
        V_peak = maximum(absVc)
        trough_threshold = V_peak * trough_thresh
        in_trough = absVc .< trough_threshold
        i = 1
        while i <= N
            while i <= N && in_trough[i]
                i += 1
            end
            i > N && break
            pstart = i
            while i <= N && !in_trough[i]
                i += 1
            end
            pend = i > N ? N : i - 1
            seg = pstart:pend
            if length(seg) >= min_pulse_len
                push!(pulses, seg)
            end
        end
    end

    # Per-pulse polarity check
    polarity_info = NamedTuple{(:pulse_idx, :V_avg, :V_sign),
        Tuple{Int,Float64,Int}}[]
    total = 0
    for (pi, pulse) in enumerate(pulses)
        half = pulse[1:end÷2]
        V_avg = mean(Vc[half])
        if abs(V_avg) > 0
            total += 1
            sv = Int(sign(V_avg))
            push!(polarity_info, (pulse_idx=pi, V_avg=V_avg, V_sign=sv))
        else
            push!(polarity_info, (pulse_idx=pi, V_avg=V_avg, V_sign=0))
        end
    end

    n_groups = length(pulses) ÷ 5
    remainder = length(pulses) % 5

    return (;
        t, V, I=I_shifted, dV, dV_threshold, dV_max,
        pulse_mask, expanded_mask,
        pulses, detection_method,
        polarity_info, total, n_groups, remainder,
        min_V_amp, baseline_V=V_baseline,
    )
end

"Contiguous runs of `true` in a BitVector."
function _true_runs(mask)
    runs = UnitRange{Int}[]
    start = nothing
    for i in eachindex(mask)
        if mask[i]
            start === nothing && (start = i)
        elseif start !== nothing
            push!(runs, start:i-1)
            start = nothing
        end
    end
    start !== nothing && push!(runs, start:lastindex(mask))
    return runs
end

function _pund_group_size(has_initial_poling::Bool)
    return has_initial_poling ? 5 : 4
end

function _pund_switch_codes(group_idx::Int, has_initial_poling::Bool)
    group_size = _pund_group_size(has_initial_poling)
    base = (group_idx - 1) * group_size
    if has_initial_poling
        return (P=base + 2, U=base + 3, N=base + 4, D=base + 5)
    end
    return (P=base + 1, U=base + 2, N=base + 3, D=base + 4)
end

function _center_current(I, N::Int)
    n0 = min(10, N)
    baseline_I = mean(skipmissing(filter(!isnan, I[1:n0])))
    return I .- baseline_I, baseline_I
end

function _slice_signal_dataframe(df::DataFrame, pulses::Vector{UnitRange{Int}}; from::Union{Nothing,Int}=nothing)
    isempty(pulses) && return DataFrame(df), UnitRange{Int}[]
    start_idx = from === nothing ? first(first(pulses)) : from
    stop_idx = last(last(pulses))
    row_range = start_idx:stop_idx
    sliced = DataFrame(df[row_range, :])
    offset = start_idx - 1
    shifted = [first(pulse)-offset:last(pulse)-offset for pulse in pulses]
    if hasproperty(sliced, :time) && !isempty(sliced.time)
        t0 = sliced.time[1]
        sliced.time = sliced.time .- t0
        hasproperty(sliced, :current_time) && (sliced.current_time = sliced.current_time .- t0)
        hasproperty(sliced, :voltage_time) && (sliced.voltage_time = sliced.voltage_time .- t0)
    end
    return sliced, shifted
end

function _pulse_voltage_sign(Vc::AbstractVector, pulse::UnitRange{Int})
    values = @view Vc[pulse]
    isempty(values) && return 0
    return sign(values[argmax(abs.(values))])
end

"""
Perform PUND analysis on the DataFrame
 and !DEBUG
    # 1. identify triangular voltage pulse train
    # 2. get current in P_n and U_n pulses
    # 3. subtract I_Pn - I_Un = I_FEn
    # 4. integrate I_FEn(t) dt = Q_FEn
    # 5. return analysis dataframe with Q_FE in each row

from this, Q_FE(V) or Q_FE(t) can be very easily extracted
"""
function analyze_pund(
    df::DataFrame;
    pulses::Union{Nothing,Vector{UnitRange{Int}}}=nothing,
    has_initial_poling::Bool=true,
    DEBUG::Bool=false,
)
    @assert all(["time", "voltage", "current"] .∈ Ref(names(df))) "columns :time, :voltage, :current must exist"
    t, V, I = df[!, :time], df[!, :voltage], df[!, :current]
    N = length(t)

    if pulses === nothing
        # ---- full pulse detection & polarity check ---------------------------------
        det = detect_pund_pulses(t, V, I)
        pulses = det.pulses
        I = det.I  # baseline-shifted
        I = -I # invert current to match expected polarity (positive for positive voltage pulses)

        if DEBUG
            n0 = min(10, N)
            baseline_I = mean(skipmissing(filter(!isnan, df[!, :current][1:n0])))
            @info "analyze_pund: baseline current offset" n0 baseline = baseline_I
            # @info "analyze_pund: detection" method = det.detection_method n_pulses = length(pulses) verdict = det.verdict
        end

        if !DEBUG
            @assert !isempty(pulses) "no valid pulses found; adjust derivative threshold or filtering parameters"
        end
    else
        # ---- pulse boundaries provided, skip detection -----------------------------
        I, baseline_I = _center_current(I, N)

        # Quick polarity check on first pulse (use centered voltage)
        n_vbl2 = min(9, N)
        V_baseline2 = mean(V[1:n_vbl2])
        Vc2 = V .- V_baseline2
        half = pulses[1][1:end÷2]
        V_avg = mean(Vc2[half])
        I_avg = mean(I[half])
        if abs(V_avg) > 0 && abs(I_avg) > 0 && sign(V_avg) != sign(I_avg)
            I = -I # commented out to disable polarity correction
        end
    end

    # ---- voltage baseline for sign checks (may be non-zero, e.g. 1 V DC offset) --------
    n_vbl = min(9, N)
    V_baseline = mean(V[1:n_vbl])
    Vc = V .- V_baseline

    # ---- consistency check and grouping ------------------------------------------------
    group_size = _pund_group_size(has_initial_poling)
    grouped_pulses = Vector{UnitRange{Int}}[]
    for i in 1:group_size:length(pulses)-group_size+1
        group = pulses[i:i+group_size-1]
        push!(grouped_pulses, group)
    end

    for group in grouped_pulses
        if has_initial_poling
            poling, P, U, Np, D = Tuple(group)
            sP, sPol = _pulse_voltage_sign(Vc, P), _pulse_voltage_sign(Vc, poling)
            # if !DEBUG
            #     @assert sPol == -sP &&
            #             _pulse_voltage_sign(Vc, U) == sP &&
            #             _pulse_voltage_sign(Vc, Np) == -sP &&
            #             _pulse_voltage_sign(Vc, D) == -sP "unexpected pulse ordering"
            # end
        else
            P, U, Np, D = Tuple(group)
            sP = _pulse_voltage_sign(Vc, P)
            # if !DEBUG
            #     @assert _pulse_voltage_sign(Vc, U) == sP &&
            #             _pulse_voltage_sign(Vc, Np) == -sP &&
            #             _pulse_voltage_sign(Vc, D) == -sP "unexpected pulse ordering"
            # end
        end
    end

    # ---- allocate result columns -------------------------------------------------------
    polarity = zeros(Int8, N)
    pulse_idx = zeros(Int, N)
    I_FE = zeros(eltype(I), N)
    Q_FE = fill!(similar(I), NaN)

    # ---- sample-wise subtraction helper -----------------------------------------------
    subtract_aligned(A, B) =
        length(A) == length(B) ?
        A .- B :
        A .- @view(B[round.(Int, LinRange(1, length(B), length(A)))])

    # ---- process every PUND group ------------------------------------------------------
    for (group_idx, group) in enumerate(grouped_pulses)
        if has_initial_poling
            poling, P, U, Np, D = Tuple(group)
        else
            P, U, Np, D = Tuple(group)
        end
        codes = _pund_switch_codes(group_idx, has_initial_poling)

        # polarity flags
        polarity[P] .= 1
        polarity[U] .= 1
        polarity[Np] .= -1
        polarity[D] .= -1
        if has_initial_poling
            pulse_idx[poling] .= (group_idx - 1) * group_size + 1
        end
        pulse_idx[P] .= codes.P
        pulse_idx[U] .= codes.U
        pulse_idx[Np] .= codes.N
        pulse_idx[D] .= codes.D

        # positive switching
        I_FE_P = subtract_aligned(I[P], I[U])
        I_FE[P] = I_FE_P

        # negative switching
        I_FE_N = subtract_aligned(I[Np], I[D])
        I_FE[Np] = I_FE_N
    end

    # ---- cumulative integral: trapezoidal rule ----------------------------------------
    dt = [zero(t[1]); diff(t)]
    dQ = I_FE .* dt
    q = 0.0
    switching_codes = Set{Int}()
    for group_idx in eachindex(grouped_pulses)
        codes = _pund_switch_codes(group_idx, has_initial_poling)
        push!(switching_codes, codes.P)
        push!(switching_codes, codes.N)
    end
    for i ∈ eachindex(I_FE)
        # Only integrate during P and N switching pulses
        if pulse_idx[i] in switching_codes
            q += dQ[i]
            Q_FE[i] = q
        else
            Q_FE[i] = NaN  # Set NaN outside P and N pulses
        end
    end

    # ---- assemble and return -----------------------------------------------------------
    df[!, :current] .= I # fix polarity
    return hcat(df, DataFrame(polarity=polarity, pulse_idx=pulse_idx, I_FE=I_FE, Q_FE=Q_FE))
end

"""
Analyze a PN segment by integrating the measured current directly over the cycle.
"""
function analyze_pn(df::DataFrame; DEBUG::Bool=false)
    @assert all(["time", "voltage", "current"] .∈ Ref(names(df))) "columns :time, :voltage, :current must exist"
    t, V, I_raw = df[!, :time], df[!, :voltage], df[!, :current]
    N = length(t)

    # Use the quiescent base section before the pulses start for baseline subtraction.
    # Threshold at 0.5% of peak to stop the baseline before the ramp begins.
    Vmax = isempty(V) ? 0.0 : maximum(abs.(V))
    first_pulse = Vmax > 0 ? findfirst(v -> abs(v) > 0.005 * Vmax, V) : nothing
    bl_end = first_pulse === nothing ? min(9, N) : first_pulse - 1
    bl_end < 1 && (bl_end = min(9, N))
    V_baseline = mean(skipmissing(filter(!isnan, V[1:bl_end])))
    baseline_I = mean(skipmissing(filter(!isnan, I_raw[1:bl_end])))
    Vc = V .- V_baseline
    I = I_raw .- baseline_I

    polarity = Int8.(sign.(Vc))
    pulse_idx = ones(Int, N)
    I_FE = copy(I)
    Q_FE = similar(I)

    dt = [zero(t[1]); diff(t)]
    dQ = I_FE .* dt
    q = 0.0
    for i in eachindex(I_FE)
        q += dQ[i]
        Q_FE[i] = q
    end

    df[!, :current] .= I
    return hcat(df, DataFrame(polarity=polarity, pulse_idx=pulse_idx, I_FE=I_FE, Q_FE=Q_FE))
end

"""
Analyze a waveform that may contain ordinary PUND, PN, or PN followed by PUND.
"""
function analyze_pund_and_pn(df::DataFrame; DEBUG::Bool=false)
    @assert all(["time", "voltage", "current"] .∈ Ref(names(df))) "columns :time, :voltage, :current must exist"
    det = detect_pund_pulses(df.time, df.voltage, df.current)
    pulses = det.pulses
    n = length(pulses)

    if n == 6
        # Split the waveform between PN and PUND but keep the pre-pulse baseline
        # at the start of each segment so baseline subtraction works downstream.
        split_row = last(pulses[2])
        pn_df = DataFrame(df[1:split_row, :])
        pund_df, pund_pulses = _slice_signal_dataframe(df, pulses[3:6]; from=split_row + 1)
        return (
            kind=:pn_pund,
            pn=analyze_pn(pn_df; DEBUG=DEBUG),
            pund=analyze_pund(pund_df; pulses=pund_pulses, has_initial_poling=false, DEBUG=DEBUG),
            pn_group_size=1,
            pund_group_size=4,
            detection=det,
        )
    elseif n == 5
        return (
            kind=:pund,
            pn=nothing,
            pund=analyze_pund(df; pulses=pulses, DEBUG=DEBUG),
            pn_group_size=2,
            pund_group_size=5,
            detection=det,
        )
    elseif n == 4
        return (
            kind=:pund,
            pn=nothing,
            pund=analyze_pund(df; pulses=pulses, has_initial_poling=false, DEBUG=DEBUG),
            pn_group_size=2,
            pund_group_size=4,
            detection=det,
        )
    elseif n == 2
        return (
            kind=:pn,
            pn=analyze_pn(df; DEBUG=DEBUG),
            pund=nothing,
            pn_group_size=1,
            pund_group_size=4,
            detection=det,
        )
    end

    if DEBUG
        return (
            kind=:unknown,
            pn=nothing,
            pund=analyze_pund(df; pulses=pulses, DEBUG=true),
            pn_group_size=2,
            pund_group_size=5,
            detection=det,
        )
    end
    error("Unsupported PUND/PN waveform: detected $n pulses")
end

"""
Analyze PUND fatigue across multiple entries, producing:
- traces: vector of NamedTuples (cycles, x, y, label, file_index, rep_index, rep_count) for top-panel overlapped curves
- x_label, y_label: suggested axis labels for top panel
- pr_points: vector of NamedTuples (cycles, Pr, Ec_plus, Ec_minus, C_F, tau_s, R_Ohm, file_index, rep_index, rep_count) for bottom panel

Input:
entries :: Vector of NamedTuples with fields:
    kind::Symbol              # :pund
    df::DataFrame             # raw data frame for the entry
    params::Dict{Symbol,Any}  # device-level params (expects :area_um2, :t_HZO_nm)
    timestamp::Any            # optional; used for chronological sorting

Output pr_points fields:
- cycles: fatigue cycle number
- Pr: remnant polarization (μC/cm²)
- Ec_plus/Ec_minus: coercive field (V or MV/cm) at positive/negative zero crossing
- C_F: capacitance in Farads, extracted from U/D pulse RC fit
- tau_s: RC time constant in seconds
- R_Ohm: series resistance in Ohms

Notes:
- Only the exact parameter names are considered for geometry: area_um2 (μm²) and t_HZO_nm (nm).
- Each PUND entry is analyzed and each repetition (quintuple) is treated as one fatigue cycle.
- Remnant polarization is computed as Pr = (max(P) - min(P)) / 2 from each repetition's P–E trace (only if area_um2 is present).
- Coercive field: Ec_plus/Ec_minus are the x-values (V or MV/cm) where P-E crosses zero during positive/negative switching.
- RC parameters are extracted from U/D (non-switching) pulses by fitting I(t) = I_ss * (1 - exp(-t/τ)) during voltage ramps.
"""
function analyze_pund_fatigue_combined(entries::Vector)
    # Results
    traces = NamedTuple{(:cycles, :x, :y, :label, :file_index, :rep_index, :rep_count),
        Tuple{Int,Vector{Float64},Vector{Float64},String,Int,Int,Int}}[]
    # pr_points includes: Pr, Ec, and RC parameters (C_F, tau_s, R_Ohm)
    pr_points = NamedTuple{(:cycles, :Pr, :Ec_plus, :Ec_minus, :C_F, :tau_s, :R_Ohm, :file_index, :rep_index, :rep_count),
        Tuple{Float64,Float64,Union{Float64,Nothing},Union{Float64,Nothing},
            Union{Float64,Nothing},Union{Float64,Nothing},Union{Float64,Nothing},Int,Int,Int}}[]

    # Axis label decisions
    any_thickness = false
    all_thickness = true
    any_area = false
    all_area = true

    # Running counters
    cycles_accum = 0
    file_counter = 0

    # Cached pulse boundaries — reused across fatigue cycles with identical waveforms
    cached_pulses::Union{Nothing,Vector{UnitRange{Int}}} = nothing

    # Sort chronologically when possible (by :timestamp or params[:timestamp]); fallback to given order
    sorted_entries = sort(entries; by=e -> begin
        if hasproperty(e, :timestamp)
            return e.timestamp
        elseif hasproperty(e, :params) && haskey(e.params, :timestamp)
            return e.params[:timestamp]
        else
            return 0
        end
    end)

    for entry in sorted_entries
        kind = hasproperty(entry, :kind) ? entry.kind : :pund
        df = hasproperty(entry, :df) ? entry.df : DataFrame()
        params = hasproperty(entry, :params) ? entry.params : Dict{Symbol,Any}()

        if kind === :pund
            nrows = nrow(df)
            nrows == 0 && continue

            # Device parameters (exact keys expected)
            area_um2 = haskey(params, :area_um2) ? Float64(params[:area_um2]) : NaN
            thickness_nm = haskey(params, :t_HZO_nm) ? Float64(params[:t_HZO_nm]) : NaN

            # Track availability for axis labels
            if isfinite(thickness_nm) && thickness_nm > 0
                any_thickness = true
            else
                all_thickness = false
            end
            if isfinite(area_um2) && area_um2 > 0
                any_area = true
            else
                all_area = false
            end

            # Analyze single PUND to get FE charge, pulse indices, etc.
            # Detect pulses on first cycle, reuse for subsequent cycles
            if cached_pulses === nothing
                det = detect_pund_pulses(df.time, df.voltage, df.current)
                cached_pulses = det.pulses
            end
            df_an = analyze_pund(df; pulses=cached_pulses)

            # Y-axis base (centered) and unit handling
            Q_FE = df_an.Q_FE
            finite_Q = filter(!isnan, Q_FE)
            q_center = isempty(finite_Q) ? 0.0 : mean(finite_Q)
            Qc = Q_FE .- q_center

            y_all = if isfinite(area_um2) && area_um2 > 0
                area_cm2 = area_um2 / 1e8      # μm² -> cm²
                (Qc ./ area_cm2) .* 1e6        # μC/cm²
            else
                Qc .* 1e12                     # pC
            end

            # X-axis base (E if thickness given, else V)
            V = df_an.voltage
            x_all = if isfinite(thickness_nm) && thickness_nm > 0
                (V ./ (thickness_nm * 1e-7)) ./ 1e6   # V/(nm->cm) => V/cm => MV/cm
            else
                V
            end

            # Per-PUND repetition handling: each quintuple (poling,P,U,N,D) is one cycle increment
            pid = df_an.pulse_idx
            maxpid = maximum(pid)
            n_groups = maxpid ÷ 5
            file_counter += 1

            # Analyze U/D pulses for RC parameters
            ud_analysis = analyze_ud_pulses(df_an.time, df_an.voltage, df_an.current, cached_pulses, n_groups)

            # If this entry has an explicit fatigue_cycle, use it directly
            explicit_cycle = haskey(params, :fatigue_cycle) ? Int(params[:fatigue_cycle]) : nothing

            for rep in 1:n_groups
                # Use only P and N pulses within this repetition
                P_code = rep * 5 - 3
                N_code = rep * 5 - 1
                rep_mask = (pid .== P_code) .| (pid .== N_code)

                valid = rep_mask .& isfinite.(x_all) .& isfinite.(y_all) .& .!isnan.(y_all)
                if any(valid)
                    cyc = explicit_cycle !== nothing ? explicit_cycle : cycles_accum + 1

                    xv = collect(x_all[valid])
                    yv = collect(y_all[valid])

                    # Top panel trace
                    push!(traces, (cycles=cyc, x=xv, y=yv, label=string(cyc),
                        file_index=file_counter, rep_index=rep, rep_count=n_groups))

                    # Extract P and N pulse data separately for coercive field
                    P_valid = (pid .== P_code) .& isfinite.(x_all) .& isfinite.(y_all) .& .!isnan.(y_all)
                    N_valid = (pid .== N_code) .& isfinite.(x_all) .& isfinite.(y_all) .& .!isnan.(y_all)
                    x_P, y_P = collect(x_all[P_valid]), collect(y_all[P_valid])
                    x_N, y_N = collect(x_all[N_valid]), collect(y_all[N_valid])

                    Ec_plus = find_zero_crossing(x_P, y_P)
                    Ec_minus = find_zero_crossing(x_N, y_N)

                    # Bottom panel: Pr = (max(P) - min(P)) / 2 (only if area was provided)
                    if isfinite(area_um2) && area_um2 > 0
                        ymin = minimum(yv)
                        ymax = maximum(yv)
                        if isfinite(ymin) && isfinite(ymax)
                            Pr = 0.5 * (ymax - ymin)
                            if isfinite(Pr)
                                # Get RC parameters for this repetition (average U and D fits)
                                C_F, tau_s, R_Ohm = nothing, nothing, nothing
                                if rep <= length(ud_analysis.groups)
                                    grp = ud_analysis.groups[rep]
                                    C_vals = filter(!isnothing, [grp.C_U, grp.C_D])
                                    tau_vals = filter(!isnothing, [grp.tau_U, grp.tau_D])
                                    R_vals = filter(!isnothing, [grp.R_U, grp.R_D])
                                    C_F = isempty(C_vals) ? nothing : mean(C_vals)
                                    tau_s = isempty(tau_vals) ? nothing : mean(tau_vals)
                                    R_Ohm = isempty(R_vals) ? nothing : mean(R_vals)
                                end

                                push!(pr_points, (cycles=float(cyc), Pr=Pr,
                                    Ec_plus=Ec_plus, Ec_minus=Ec_minus,
                                    C_F=C_F, tau_s=tau_s, R_Ohm=R_Ohm,
                                    file_index=file_counter, rep_index=rep, rep_count=n_groups))
                            end
                        end
                    end

                    # Update cumulative cycle counter
                    cycles_accum = explicit_cycle !== nothing ? explicit_cycle : cycles_accum + 1
                end
            end
        end
    end

    # Decide axis labels
    x_label = all_thickness ? "E (MV/cm)" : "Voltage (V)"
    y_label = all_area ? "Polarization (μC/cm²)" : "Switching charge (pC)"

    return Dict(
        :traces => traces,
        :x_label => x_label,
        :y_label => y_label,
        :pr_points => pr_points,
    )
end
