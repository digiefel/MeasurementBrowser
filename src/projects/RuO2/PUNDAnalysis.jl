using DataFrames: DataFrame, nrow
using SmoothData: smoothdata
using Statistics: mean

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

    if n == 7
        split_row = last(pulses[2])
        pn_df = DataFrame(df[1:split_row, :])
        pund_df, pund_pulses = _slice_signal_dataframe(df, pulses[3:7]; from=split_row + 1)
        return (
            kind=:pn_pund,
            pn=analyze_pn(pn_df; DEBUG=DEBUG),
            pund=analyze_pund(pund_df; pulses=pund_pulses, DEBUG=DEBUG),
            pn_group_size=1,
            pund_group_size=5,
            detection=det,
        )
    elseif n == 6
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
- pr_points: vector of NamedTuples (cycles, Pr, Ec_plus, Ec_minus, file_index, rep_index, rep_count) for bottom panel

Input:
entries :: Vector of NamedTuples with fields:
    kind::Symbol              # :pund
    df::DataFrame             # raw data frame for the entry
    params::Dict{Symbol,Any}  # merged device and measurement params; requires :fatigue_count
    timestamp::Any            # used for chronological sorting

Output pr_points fields:
- cycles: fatigue cycle number
- Pr: remnant polarization (μC/cm²)
- Ec_plus/Ec_minus: coercive field (V or MV/cm) at positive/negative zero crossing

Notes:
- Only the exact parameter names are considered for geometry: area_um2 (μm²) and t_HZO_nm (nm).
- Each PUND entry must already identify its fatigue count; zero is the unfatigued baseline.
- Remnant polarization is computed as Pr = (max(P) - min(P)) / 2 from each repetition's P–E trace (only if area_um2 is present).
- Coercive field: Ec_plus/Ec_minus are the x-values (V or MV/cm) where P-E crosses zero during positive/negative switching.
"""
function analyze_pund_fatigue_combined(entries::Vector)
    # Results
    traces = NamedTuple{(:cycles, :x, :y, :label, :file_index, :rep_index, :rep_count),
        Tuple{Int,Vector{Float64},Vector{Float64},String,Int,Int,Int}}[]
    pr_points = NamedTuple{(:cycles, :Pr, :Ec_plus, :Ec_minus, :file_index, :rep_index, :rep_count),
        Tuple{Float64,Float64,Union{Float64,Nothing},Union{Float64,Nothing},
            Int,Int,Int}}[]

    # Axis label decisions
    all_thickness = true
    all_area = true

    file_counter = 0

    # Cached pulse boundaries — reused across fatigue cycles with identical waveforms
    cached_pulses::Union{Nothing,Vector{UnitRange{Int}}} = nothing

    sorted_entries = sort(entries; by=e -> e.timestamp)

    for entry in sorted_entries
        kind = entry.kind
        df = entry.df
        params = entry.params

        kind === :pund || error("PUND fatigue analysis only supports :pund entries")
        nrow(df) > 0 || error("PUND fatigue analysis received an empty dataframe")

        # Device parameters (exact keys expected)
        area_um2 = haskey(params, :area_um2) ? Float64(params[:area_um2]) : NaN
        thickness_nm = haskey(params, :t_HZO_nm) ? Float64(params[:t_HZO_nm]) : NaN

        all_thickness &= isfinite(thickness_nm) && thickness_nm > 0
        all_area &= isfinite(area_um2) && area_um2 > 0

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

        pid = df_an.pulse_idx
        maxpid = maximum(pid)
        n_groups = maxpid ÷ 5
        file_counter += 1

        fatigue_count = Int(params[:fatigue_count])

        for rep in 1:n_groups
            # Use only P and N pulses within this repetition
            P_code = rep * 5 - 3
            N_code = rep * 5 - 1
            rep_mask = (pid .== P_code) .| (pid .== N_code)

            valid = rep_mask .& isfinite.(x_all) .& isfinite.(y_all) .& .!isnan.(y_all)
            if any(valid)
                cyc = fatigue_count

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
                            push!(pr_points, (cycles=float(cyc), Pr=Pr,
                                Ec_plus=Ec_plus, Ec_minus=Ec_minus,
                                file_index=file_counter, rep_index=rep, rep_count=n_groups))
                        end
                    end
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
