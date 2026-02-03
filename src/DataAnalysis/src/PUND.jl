"""
    detect_pund_pulses(t, V, I; smooth_window=9, dV_thresh_mult=5.0,
        expand_win=5, min_pulse_len=20, min_V_amp_mult=5.0,
        trough_thresh=0.05) -> NamedTuple

Detect triangular voltage pulses in a PUND waveform and check polarity
consistency. Returns a diagnostic NamedTuple with all intermediate results;
never errors on inconsistent polarity.

Detection has two stages (fallback chain):
  1. Derivative-based: finds pulses via smoothed dV/dt threshold + expansion.
     Works when there are flat inter-pulse gaps.
  2. Trough-based (fallback): splits at near-zero |V| regions.
     Used for continuous fatigue waveforms with no flat gaps.
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

    # Smoothed derivative
    baseline_dV = std(diff(V)[1:max(20, length(V) ÷ 10)])
    dV_threshold = baseline_dV * dV_thresh_mult

    dV = smoothdata([0.0; diff(V)], :movmedian, smooth_window)
    baseline_V = mean(abs.(V[1:min(9, length(V))]))
    min_V_amp = baseline_V * min_V_amp_mult

    # Stage 1: derivative-based detection
    pulses = UnitRange{Int}[]
    detection_method = :none
    pulse_mask = abs.(dV) .> dV_threshold
    expanded_mask = copy(pulse_mask)

    if baseline_dV > 1e-12
        detection_method = :derivative
        for i in (expand_win+1):(length(pulse_mask)-expand_win)
            if any(pulse_mask[(i-expand_win):(i+expand_win)])
                expanded_mask[i] = true
            end
        end
        all_pulses = _true_runs(expanded_mask)

        for pulse in all_pulses
            if length(pulse) >= min_pulse_len
                pulse_V_range = maximum(abs.(V[pulse]))
                if pulse_V_range >= min_V_amp
                    push!(pulses, pulse)
                end
            end
        end
    end

    # Stage 2: trough-based fallback (only if stage 1 found nothing)
    if isempty(pulses)
        detection_method = :trough_fallback
        absV = abs.(V)
        V_peak = maximum(absV)
        trough_threshold = V_peak * trough_thresh
        in_trough = absV .< trough_threshold
        i = 1
        while i <= N
            while i <= N && in_trough[i]; i += 1; end
            i > N && break
            pstart = i
            while i <= N && !in_trough[i]; i += 1; end
            pend = i > N ? N : i - 1
            seg = pstart:pend
            if length(seg) >= min_pulse_len
                push!(pulses, seg)
            end
        end
    end

    # Per-pulse polarity check
    polarity_info = NamedTuple{(:pulse_idx, :V_avg, :I_avg, :V_sign, :I_sign, :match),
        Tuple{Int,Float64,Float64,Int,Int,Bool}}[]
    mismatches = 0
    total = 0
    for (pi, pulse) in enumerate(pulses)
        half = pulse[1:end÷2]
        V_avg = mean(V[half])
        I_avg = mean(I_shifted[half])
        if abs(V_avg) > 0 && abs(I_avg) > 0
            total += 1
            sv = Int(sign(V_avg))
            si = Int(sign(I_avg))
            matched = sv == si
            mismatches += !matched
            push!(polarity_info, (pulse_idx=pi, V_avg=V_avg, I_avg=I_avg,
                V_sign=sv, I_sign=si, match=matched))
        else
            push!(polarity_info, (pulse_idx=pi, V_avg=V_avg, I_avg=I_avg,
                V_sign=0, I_sign=0, match=true))
        end
    end

    verdict = if total == 0
        :no_data
    elseif mismatches == 0
        :ok
    elseif mismatches == total
        :all_flipped
    else
        :inconsistent
    end

    n_groups = length(pulses) ÷ 5
    remainder = length(pulses) % 5

    return (;
        t, V, I=I_shifted, dV, dV_threshold, baseline_dV,
        pulse_mask, expanded_mask,
        pulses, detection_method,
        polarity_info, mismatches, total, verdict,
        n_groups, remainder,
        min_V_amp, baseline_V,
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

"""
Perform PUND analysis on the DataFrame

    # 1. identify triangular voltage pulse train
    # 2. get current in P_n and U_n pulses
    # 3. subtract I_Pn - I_Un = I_FEn
    # 4. integrate I_FEn(t) dt = Q_FEn
    # 5. return analysis dataframe with Q_FE in each row

from this, Q_FE(V) or Q_FE(t) can be very easily extracted
"""
function analyze_pund(df::DataFrame; DEBUG::Bool=false)
    @assert all(["time", "voltage", "current"] .∈ Ref(names(df))) "columns :time, :voltage, :current must exist"
    t, V, I = df[!, :time], df[!, :voltage], df[!, :current]
    N = length(t)

    # ---- pulse detection & polarity check via shared helper --------------------
    det = detect_pund_pulses(t, V, I)
    pulses = det.pulses
    I = det.I  # baseline-shifted

    if DEBUG
        n0 = min(10, N)
        baseline_I = mean(skipmissing(filter(!isnan, df[!, :current][1:n0])))
        @info "analyze_pund: baseline current offset" n0 baseline=baseline_I
        @info "analyze_pund: detection" method=det.detection_method n_pulses=length(pulses) verdict=det.verdict
    end

    @assert !isempty(pulses) "no valid pulses found; adjust derivative threshold or filtering parameters"

    # ---- polarity alignment --------------------------------------------------------
    if det.verdict == :all_flipped
        I = -I
    elseif det.verdict == :inconsistent
        error("Inconsistent polarity: $(det.mismatches)/$(det.total) pulses misaligned")
    end

    # ---- consistency check and grouping into quintuples --------------------------------
    groups = [(pulses[i], pulses[i+1], pulses[i+2], pulses[i+3], pulses[i+4]) for i in 1:5:length(pulses)-4]

    for g in groups
        poling, P, U, Np, D = Tuple(g)
        sP, sPol = sign(sum(V[P])), sign(sum(V[poling]))
        @assert sPol == -sP &&
                sign(sum(V[U])) == sP &&
                sign(sum(V[Np])) == -sP &&
                sign(sum(V[D])) == -sP "unexpected pulse ordering"
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
    for (group_idx, (poling, P, U, Np, D)) in enumerate(groups)
        # polarity flags
        polarity[P] .= 1
        polarity[U] .= 1
        polarity[Np] .= -1
        polarity[D] .= -1
        pulse_idx[poling] .= group_idx * 5 - 4
        pulse_idx[P] .= group_idx * 5 - 3
        pulse_idx[U] .= group_idx * 5 - 2
        pulse_idx[Np] .= group_idx * 5 - 1
        pulse_idx[D] .= group_idx * 5

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
    for i ∈ eachindex(I_FE)
        # Only integrate during P and N switching pulses
        if (pulse_idx[i] % 5 == 2 || pulse_idx[i] % 5 == 4)
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
Analyze PUND fatigue across multiple entries (FE_PUND and Wakeup), producing:
- traces: vector of NamedTuples (cycles, x, y, label, file_index, rep_index, rep_count) for top-panel overlapped curves
- x_label, y_label: suggested axis labels for top panel
- pr_points: vector of NamedTuples (cycles, Pr, file_index, rep_index, rep_count) for bottom panel

Input:
entries :: Vector of NamedTuples with fields:
    kind::Symbol              # :pund or :wakeup
    df::DataFrame             # raw data frame for the entry
    params::Dict{Symbol,Any}  # device-level params (expects :area_um2, :t_HZO_nm)
    timestamp::Any            # optional; used for chronological sorting

Notes:
- Only the exact parameter names are considered for geometry: area_um2 (μm²) and t_HZO_nm (nm).
- Wakeup entries contribute their pulse_count to the cumulative cycle count.
- FE_PUND entries are analyzed and each repetition (quintuple) is treated as one fatigue cycle.
- Remnant polarization is computed as Pr = (max(P) - min(P)) / 2 from each repetition's P–E trace (only if area_um2 is present).
"""
function analyze_pund_fatigue_combined(entries::Vector)
    # Results
    traces = NamedTuple{(:cycles, :x, :y, :label, :file_index, :rep_index, :rep_count),
        Tuple{Int,Vector{Float64},Vector{Float64},String,Int,Int,Int}}[]
    pr_points = NamedTuple{(:cycles, :Pr, :file_index, :rep_index, :rep_count),
        Tuple{Float64,Float64,Int,Int,Int}}[]

    # Axis label decisions
    any_thickness = false
    all_thickness = true
    any_area = false
    all_area = true

    # Running counters
    cycles_accum = 0
    file_counter = 0

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

        if kind === :wakeup
            # Add wakeup pulses to cumulative cycles
            pulses = 0
            if "pulse_count" ∈ names(df)
                pulses = Int(df.pulse_count[1])
            elseif haskey(params, :pulse_count)
                pulses = Int(params[:pulse_count])
            end
            cycles_accum += pulses
            continue
        elseif kind === :pund
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
            df_an = analyze_pund(df)

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

                    # Bottom panel: Pr = (max(P) - min(P)) / 2 (only if area was provided)
                    if isfinite(area_um2) && area_um2 > 0
                        ymin = minimum(yv)
                        ymax = maximum(yv)
                        if isfinite(ymin) && isfinite(ymax)
                            Pr = 0.5 * (ymax - ymin)
                            if isfinite(Pr)
                                push!(pr_points, (cycles=float(cyc), Pr=Pr,
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
