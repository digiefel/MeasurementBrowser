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

    # offset current data to zero
    n0 = min(10, N)
    baseline_I = mean(skipmissing(filter(!isnan, I[1:n0])))
    if DEBUG
        @info "analyze_pund: baseline current offset" n0 = n0 baseline = baseline_I
    end
    I .-= baseline_I

    # ---- helper: contiguous true–runs --------------------------------------------------
    function true_runs(mask::BitVector)
        runs = UnitRange{Int}[]
        start = nothing
        for i ∈ eachindex(mask)
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

    # ---- pulse detection ---------------------------------------------------------------
    # Use smoothed derivative detection to capture full triangular pulses
    dV = smoothdata([0.0; diff(V)], :movmedian, 9)
    baseline_dV = std(dV[1:min(10, length(dV) ÷ 10)])
    dV_threshold = baseline_dV * 5
    # Find regions with significant voltage changes
    pulse_mask = abs.(dV) .> dV_threshold

    # expand pulse regions to capture full triangular waves
    expanded_mask = copy(pulse_mask)
    safe_win = 5  # Larger expansion window
    for i in (safe_win+1):(length(pulse_mask)-safe_win)
        if any(pulse_mask[(i-safe_win):(i+safe_win)])
            expanded_mask[i] = true
        end
    end
    all_pulses = true_runs(expanded_mask)

    # Filter out short pulses and pulses with small voltage amplitudes
    baseline_V = mean(abs.(V[1:min(9, length(V))]))
    min_pulse_length = 20  # minimum points for a valid pulse
    min_voltage_amplitude = baseline_V * 5  # minimum voltage amplitude

    pulses = UnitRange{Int}[]
    for pulse in all_pulses
        if length(pulse) >= min_pulse_length
            pulse_V_range = maximum(abs.(V[pulse]))
            if pulse_V_range >= min_voltage_amplitude
                push!(pulses, pulse)
            end
        end
    end

    @assert !isempty(pulses) "no valid pulses found; adjust derivative threshold or filtering parameters"

    # ---- polarity alignment --------------------------------------------------------
    mismatches = 0
    total = 0
    for pulse in pulses
        # take only the first half of the pulse in case it is dominated by I=dV/dt
        # (and thus mean(I) ~ 0)
        V_avg, I_avg = mean(V[pulse[1:end÷2]]), mean(I[pulse[1:end÷2]])
        if abs(V_avg) > 0 && abs(I_avg) > 0
            total += 1
            mismatches += sign(V_avg) != sign(I_avg)
        end
    end

    if total > 0 && mismatches == total
        I = -I
    elseif mismatches > 0
        error("Inconsistent polarity: $(mismatches)/$(total) pulses misaligned")
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

            for rep in 1:n_groups
                # Use only P and N pulses within this repetition
                P_code = rep * 5 - 3
                N_code = rep * 5 - 1
                rep_mask = (pid .== P_code) .| (pid .== N_code)

                valid = rep_mask .& isfinite.(x_all) .& isfinite.(y_all) .& .!isnan.(y_all)
                if any(valid)
                    cyc = cycles_accum + 1

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

                    # Increment cumulative cycles after this repetition
                    cycles_accum += 1
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
