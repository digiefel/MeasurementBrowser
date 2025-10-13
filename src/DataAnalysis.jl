module Analysis

using DataFrames
using Dates
using Statistics, SmoothData
using Printf

export analyze_breakdown, analyze_pund, extract_tlm_geometry_from_params, analyze_tlm_combined, analyze_pund_fatigue_combined

"""
Simple breakdown analysis for I-V data
"""
function analyze_breakdown(df, threshold_current=1e-6)
    if nrow(df) == 0
        return Dict()
    end

    abs_current = abs.(df.current1)
    breakdown_idx = findfirst(x -> x > threshold_current, abs_current)
    breakdown_voltage = breakdown_idx !== nothing ? df.voltage[breakdown_idx] : NaN
    leakage_current = breakdown_idx !== nothing ? mean(abs_current[1:max(1, breakdown_idx - 5)]) : mean(abs_current)

    return Dict(
        "breakdown_voltage" => breakdown_voltage,
        "leakage_current" => leakage_current,
        "max_current" => maximum(abs_current),
        "voltage_range" => (minimum(df.voltage), maximum(df.voltage))
    )
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

    # Expand pulse regions to capture full triangular waves
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
Validate TLM measurement dataframe for required columns and data quality

Returns (is_valid, issues) where issues is a vector of warning strings
"""
function validate_tlm_dataframe(df::DataFrame, filepath::String="")
    issues = String[]

    # Check required columns
    required_cols = ["current_source", "v_gnd"]
    for col in required_cols
        if !(col in names(df))
            push!(issues, "Missing required column: $col")
        end
    end

    if !isempty(issues)
        return (false, issues)
    end

    # Check for empty data
    if nrow(df) == 0
        push!(issues, "Empty dataframe")
        return (false, issues)
    end

    # Check for all-zero currents (would cause division by zero)
    if all(df.current_source .== 0)
        push!(issues, "All source currents are zero - cannot calculate resistance")
    end

    # Check for excessive NaN/missing values
    nan_fraction = sum(ismissing.(df.current_source) .| isnan.(df.current_source)) / nrow(df)
    if nan_fraction > 0.5
        push!(issues, "More than 50% of current data is missing or NaN")
    end

    voltage_nan_fraction = sum(ismissing.(df.v_gnd) .| isnan.(df.v_gnd)) / nrow(df)
    if voltage_nan_fraction > 0.5
        push!(issues, "More than 50% of voltage data is missing or NaN")
    end

    return (isempty(issues), issues)
end

"""
Calculate sheet resistance and specific contact resistivity (long-contact TLM).

Given analyzed TLM data with columns (:length_um, :width_um, :resistance_ohm),
fits the width-invariant relation  R*W = 2R_c' + R_sheet * L
and returns (sheet_resistance_ohm_per_square,
             contact_resistance_per_width_ohm_cm,
             contact_resistivity_ohm_cm2,
             r_squared).

Notes:
- Uses long-contact approximation (L_c ≫ L_T):  ρ_c = (R_c')^2 / R_sheet.
- Converts μm → cm internally. Width must be > 0.
"""
function calculate_sheet_resistance(analysis_df::DataFrame)
    if nrow(analysis_df) == 0
        @warn "Empty analysis dataframe for sheet/contact calculation"
        return (NaN, NaN, NaN, NaN)
    end

    # Average repeated measurements per geometry
    g = combine(groupby(analysis_df, [:length_um, :width_um]),
        :resistance_ohm => (x -> mean(filter(isfinite, x))) => :R_ohm)

    # Validate inputs
    mask = isfinite.(g.R_ohm) .& isfinite.(g.length_um) .& isfinite.(g.width_um) .& (g.width_um .> 0)
    d = g[mask, :]
    if nrow(d) < 2
        @warn "Need at least 2 valid geometry points"
        return (NaN, NaN, NaN, NaN)
    end

    # Build regression: y = R*W (Ω·cm), x = L (cm)
    L_cm = d.length_um .* 1e-4
    W_cm = d.width_um .* 1e-4
    y = d.R_ohm .* W_cm           # Ω·cm
    x = L_cm                      # cm

    n = length(x)
    sx, sy = sum(x), sum(y)
    sxx, sxy = sum(x .^ 2), sum(x .* y)
    denom = n * sxx - sx^2
    if abs(denom) < 1e-12
        @warn "Cannot fit: insufficient variation in lengths"
        return (NaN, NaN, NaN, NaN)
    end

    R_sheet = (n * sxy - sx * sy) / denom     # Ω/□ (slope)
    two_Rcprime = (sy - R_sheet * sx) / n     # 2 * R_c' (Ω·cm) (intercept)
    Rcprime = 0.5 * two_Rcprime               # Ω·cm

    # R^2 on the y = R*W regression
    y_pred = two_Rcprime .+ R_sheet .* x
    ss_res = sum((y .- y_pred) .^ 2)
    ss_tot = sum((y .- mean(y)) .^ 2)
    R2 = 1 - ss_res / ss_tot

    rho_c = (Rcprime^2) / R_sheet             # Ω·cm²  (long-contact limit)

    return (R_sheet, Rcprime, rho_c, R2)
end

"""
Extract geometry information from device parameters

Returns (length_um, width_um) or (NaN, NaN) if not found
Expects device_params to contain :length_um and :width_um keys
"""
function extract_tlm_geometry_from_params(device_params::Dict{Symbol,Any}, filepath::String="")
    length_um = get(device_params, :length_um, NaN)
    width_um = get(device_params, :width_um, NaN)

    # Try alternative key names
    if isnan(length_um)
        length_um = get(device_params, :length, NaN)
    end
    if isnan(width_um)
        width_um = get(device_params, :width, NaN)
    end

    # If still not found, try fallback filename parsing
    if isnan(length_um) || isnan(width_um)
        @info "Geometry not found in device parameters, trying filename parsing for: $filepath"

        # Remove path and extension
        basename_file = basename(filepath)
        name_part = replace(basename_file, r"\.(csv|txt)$" => "")

        # Pattern to match TLML<length>W<width>
        pattern = r"TLML(\d+)W(\d+)"
        m = match(pattern, name_part)

        if m !== nothing
            length_um = parse(Float64, m.captures[1])
            width_um = parse(Float64, m.captures[2])
            @info "Extracted geometry from filename: L=$(length_um)μm, W=$(width_um)μm"
        else
            @warn "Could not extract geometry from device parameters or filename: $filepath" available_keys = keys(device_params)
            return (NaN, NaN)
        end
    end

    return (Float64(length_um), Float64(width_um))
end

"""
Analyze multiple TLM measurements for combined plotting

Takes a vector of (filepath, dataframe, device_params) tuples and returns a combined analysis DataFrame
suitable for plotting width-normalized resistance vs length.

device_params should contain :length_um and :width_um keys with geometry information.

Returns DataFrame with columns:
- filepath: original file path
- length_um: extracted length in micrometers
- width_um: extracted width in micrometers
- resistance_ohm: calculated resistance (V/I at each current point)
- resistance_normalized: resistance * width (Ω·μm)
- current_source: source current values
- voltage: measured voltage values
"""
function analyze_tlm_combined(files_data_params::Vector{Tuple{String,DataFrame,Dict{Symbol,Any}}}; Vmin=0.0003, Imin=1e-15)
    if isempty(files_data_params)
        @warn "No TLM data provided for combined analysis"
        return DataFrame()
    end

    combined_data = DataFrame()
    processed_files = 0

    for (filepath, df, device_params) in files_data_params
        # Validate the dataframe
        is_valid, issues = validate_tlm_dataframe(df, filepath)
        if !is_valid
            @warn "Skipping invalid TLM file: $filepath" issues = issues
            continue
        end

        if !isempty(issues)
            @warn "TLM data quality issues in $filepath" issues = issues
        end

        # Extract geometry from device parameters
        length_um, width_um = extract_tlm_geometry_from_params(device_params, filepath)

        if isnan(length_um) || isnan(width_um) || length_um <= 0 || width_um <= 0
            @warn "Skipping file with invalid geometry: $filepath (L=$length_um, W=$width_um)"
            continue
        end

        # Calculate resistance from V/I, handling division by zero
        resistance_ohm = similar(df.current_source, Float64)
        for i in eachindex(df.current_source)
            if abs(df.current_source[i]) < Imin  # Avoid division by very small numbers
                resistance_ohm[i] = NaN
            elseif abs(df.v_gnd[i]) < Vmin # avoid inaccurate values
                resistance_ohm[i] = NaN
            else
                resistance_ohm[i] = df.v_gnd[i] / df.current_source[i]
            end
        end

        # Width-normalize the resistance
        resistance_normalized = resistance_ohm .* width_um

        # Create a dataframe for this file
        file_data = DataFrame(
            filepath=fill(filepath, nrow(df)),
            length_um=fill(length_um, nrow(df)),
            width_um=fill(width_um, nrow(df)),
            resistance_ohm=resistance_ohm,
            resistance_normalized=resistance_normalized,
            current_source=df.current_source,
            voltage=df.v_gnd
        )

        # Add device name for plotting
        # device_name = "L$(Int(length_um))W$(Int(width_um))"
        device_name = @sprintf "L%.4gW%.2g" length_um width_um
        file_data.device_name = fill(device_name, nrow(df))

        # Append to combined data
        combined_data = vcat(combined_data, file_data)
        processed_files += 1
    end

    if nrow(combined_data) == 0
        @warn "No valid TLM data after processing all files" attempted_files = length(files_data_params)
    else
        @info "TLM combined analysis completed" n_files = processed_files n_total_files = length(files_data_params) n_points = nrow(combined_data)

        # Calculate and report sheet resistance if we have enough data
        R_sheet, R_cprime, rho_c, R2 = calculate_sheet_resistance(combined_data)
        if isfinite(R_sheet)
            @info "Sheet/contact analysis" sheet_resistance_ohm_per_sq = R_sheet contact_resistance_per_width_ohm_cm = R_cprime contact_resistivity_ohm_cm2 = rho_c r_squared = R2
        end
    end

    return combined_data
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
    params::Dict{Symbol,Any}  # device-level params (expects :area_um2, :thickness_nm)
    timestamp::Any            # optional; used for chronological sorting

Notes:
- Only the exact parameter names are considered for geometry: area_um2 (μm²) and thickness_nm (nm).
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
            thickness_nm = haskey(params, :thickness_nm) ? Float64(params[:thickness_nm]) : NaN

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

end # module
