"""
Validate TLM measurement dataframe for required columns and data quality

Returns (is_valid, issues) where issues is a vector of warning strings
"""
function validate_tlm_dataframe(df::DataFrame, filepath::String="")
    issues = String[]

    # Check required columns
    required_cols = ["current_source", "voltage_drop"]
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

    voltage_nan_fraction = sum(ismissing.(df.voltage_drop) .| isnan.(df.voltage_drop)) / nrow(df)
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
function extract_tlm_geometry_from_params(device_params::Dict{Symbol,Any})
    length_um = get(device_params, :length_um, NaN)
    width_um = get(device_params, :width_um, NaN)

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
        length_um, width_um = extract_tlm_geometry_from_params(device_params)

        if isnan(length_um) || isnan(width_um) || length_um <= 0 || width_um <= 0
            @warn "Skipping file with invalid geometry: $filepath (L=$length_um, W=$width_um)"
            continue
        end

        # Calculate resistance from V/I, handling division by zero
        resistance_ohm = similar(df.current_source, Float64)
        for i in eachindex(df.current_source)
            if abs(df.current_source[i]) < Imin  # Avoid division by very small numbers
                resistance_ohm[i] = NaN
            elseif abs(df.voltage_drop[i]) < Vmin # avoid inaccurate values
                resistance_ohm[i] = NaN
            else
                resistance_ohm[i] = df.voltage_drop[i] / df.current_source[i]
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
            voltage=df.voltage_drop
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
