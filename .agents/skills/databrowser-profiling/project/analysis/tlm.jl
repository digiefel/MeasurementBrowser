using DataFrames: DataFrame, combine, groupby, nrow
using Printf: @sprintf
using Statistics: mean

tlm_geometry(params::AbstractDict)::Tuple{Float64,Float64} = (
    Float64(get(params, :length_um, NaN)),
    Float64(get(params, :width_um, NaN)),
)

function tlm_combined_table(
    measurements,
    data;
    Vmin::Real=0.0003,
    Imin::Real=1e-15,
)::DataFrame
    length(measurements) == length(data) || error("Measurement and data counts must match")
    combined = DataFrame()

    for (measurement, df) in zip(measurements, data)
        params = metadata(measurement)
        length_um, width_um = tlm_geometry(params)
        isfinite(length_um) && isfinite(width_um) && length_um > 0 && width_um > 0 || continue

        resistance = [
            abs(i) < Imin || abs(v) < Vmin ? NaN : v / i
            for (i, v) in zip(df.i, df.v)
        ]
        rows = nrow(df)
        append!(
            combined,
            DataFrame(
                filepath=fill(String(params[:filepath]), rows),
                length_um=fill(length_um, rows),
                width_um=fill(width_um, rows),
                resistance_ohm=resistance,
                resistance_normalized=resistance .* width_um,
                i=df.i,
                v=df.v,
                device_name=fill(@sprintf("L%.4gW%.2g", length_um, width_um), rows),
            ),
        )
    end

    return combined
end

function fit_tlm_sheet_resistance(table::DataFrame)::NamedTuple
    nrow(table) == 0 && return (R_sheet=NaN, R_cprime=NaN, rho_c=NaN, r_squared=NaN)

    grouped = combine(
        groupby(table, [:length_um, :width_um]),
        :resistance_ohm => (values -> mean(filter(isfinite, values))) => :R_ohm,
    )
    valid = isfinite.(grouped.R_ohm) .&
        isfinite.(grouped.length_um) .&
        isfinite.(grouped.width_um) .&
        (grouped.width_um .> 0)
    data = grouped[valid, :]
    nrow(data) < 2 && return (R_sheet=NaN, R_cprime=NaN, rho_c=NaN, r_squared=NaN)

    x = data.length_um .* 1e-4
    width_cm = data.width_um .* 1e-4
    y = data.R_ohm .* width_cm
    n = length(x)
    sx = sum(x)
    sy = sum(y)
    sxx = sum(x .^ 2)
    sxy = sum(x .* y)
    denom = n * sxx - sx^2
    abs(denom) < 1e-12 && return (R_sheet=NaN, R_cprime=NaN, rho_c=NaN, r_squared=NaN)

    R_sheet = (n * sxy - sx * sy) / denom
    two_R_cprime = (sy - R_sheet * sx) / n
    R_cprime = 0.5 * two_R_cprime
    predicted = two_R_cprime .+ R_sheet .* x
    residual = sum((y .- predicted) .^ 2)
    spread = sum((y .- mean(y)) .^ 2)
    r_squared = spread == 0 ? NaN : 1 - residual / spread
    rho_c = R_cprime^2 / R_sheet

    return (; R_sheet, R_cprime, rho_c, r_squared)
end
