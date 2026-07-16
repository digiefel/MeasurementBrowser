using CSV
using DataFrames

"""
Normalize a RuO2 CV export to frequency, bias, capacitance, impedance, and status. CV files carry
'#'-prefixed procedure metadata before the data table, skipped in one pass via `comment="#"`.
Impedance is read from `Z (Ohm)`, or derived from `G (S)` / `Rp (Ohm)` when only those are present.
A file whose schema lacks the expected columns errors out rather than being silently skipped.
"""
function read_cv_sweep(filename, workdir=".")
    filepath = joinpath(workdir, filename)
    df = CSV.read(filepath, DataFrame; comment="#", ntasks=1)
    columns = names(df)

    f = Float64.(df[!, "Frequency_Hz"])
    bias = Float64.(df[!, "Bias_V"])
    cp = Float64.(df[!, "Cp (F)"])
    t = Float64.(df[!, "Time_sec"])
    z = if "Z (Ohm)" in columns
        Float64.(df[!, "Z (Ohm)"])
    elseif "G (S)" in columns
        _z_from_admittance(f, cp, Float64.(df[!, "G (S)"]))
    elseif "Rp (Ohm)" in columns
        _z_from_admittance(f, cp, 1.0 ./ Float64.(df[!, "Rp (Ohm)"]))
    else
        error("CVSweep $filepath has no Z, G, or Rp column")
    end

    n = nrow(df)
    status_cp = "Status_Cp" in columns ? df[!, "Status_Cp"] : zeros(Int, n)
    status_combined = "Status_Combined" in columns ? df[!, "Status_Combined"] : zeros(Int, n)

    return DataFrame(
        frequency_Hz=f,
        bias_V=bias,
        Cp_F=cp,
        Z_Ohm=z,
        time_s=t,
        status_cp=Int.(status_cp),
        status_combined=Int.(status_combined),
    )
end

_z_from_admittance(f, cp, g) = 1.0 ./ sqrt.(g .^ 2 .+ (2π .* f .* cp) .^ 2)
