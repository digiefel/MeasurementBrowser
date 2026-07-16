using CSV
using DataFrames

"""
Read one IVSweep-format file. IV sweeps, four-terminal IV, TLM 4-point, and breakdown all share
this layout: '#'-prefixed metadata, one column-header line, then the data table.
"""
function read_iv_sweep(filepath::AbstractString)::DataFrame
    df = CSV.read(filepath, DataFrame; comment="#", ntasks=1)
    cols = names(df)
    n = nrow(df)
    if "Voltage_V" in cols                          # voltage-driven IV / breakdown layout
        v = Float64.(df[!, "Voltage_V"])
        i = Float64.(df[!, "Current_High_A"])
        i_low = ("Current_Low_A" in cols) ? Float64.(df[!, "Current_Low_A"]) : fill(NaN, n)
    else
        i = Float64.(df[!, "Current_A"])
        vlow = ("VoltageLow_V" in cols) ? Float64.(df[!, "VoltageLow_V"]) : zeros(n)
        v = Float64.(df[!, "VoltageHigh_V"]) .- vlow
        i_low = fill(NaN, n)
    end
    time_s = ("Time_sec" in cols) ? Float64.(df[!, "Time_sec"]) : fill(NaN, n)
    return DataFrame(v=v, i=i, i_low=i_low, time_s=time_s)
end
