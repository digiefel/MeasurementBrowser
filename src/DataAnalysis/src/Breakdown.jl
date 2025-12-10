"""
Simple breakdown analysis for I-V data
"""
function analyze_breakdown(df::DataFrame, threshold_current=1e-6)
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
