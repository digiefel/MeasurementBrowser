export read_cv_sweep

"""
Read a RuO2 CVSweep CSV file with an explicit `Frequency_Hz` column.
Returns a NamedTuple with a normalized DataFrame and `secondary_kind`.
"""
function read_cv_sweep(filename, workdir=".")
    filepath = joinpath(workdir, filename)
    header = open(filepath, "r") do io
        eof(io) && error("Empty CVSweep file: $filepath")
        return chomp(readline(io))
    end

    columns = split(header, ',')
    columns[1] == "Frequency_Hz" || error("Unsupported CVSweep without Frequency_Hz header: $filepath")

    required = ("Frequency_Hz", "Bias_V", "Cp (F)", "Time_sec")
    all(col -> col in columns, required) || error("Unsupported CVSweep schema: $filepath")

    secondary_kind = if "G (S)" in columns
        :conductance
    elseif "Rp (Ohm)" in columns
        :parallel_resistance
    else
        error("Unsupported CVSweep secondary column in $filepath")
    end

    df = CSV.read(filepath, DataFrame)
    rename_map = Dict(
        "Frequency_Hz" => :frequency_Hz,
        "Bias_V" => :bias_V,
        "Cp (F)" => :cp_F,
        "Time_sec" => :time_s,
    )
    "Status_Cp" in columns && (rename_map["Status_Cp"] = :status_cp)
    "Status_Combined" in columns && (rename_map["Status_Combined"] = :status_combined)
    if secondary_kind === :conductance
        rename_map["G (S)"] = :secondary_value
        "Status_G" in columns && (rename_map["Status_G"] = :status_secondary)
    else
        rename_map["Rp (Ohm)"] = :secondary_value
        "Status_Rp" in columns && (rename_map["Status_Rp"] = :status_secondary)
    end
    rename!(df, rename_map)
    for status_column in (:status_cp, :status_secondary, :status_combined)
        hasproperty(df, status_column) || (df[!, status_column] = zeros(Int, nrow(df)))
    end

    return (
        df=select(
            df,
            :frequency_Hz,
            :bias_V,
            :cp_F,
            :secondary_value,
            :time_s,
            :status_cp,
            :status_secondary,
            :status_combined,
        ),
        secondary_kind=secondary_kind,
    )
end
