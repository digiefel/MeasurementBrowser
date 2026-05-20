using CSV
using DataFrames

const _CV_REQUIRED_BASE = ("Frequency_Hz", "Bias_V", "Cp (F)", "Time_sec")

"""
CV files may start with '#'-prefixed procedure metadata before the data table.
"""
function _cv_data_header_line(filepath::AbstractString)
    open(filepath, "r") do io
        while !eof(io)
            line = chomp(readline(io))
            startswith(line, '#') && continue
            return line
        end
        return nothing
    end
end

"""
Recognize RuO2 CV exports without accepting unrelated capacitor/impedance CSV files.
"""
function cv_sweep_has_schema(filepath::AbstractString)
    header = _cv_data_header_line(filepath)
    header === nothing && return false
    columns = split(header, ',')
    all(col -> col in columns, _CV_REQUIRED_BASE) || return false
    return ("Z (Ohm)" in columns) || ("G (S)" in columns) || ("Rp (Ohm)" in columns)
end

"""
Normalize supported RuO2 CV exports to frequency, bias, capacitance, impedance, and status.
"""
function read_cv_sweep(filename, workdir=".")
    filepath = joinpath(workdir, filename)
    header = _cv_data_header_line(filepath)
    header === nothing && error("Empty CVSweep file: $filepath")
    columns = split(header, ',')
    all(col -> col in columns, _CV_REQUIRED_BASE) ||
        error("Unsupported CVSweep schema: $filepath")

    skipto = _cv_data_skipto(filepath)
    df = CSV.read(filepath, DataFrame; skipto=skipto, header=skipto - 1)

    f = df[!, "Frequency_Hz"]
    bias = df[!, "Bias_V"]
    cp = df[!, "Cp (F)"]
    t = df[!, "Time_sec"]
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
        frequency_Hz=Float64.(f),
        bias_V=Float64.(bias),
        Cp_F=Float64.(cp),
        Z_Ohm=z,
        time_s=Float64.(t),
        status_cp=Int.(status_cp),
        status_combined=Int.(status_combined),
    )
end

_z_from_admittance(f, cp, g) = 1.0 ./ sqrt.(g .^ 2 .+ (2π .* f .* cp) .^ 2)

function _cv_data_skipto(filepath::AbstractString)
    open(filepath, "r") do io
        line_no = 0
        while !eof(io)
            line_no += 1
            line = readline(io)
            startswith(line, '#') || return line_no + 1
        end
        return line_no + 1
    end
end
