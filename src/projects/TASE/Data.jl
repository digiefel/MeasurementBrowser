"""Read one TASE four-terminal I-V source file as source current `:i` and measured voltage `:v`."""
function read_tase_four_terminal_iv(
    filename::AbstractString,
    workdir::AbstractString=".",
)::DataFrame
    filepath = joinpath(workdir, filename)
    lines = readlines(filepath)
    data_start = 1
    has_header = false

    for (i, line) in enumerate(lines)
        if occursin("Current_A,VoltageHigh_V,VoltageLow_V", line)
            data_start = i + 1
            has_header = true
            break
        end
        if occursin(r"^-?\d+\.?\d*[eE]?-?\d*,(-?\d+\.?\d*[eE]?-?\d*,){2,}", line)
            data_start = i
            break
        end
    end

    current = Float64[]
    voltage = Float64[]
    for line in lines[data_start:end]
        isempty(strip(line)) && continue
        parts = split(line, ',')
        try
            push!(current, parse(Float64, parts[1]))
            if has_header
                push!(voltage, parse(Float64, parts[2]) - parse(Float64, parts[3]))
            else
                values = filter(part -> !isempty(strip(part)), parts)
                push!(voltage, parse(Float64, values[2]))
            end
        catch
            continue
        end
    end

    return DataFrame(i=current, v=voltage)
end
