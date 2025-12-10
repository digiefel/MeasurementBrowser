export read_fe_pund, read_wakeup

"""
Read FE PUND data from CSV file
"""
function read_fe_pund(filename, workdir=".")
    filepath = joinpath(workdir, filename)
    lines = readlines(filepath)
    data_start = 1

    for (i, line) in enumerate(lines)
        if occursin("Time,MeasResult1_value,MeasResult2_value", line)
            data_start = i + 1
            break
        end
    end

    if data_start == 1
        return DataFrame()
    end

    data_lines = lines[data_start:end]
    time = Float64[]
    current = Float64[]
    voltage = Float64[]
    current_time = Float64[]
    voltage_time = Float64[]

    for line in data_lines
        if !isempty(strip(line))
            parts = split(line, ',')
            if length(parts) >= 5
                try
                    push!(time, parse(Float64, parts[1]))
                    push!(current, parse(Float64, parts[2]))
                    push!(voltage, parse(Float64, parts[3]))
                    push!(current_time, parse(Float64, parts[4]))
                    push!(voltage_time, parse(Float64, parts[5]))
                catch
                    continue
                end
            end
        end
    end

    return DataFrame(time=time, current=current, voltage=voltage,
        current_time=current_time, voltage_time=voltage_time)
end

"""
Read wakeup data from CSV file - counts lines and extracts amplitude from filename
Returns DataFrame with pulse_count and amplitude
"""
function read_wakeup(filename, workdir=".")
    filepath = joinpath(workdir, filename)
    lines = readlines(filepath)

    # Find the header line that contains "Time,MeasResult1_value,MeasResult2_value"
    data_start = 1
    for (i, line) in enumerate(lines)
        if occursin("Time,MeasResult1_value,MeasResult2_value", line)
            data_start = i + 1
            break
        end
    end

    # Count actual data lines (non-empty lines with commas after the header)
    data_lines = 0
    for line in lines[data_start:end]
        if !isempty(strip(line)) && occursin(',', line)
            data_lines += 1
        end
    end

    # Extract amplitude from filename (pattern like "3V", "10V", etc.)
    amplitude_match = match(r"(\d+(?:\.\d+)?)V", filename)
    amplitude = amplitude_match !== nothing ? parse(Float64, amplitude_match.captures[1]) : 0.0

    return DataFrame(pulse_count=data_lines, amplitude=amplitude)
end
