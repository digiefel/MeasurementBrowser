module DataLoader

using CSV
using DataFrames
using Dates
using PrecompileTools: @setup_workload, @compile_workload

export find_files, get_file_patterns, read_iv_sweep, read_fe_pund, read_tlm_4p

"""
Find files matching a pattern in the specified directory
"""
function find_files(pattern, workdir=".")
    all_files = readdir(workdir)
    csv_files = filter(f -> endswith(f, ".csv"), all_files)
    return filter(f -> occursin(pattern, f), csv_files)
end


"""
Get standard file patterns for different measurement types
"""
function get_file_patterns()
    return (
        iv_sweep=r"I_V Sweep",
        fe_pund=r"FE PUND",
        tlm_4p=r"TLM_4P",
        breakdown=r"Break.*oxide",
        wakeup=r"Wakeup"
    )
end

"""
Read I-V sweep data from CSV file, skipping header metadata
"""
function read_iv_sweep(filename, workdir=".")
    filepath = joinpath(workdir, filename)
    lines = readlines(filepath)
    data_start = 1

    for (i, line) in enumerate(lines)
        if occursin(r"^-?\d+\.?\d*,-?\d+\.?\d*[eE]?-?\d*,-?\d+\.?\d*[eE]?-?\d*", line)
            data_start = i
            break
        end
    end

    header_cols = split(lines[data_start-1], ',')
    data_lines = lines[data_start:end]
    voltage = Float64[]
    current = Float64[]

    # check if header line contains expected columns
    # it's enough that there's a column containing "V" and one containing "I"
    if !any(occursin.("V", header_cols)) && !any(occursin.("I", header_cols))
        error("Invalid column names")
    end

    v_idx = findfirst(contains.(header_cols, "V"))
    i_idx = findfirst(contains.(header_cols, "I"))

    # parse data lines
    for line in data_lines
        if !isempty(strip(line))
            parts = split(line, ',')
            if length(parts) >= 2
                try
                    push!(voltage, parse(Float64, parts[v_idx]))
                    push!(current, parse(Float64, parts[i_idx]))
                catch
                    continue
                end
            end
        end
    end

    return DataFrame(v=voltage, i=current)
end

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
Read TLM 4-point data from CSV file
"""
function read_tlm_4p(filename, workdir=".")
    filepath = joinpath(workdir, filename)
    lines = readlines(filepath)
    data_start = 1
    is_new_format = false

    for (i, line) in enumerate(lines)
        if occursin("Current_A,VoltageHigh_V,VoltageLow_V", line)
            data_start = i + 1
            is_new_format = true
            break
        end
        if occursin(r"^-?\d+\.?\d*[eE]?-?\d*,(-?\d+\.?\d*[eE]?-?\d*,){2,}", line)
            data_start = i
            break
        end
    end

    data_lines = lines[data_start:end]
    current_source = Float64[]
    v_gnd = Float64[]

    if is_new_format
        for line in data_lines
            if !isempty(strip(line))
                parts = split(line, ',')
                if length(parts) >= 3
                    try
                        curr = parse(Float64, parts[1])
                        v_high = parse(Float64, parts[2])
                        v_low = parse(Float64, parts[3])
                        push!(current_source, curr)
                        push!(v_gnd, v_high - v_low) # Calculate voltage drop
                    catch
                        continue
                    end
                end
            end
        end
    else
        for line in data_lines
            if !isempty(strip(line))
                parts = split(line, ',')
                # Filter out empty parts
                valid_parts = filter(p -> !isempty(strip(p)), parts)

                if length(valid_parts) >= 3
                    try
                        push!(current_source, parse(Float64, valid_parts[1]))
                        push!(v_gnd, parse(Float64, valid_parts[2]))  # voltage
                    catch e
                        println("Error parsing line: $line, error: $e")
                        continue
                    end
                end
            end
        end
    end

    return DataFrame(current_source=current_source, v_gnd=v_gnd)
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

"""
Extract datetime from filename in format: [... ; YYYY-MM-DD HH_MM_SS].csv
Returns DateTime object or nothing if parsing fails
"""
function extract_datetime_from_filename(filename)
    # Look for pattern like "; 2025-08-06 12_20_59]"
    datetime_match = match(r"; (\d{4}-\d{2}-\d{2}) (\d{2})_(\d{2})_(\d{2})", filename)
    if datetime_match !== nothing
        date_part = datetime_match.captures[1]
        hour = datetime_match.captures[2]
        minute = datetime_match.captures[3]
        second = datetime_match.captures[4]

        try
            # Construct full datetime string
            datetime_str = "$(date_part) $(hour):$(minute):$(second)"
            return DateTime(datetime_str, "yyyy-mm-dd HH:MM:SS")
        catch e
            println("Warning: Could not parse datetime from '$filename': $e")
            return nothing
        end
    end
    return nothing
end

@setup_workload begin
    # TODO
    @compile_workload begin
        # TODO
    end
end


end # module
