"""
Deprecated shared CSV-reader package.

This package is retained only for existing generic readers that have not yet moved to project-owned
source loading. Do not add new project readers here. New source-file interpretation should live with
the project that understands the file format.
"""
module DataLoader

using CSV
using DataFrames
using Dates

export find_files, get_file_patterns, read_iv_sweep, read_tlm_4p

"""
Read one voltage-driven I-V sweep and return its measured current and voltage as `:i` and `:v`.
"""
function read_iv_sweep(
    filename::AbstractString,
    workdir::AbstractString=".",
)::DataFrame
    filepath = joinpath(workdir, filename)

    # 1. Find the header line index
    lines = readlines(filepath)
    header_line = 1
    for (i, line) in enumerate(lines)
        # A header line typically contains "Voltage" or "V" AND "Current" or "I"
        # and is NOT a data line (doesn't start with a number)
        if (occursin("Voltage", line) || occursin("V", line)) &&
           (occursin("Current", line) || occursin("I", line)) &&
           !occursin(r"^-?\d", line) &&
           !occursin(' ', line)
            header_line = i
            break
        end
    end

    # 2. Use CSV.read
    # silence warnings about empty lines or metadata
    df = CSV.read(filepath, DataFrame; header=header_line, silencewarnings=true)

    # 3. Identify columns robustly
    cols = names(df)

    # Helper to find column matching candidates
    function find_col(candidates)
        for cand in candidates
            # Exact match or starts with
            match = findfirst(c -> c == cand || startswith(c, cand), cols)
            if match !== nothing
                return cols[match]
            end
        end
        # Fallback: contains
        for cand in candidates
            match = findfirst(c -> occursin(cand, c), cols)
            if match !== nothing
                return cols[match]
            end
        end
        return nothing
    end

    v_col = find_col(["Voltage", "V", "Voltage_V", "VoltageHigh_V"])
    i_col = find_col(["Current", "I", "Current_A", "Current_High_A", "I1"])

    if v_col === nothing || i_col === nothing
        # Last ditch effort: look for "V" and "I" anywhere, excluding common non-data words
        if v_col === nothing
            idx = findfirst(c -> occursin("V", c) && !occursin("Time", c), cols)
            v_col = idx !== nothing ? cols[idx] : nothing
        end
        if i_col === nothing
            idx = findfirst(c -> occursin("I", c) && !occursin("Time", c) && !occursin("Info", c) && !occursin("Index", c), cols)
            i_col = idx !== nothing ? cols[idx] : nothing
        end

        (v_col === nothing || i_col === nothing) &&
            error("Could not identify current and voltage columns in $filepath. Columns: $cols")
    end

    return DataFrame(i=df[!, i_col], v=df[!, v_col])
end


"""
Read one four-terminal I-V measurement and return source current and measured voltage as `:i` and
`:v`.
"""
function read_tlm_4p(
    filename::AbstractString,
    workdir::AbstractString=".",
)::DataFrame
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
    i = Float64[]
    v = Float64[]

    if is_new_format
        for line in data_lines
            if !isempty(strip(line))
                parts = split(line, ',')
                try
                    curr = parse(Float64, parts[1])
                    v_high = parse(Float64, parts[2])
                    v_low = parse(Float64, parts[3])
                    push!(i, curr)
                    push!(v, v_high - v_low) # Use measured drop (high - low)
                catch
                    continue
                end
            end
        end
    else
        for line in data_lines
            if !isempty(strip(line))
                parts = split(line, ',')
                # Filter out empty parts
                valid_parts = filter(p -> !isempty(strip(p)), parts)
                try
                    curr = parse(Float64, valid_parts[1])
                    dv = parse(Float64, valid_parts[2])
                    push!(i, curr)
                    push!(v, dv)
                catch e
                    println("Error parsing line: $line, error: $e")
                    continue
                end
            end
        end
    end

    return DataFrame(i=i, v=v)
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
        cv_sweep=r"CVSweep",
    )
end


end # module DataLoader
