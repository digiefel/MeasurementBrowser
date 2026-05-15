export read_fe_pund, read_pund_wakeup_amplitude, read_pund_wakeup_reps

"""
Read FE PUND data from CSV file
"""
function read_fe_pund(filename, workdir=".")
    filepath = joinpath(workdir, filename)
    lines = readlines(filepath)
    data_start = 0
    column_map = nothing

    for (i, line) in enumerate(lines)
        columns = strip.(split(line, ','))
        if "Time" in columns &&
           "MeasResult1_value" in columns &&
           "MeasResult2_value" in columns
            indices = Dict(column => index for (index, column) in pairs(columns))
            column_map = (
                time=indices["Time"],
                current=indices["MeasResult1_value"],
                voltage=indices["MeasResult2_value"],
                current_time=get(indices, "MeasResult1_time", indices["Time"]),
                voltage_time=get(indices, "MeasResult2_time", indices["Time"]),
            )
            data_start = i + 1
            break
        elseif "Time_s" in columns &&
               "Current_A" in columns &&
               "Voltage_V" in columns
            indices = Dict(column => index for (index, column) in pairs(columns))
            column_map = (
                time=indices["Time_s"],
                current=indices["Current_A"],
                voltage=indices["Voltage_V"],
                current_time=indices["Time_s"],
                voltage_time=indices["Time_s"],
            )
            data_start = i + 1
            break
        end
    end

    if data_start == 0 || column_map === nothing
        error("Could not find FE PUND data header in $filepath")
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
            if length(parts) >= maximum(column_map)
                try
                    push!(time, parse(Float64, parts[column_map.time]))
                    push!(current, parse(Float64, parts[column_map.current]))
                    push!(voltage, parse(Float64, parts[column_map.voltage]))
                    push!(current_time, parse(Float64, parts[column_map.current_time]))
                    push!(voltage_time, parse(Float64, parts[column_map.voltage_time]))
                catch
                    continue
                end
            end
        end
    end

    isempty(time) && error("No numeric FE PUND rows could be parsed from $filepath")

    return DataFrame(time=time, current=current, voltage=voltage,
        current_time=current_time, voltage_time=voltage_time)
end

"""
Read a single amplitude's data from a PUND_Wakeup CSV file.
Columns: Vmax_V, Time_s, Voltage_V, Current_A (header rows prefixed with #).
Returns DataFrame with :time, :current, :voltage matching analyze_pund_and_pn expectations.
rep selects which repetition group to return (1-based). Groups are separated by time gaps
larger than REP_GAP_S (default 1 ms), which distinguishes consecutive read events separated
by fatigue cycling from individual sample points within one read.
"""
const _WAKEUP_REP_GAP_S = 1e-3

function read_pund_wakeup_amplitude(filename, workdir, amplitude_V::Float64, rep::Int=1)
    filepath = joinpath(workdir, filename)
    groups = Vector{NTuple{3,Float64}}[]
    current_group = NTuple{3,Float64}[]
    last_t = -Inf
    open(filepath, "r") do io
        for line in eachline(io)
            startswith(line, '#') && continue
            isempty(strip(line)) && continue
            parts = split(line, ',')
            length(parts) >= 4 || continue
            v = tryparse(Float64, parts[1])
            (v === nothing || v != amplitude_V) && continue
            t   = tryparse(Float64, parts[2])
            vol = tryparse(Float64, parts[3])
            c   = tryparse(Float64, parts[4])
            (t === nothing || vol === nothing || c === nothing) && continue
            if last_t != -Inf && t - last_t > _WAKEUP_REP_GAP_S
                push!(groups, current_group)
                current_group = NTuple{3,Float64}[]
            end
            push!(current_group, (t, vol, c))
            last_t = t
        end
    end
    isempty(current_group) || push!(groups, current_group)
    isempty(groups) && return DataFrame(time=Float64[], current=Float64[], voltage=Float64[])
    grp = groups[clamp(rep, 1, length(groups))]
    return DataFrame(time=[x[1] for x in grp], voltage=[x[2] for x in grp], current=[x[3] for x in grp])
end

"""
Count repetition groups per amplitude in a PUND_Wakeup CSV.
Returns Dict{Float64,Int} mapping each amplitude to its rep count.
"""
function read_pund_wakeup_reps(filename, workdir=".")
    filepath = joinpath(workdir, filename)
    amp_last_t = Dict{Float64,Float64}()
    amp_reps   = Dict{Float64,Int}()
    open(filepath, "r") do io
        for line in eachline(io)
            startswith(line, '#') && continue
            isempty(strip(line)) && continue
            parts = split(line, ',')
            length(parts) >= 4 || continue
            v = tryparse(Float64, parts[1])
            v === nothing && continue
            t = tryparse(Float64, parts[2])
            t === nothing && continue
            last_t = get(amp_last_t, v, -Inf)
            if last_t == -Inf
                amp_reps[v] = 1
            elseif t - last_t > _WAKEUP_REP_GAP_S
                amp_reps[v] = get(amp_reps, v, 1) + 1
            end
            amp_last_t[v] = t
        end
    end
    return amp_reps
end
