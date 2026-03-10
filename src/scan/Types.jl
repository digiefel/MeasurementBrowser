using Dates

struct IndexedCsvFile
    id::String
    filepath::String
    filename::String
    timestamp::Union{DateTime,Nothing}
    header_summary::Dict{Symbol,Any}
end

struct MeasurementItem
    id::String
    source_file_id::String
    filepath::String
    kind::Symbol
    device_path::Vector{String}
    timestamp::Union{DateTime,Nothing}
    device_parameters::Dict{Symbol,Any}
    parameters::Dict{Symbol,Any}
    title::String
end

file_id(path::AbstractString) = abspath(String(path))

function item_id(
    source_file_id::AbstractString;
    cycle::Union{Nothing,Integer}=nothing,
    split::Union{Nothing,AbstractString}=nothing,
)
    base = String(source_file_id)
    if cycle !== nothing
        return base * "#cycle=$(Int(cycle))"
    elseif split !== nothing
        return base * "#split=$(String(split))"
    end
    return base
end
