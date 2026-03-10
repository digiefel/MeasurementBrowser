using DataFrames

struct AnalysisResult
    key::Symbol
    label::String
    row_kind::Symbol
    table::DataFrame
    view_presets::Vector{NamedTuple}
    meta::Dict{Symbol,Any}
end
