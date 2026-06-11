using GLMakie
using DataFrames
using Statistics

using DataAnalysis: analyze_tlm_combined, calculate_sheet_resistance, extract_tlm_geometry_from_params

function _plot_title(path::AbstractString)
    return strip(replace(basename(String(path)), r"\.csv$" => ""))
end

function _format_frequency_label(freq_Hz::Real)
    freq_Hz >= 1e6 && return "$(round(freq_Hz / 1e6; digits=3)) MHz"
    freq_Hz >= 1e3 && return "$(round(freq_Hz / 1e3; digits=3)) kHz"
    return "$(round(freq_Hz; digits=3)) Hz"
end

function _extract_temperature_K(params::Dict{Symbol,Any}, filepath::String)
    if haskey(params, :temperature_K)
        return try
            Float64(params[:temperature_K])
        catch
            NaN
        end
    end

    if (m = match(r"(\d+(?:\.\d+)?)K", basename(filepath))) !== nothing
        return try
            parse(Float64, m.captures[1])
        catch
            NaN
        end
    end

    return NaN
end

function _extract_site_label(params::Dict{Symbol,Any}, filepath::String)
    for key in (:site, :chip, :device_id, :die, :location)
        if haskey(params, key)
            return string(params[key])
        end
    end

    name_part = replace(basename(filepath), r"\.(csv|txt)$" => "")
    tokens = split(name_part, '_')
    device_tokens = String[]
    for token in tokens
        if occursin(r"^\d{6,}$", token) ||
           occursin(r"^\d{4}-\d{2}-\d{2}$", token) ||
           occursin(r"^\d+K$", token) ||
           occursin(r"^L\d+W\d+$", token)
            break
        end
        push!(device_tokens, token)
    end
    return isempty(device_tokens) ? name_part : join(device_tokens, "_")
end

function _extract_oxygen_percent(params::Dict{Symbol,Any}, filepath::String)
    for key in (:oxygen_percent, :o2_percent)
        if haskey(params, key)
            return try
                Float64(params[key])
            catch
                NaN
            end
        end
    end

    if (m = match(r"(\d+(?:\.\d+)?)%O2", basename(filepath))) !== nothing
        return try
            parse(Float64, m.captures[1])
        catch
            NaN
        end
    end

    return NaN
end

function _extract_oxygen_flow(params::Dict{Symbol,Any}, filepath::String)
    for key in (:oxygen_flow_sccm, :o2_flow_sccm, :oxygen_flow, :o2_flow)
        if haskey(params, key)
            return try
                Float64(params[key])
            catch
                NaN
            end
        end
    end

    if (m = match(r"(\d+(?:\.\d+)?)\s*sccm", basename(filepath))) !== nothing
        return try
            parse(Float64, m.captures[1])
        catch
            NaN
        end
    end

    return NaN
end

function _ols_ab(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})
    n = length(x)
    sx = sum(x)
    sy = sum(y)
    sxx = sum(x .^ 2)
    sxy = sum(x .* y)
    denom = n * sxx - sx^2
    if !isfinite(denom) || abs(denom) < 1e-12
        return (NaN, NaN, false)
    end
    b = (n * sxy - sx * sy) / denom
    a = (sy - b * sx) / n
    return (a, b, true)
end
