"""
Coords — spatial coordinates pulled from collection parameters.

Collection parameter rows may carry optional columns `x_um`, `y_um`, `w_um`, `h_um`. This module
reads them out into path-keyed dictionaries and computes axis-aligned bounding boxes from descendant
positions.
"""
module Coords

export read_positions, read_overrides, bounding_box, Rect

const CollectionParametersTable = Dict{Tuple{Vararg{String}},Dict{Symbol,Any}}

"""
    Rect

Axis-aligned rectangle in micrometres. `(x, y)` is the lower-left corner.
"""
struct Rect
    x::Float64
    y::Float64
    w::Float64
    h::Float64
end

function _path_key(segs::Tuple{Vararg{String}})
    return join(segs, '/')
end

function _as_float(value)
    value isa Real && return Float64(value)
    if value isa AbstractString
        s = strip(value)
        isempty(s) && return nothing
        parsed = tryparse(Float64, s)
        return parsed
    end
    return nothing
end

function _pair(params::Dict{Symbol,Any}, kx::Symbol, ky::Symbol)
    haskey(params, kx) || return nothing
    haskey(params, ky) || return nothing
    x = _as_float(params[kx])
    y = _as_float(params[ky])
    (x === nothing || y === nothing) && return nothing
    return (x, y)
end

"""
    read_positions(metadata) -> Dict{String, Tuple{Float64,Float64}}

Pull `(x_um, y_um)` for every row that supplies both columns. Keys are the
slash-joined path strings used elsewhere in the codebase.
"""
function read_positions(metadata::CollectionParametersTable)
    out = Dict{String,Tuple{Float64,Float64}}()
    for (segs, params) in metadata
        isempty(segs) && continue
        pair = _pair(params, :x_um, :y_um)
        pair === nothing && continue
        out[_path_key(segs)] = pair
    end
    return out
end

read_positions(::Nothing) = Dict{String,Tuple{Float64,Float64}}()

"""
    read_overrides(metadata) -> Dict{String, Tuple{Float64,Float64}}

Pull `(w_um, h_um)` for every row that supplies both columns. These act as
explicit bounding-box dimensions overriding any computed value.
"""
function read_overrides(metadata::CollectionParametersTable)
    out = Dict{String,Tuple{Float64,Float64}}()
    for (segs, params) in metadata
        isempty(segs) && continue
        pair = _pair(params, :w_um, :h_um)
        pair === nothing && continue
        out[_path_key(segs)] = pair
    end
    return out
end

read_overrides(::Nothing) = Dict{String,Tuple{Float64,Float64}}()

"""
    bounding_box(positions, descendant_paths; override=nothing) -> Union{Rect,Nothing}

Compute the axis-aligned rectangle covering the explicit positions of
`descendant_paths` looked up in `positions`. Paths missing from `positions`
are skipped. Returns `nothing` if no positions are available.

If `override` is supplied as `(w, h)`, it sets the rectangle's width and
height; the lower-left corner is still derived from the descendant minima.
A degenerate rectangle (zero descendants on one axis) gets a zero-width
or zero-height side; callers can pad as needed.
"""
function bounding_box(
    positions::Dict{String,Tuple{Float64,Float64}},
    descendant_paths;
    override::Union{Nothing,Tuple{Float64,Float64}}=nothing,
)
    xs = Float64[]
    ys = Float64[]
    for path in descendant_paths
        haskey(positions, path) || continue
        x, y = positions[path]
        push!(xs, x)
        push!(ys, y)
    end
    isempty(xs) && return nothing
    xmin = minimum(xs)
    ymin = minimum(ys)
    xmax = maximum(xs)
    ymax = maximum(ys)
    if override !== nothing
        w, h = override
        return Rect(xmin, ymin, w, h)
    end
    return Rect(xmin, ymin, xmax - xmin, ymax - ymin)
end

end # module
