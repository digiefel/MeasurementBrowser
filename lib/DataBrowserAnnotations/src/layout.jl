"""
Layout — per-path explicit XY positions stored at `<root>/layout.txt`.

Used for nodes the user has dragged into place (typically unpositioned
container nodes). One line per path:

    <path>\\t<x_um>\\t<y_um>

Lines starting with `#` are comments. Blank lines are ignored.
"""
module Layout

using ..DataBrowserAnnotations: AnnotationKey, collection_annotation_key

export load, save, reset!, layout_path, LAYOUT_FILENAME

const LAYOUT_FILENAME = "layout.txt"

const PositionMap = Dict{AnnotationKey,Tuple{Float64,Float64}}

layout_path(root::AbstractString) = joinpath(String(root), LAYOUT_FILENAME)

struct LayoutParseError <: Exception
    path::String
    line_number::Int
    message::String
end

function Base.showerror(io::IO, err::LayoutParseError)
    print(io, "Invalid $(basename(err.path)) at line $(err.line_number): $(err.message)")
end

function load(root::AbstractString)::PositionMap
    path = layout_path(root)
    out = PositionMap()
    isfile(path) || return out
    for (line_number, raw) in enumerate(eachline(path))
        line = strip(raw)
        isempty(line) && continue
        startswith(line, '#') && continue
        parts = split(line, '\t')
        length(parts) == 3 ||
            throw(LayoutParseError(path, line_number, "expected '<path>\\t<x>\\t<y>'"))
        key = collection_annotation_key(strip(parts[1]))
        isempty(key) && throw(LayoutParseError(path, line_number, "empty path"))
        x = tryparse(Float64, strip(parts[2]))
        y = tryparse(Float64, strip(parts[3]))
        x === nothing && throw(LayoutParseError(path, line_number, "bad x value '$(parts[2])'"))
        y === nothing && throw(LayoutParseError(path, line_number, "bad y value '$(parts[3])'"))
        out[key] = (x, y)
    end
    return out
end

function save(root::AbstractString, positions::PositionMap)
    path = layout_path(root)
    if isempty(positions)
        rm(path; force=true)
        return
    end
    open(path, "w") do io
        for key in sort!(collect(keys(positions)))
            x, y = positions[key]
            println(io, key, '\t', x, '\t', y)
        end
    end
end

"""
    reset!(positions, paths; cols=nothing, spacing_um=200.0, origin=(0.0, 0.0))

Overwrite the entries for `paths` with a default grid layout. Grid is
`cols` columns wide; if `cols === nothing` it picks `ceil(sqrt(n))`.
Spacing is in micrometres, origin is the lower-left corner of the grid.
Entries already present for paths not in `paths` are left untouched.
"""
function reset!(
    positions::PositionMap,
    paths;
    cols::Union{Nothing,Integer}=nothing,
    spacing_um::Real=200.0,
    origin::Tuple{Real,Real}=(0.0, 0.0),
)
    pvec = collect(paths)
    n = length(pvec)
    n == 0 && return positions
    ncols = cols === nothing ? max(1, Int(ceil(sqrt(n)))) : Int(cols)
    ox, oy = Float64(origin[1]), Float64(origin[2])
    sp = Float64(spacing_um)
    for (i, path) in enumerate(pvec)
        col = (i - 1) % ncols
        row = (i - 1) ÷ ncols
        positions[collection_annotation_key(path)] = (ox + col * sp, oy + row * sp)
    end
    return positions
end

end # module
