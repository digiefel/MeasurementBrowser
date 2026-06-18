"""
Tags — catalog of named tags plus per-key assignments stored at `<root>/tags.txt`.

The file has two sections, each introduced by a bracketed header:

    [catalog]
    bad\\tff3030\\t100
    todo\\t30c0ff\\t50

    [assignments]
    RuO2test/A9/VI/D1                                          bad
    /Users/davide/data/RuO2test/A9/VI/D1/3V FE PUND.csv        bad
    RuO2test/A10/VI                                            todo

Assignment rows are `<key>\\t<tag_name>`. Keys are either collection-path strings
(slash-joined segments, e.g. `RuO2test/A9/VI/D1`) or item keys.
The two namespaces never overlap, so no kind prefix is needed. Fields are
tab-separated on write; whitespace-tolerant on read. Lines starting with `#`
are comments. Malformed rows raise `TagsParseError`.

`Tags.load` and `Tags.save` only use `tags.txt`.
"""
module Tags

export TagDef, TagState, load, save, effective, dominant_color,
       tags_path, TAGS_FILENAME

const TAGS_FILENAME = "tags.txt"

tags_path(root::AbstractString) = joinpath(String(root), TAGS_FILENAME)

"""
    TagDef(name, color, priority)

A single catalog entry. `color` is a `(r, g, b)` tuple of bytes;
`priority` is an integer (higher wins when multiple tags apply).
"""
struct TagDef
    name::String
    color::NTuple{3,UInt8}
    priority::Int
end

"""
    TagState

Catalog of `TagDef`s plus a single `assignments` map.

Keys are arbitrary strings — either collection-path keys (slash-joined segments,
e.g. `"RuO2test/A9/VI/D1"`) or item keys. The two namespaces never overlap.

Only explicitly attached tags are stored; inheritance is applied at lookup time
via `effective`.
"""
struct TagState
    catalog::Vector{TagDef}
    assignments::Dict{String,Set{String}}
end

TagState() = TagState(TagDef[], Dict{String,Set{String}}())

struct TagsParseError <: Exception
    path::String
    line_number::Int
    message::String
end

function Base.showerror(io::IO, err::TagsParseError)
    print(io, "Invalid $(basename(err.path)) at line $(err.line_number): $(err.message)")
end

function _parse_color_hex(path::String, line_number::Int, raw::AbstractString)
    s = strip(raw)
    startswith(s, '#') && (s = s[2:end])
    length(s) == 6 || throw(TagsParseError(path, line_number,
        "expected 6-digit hex color, got '$raw'"))
    r = tryparse(UInt8, s[1:2]; base=16)
    g = tryparse(UInt8, s[3:4]; base=16)
    b = tryparse(UInt8, s[5:6]; base=16)
    (r === nothing || g === nothing || b === nothing) &&
        throw(TagsParseError(path, line_number, "bad hex color '$raw'"))
    return (r, g, b)
end

function _split_row(line::AbstractString)
    if occursin('\t', line)
        return [strip(p) for p in split(line, '\t') if !isempty(strip(p))]
    end
    return [strip(p) for p in split(line) if !isempty(strip(p))]
end

function load(root::AbstractString)::TagState
    path = tags_path(root)

    catalog = TagDef[]
    assignments = Dict{String,Set{String}}()

    if isfile(path)
        section = :none
        for (line_number, raw) in enumerate(eachline(path))
            body = strip(raw)
            isempty(body) && continue
            startswith(body, '#') && continue
            if startswith(body, '[') && endswith(body, ']')
                header = strip(body[2:end-1])
                if header == "catalog"
                    section = :catalog
                elseif header == "assignments"
                    section = :assignments
                else
                    throw(TagsParseError(path, line_number, "unknown section '[$header]'"))
                end
                continue
            end

            fields = _split_row(body)
            if section == :catalog
                length(fields) == 3 || throw(TagsParseError(path, line_number,
                    "expected '<name>\\t<color>\\t<priority>'"))
                name = String(fields[1])
                isempty(name) && throw(TagsParseError(path, line_number, "empty tag name"))
                color = _parse_color_hex(path, line_number, fields[2])
                priority = tryparse(Int, fields[3])
                priority === nothing && throw(TagsParseError(path, line_number,
                    "bad priority '$(fields[3])'"))
                push!(catalog, TagDef(name, color, priority))
            elseif section == :assignments
                length(fields) == 2 || throw(TagsParseError(path, line_number,
                    "expected '<key>\\t<tag_name>'"))
                key = String(fields[1])
                tag = String(fields[2])
                (isempty(key) || isempty(tag)) &&
                    throw(TagsParseError(path, line_number, "empty key or tag name"))
                set = get!(() -> Set{String}(), assignments, key)
                push!(set, tag)
            else
                throw(TagsParseError(path, line_number,
                    "row outside of any section; expected '[catalog]' or '[assignments]' header first"))
            end
        end
    end

    return TagState(catalog, assignments)
end

function _format_color(color::NTuple{3,UInt8})
    return string(
        lpad(string(color[1]; base=16), 2, '0'),
        lpad(string(color[2]; base=16), 2, '0'),
        lpad(string(color[3]; base=16), 2, '0'),
    )
end

function save(root::AbstractString, state::TagState)
    path = tags_path(root)
    if isempty(state.catalog) && isempty(state.assignments)
        rm(path; force=true)
        return
    end

    open(path, "w") do io
        println(io, "[catalog]")
        for tag in state.catalog
            println(io, tag.name, '\t', _format_color(tag.color), '\t', tag.priority)
        end

        if !isempty(state.assignments)
            println(io)
            println(io, "[assignments]")
            for key in sort!(collect(keys(state.assignments)))
                tags = sort!(collect(state.assignments[key]))
                for tag in tags
                    println(io, key, '\t', tag)
                end
            end
        end
    end
end

"""
    effective(state, path, ancestor_paths) -> Set{String}

Return tags that apply to `path`: its own assignments unioned with any
assignments on ancestor paths. The caller supplies `ancestor_paths`
(order doesn't affect the union).

To get the full applicable set for an item, call:

    effective(state, item_key, [collection_path; collection_ancestors...])

This works because all keys — collection paths and item keys — live in the
same `assignments` map and are looked up uniformly.
"""
function effective(state::TagState, path::AbstractString, ancestor_paths)::Set{String}
    out = Set{String}()
    for ancestor in ancestor_paths
        tags = get(state.assignments, String(ancestor), nothing)
        tags === nothing && continue
        union!(out, tags)
    end
    own = get(state.assignments, String(path), nothing)
    own !== nothing && union!(out, own)
    return out
end

"""
    dominant_color(state, effective_tags) -> Union{Nothing, NTuple{3,UInt8}}

Resolve `effective_tags` against the catalog and return the highest-
priority hit's color. Returns `nothing` if `effective_tags` is empty or
contains no catalog matches.
"""
function dominant_color(state::TagState, effective_tags)::Union{Nothing,NTuple{3,UInt8}}
    isempty(effective_tags) && return nothing
    best::Union{Nothing,TagDef} = nothing
    for tag in state.catalog
        in(tag.name, effective_tags) || continue
        if best === nothing || tag.priority > best.priority
            best = tag
        end
    end
    best === nothing && return nothing
    return best.color
end

end # module
