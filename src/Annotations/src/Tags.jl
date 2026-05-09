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

Assignment rows are `<key>\\t<tag_name>`. Keys are either device-path strings
(slash-joined segments, e.g. `RuO2test/A9/VI/D1`) or measurement-ID strings
(absolute filesystem paths, optionally with `#cycle=N` / `#split=X` suffixes).
The two namespaces never overlap, so no kind prefix is needed. Fields are
tab-separated on write; whitespace-tolerant on read. Lines starting with `#`
are comments. Malformed rows raise `TagsParseError`.

`Tags.load` reads both `tags.txt` and `bad_measurements` when both are present.
Entries from `bad_measurements` are merged as `bad` assignments into the single
`assignments` map: `device <path>` and `measurement <id>` lines both land in
`assignments`, prefix-stripped. If the catalog has no `bad` entry it is added.
This makes repeated `load → save` cycles idempotent as long as `bad_measurements`
has not grown.

`Tags.save` writes `tags.txt` and simultaneously mirrors the `bad`-tagged subset
to `bad_measurements`. Device-path keys (no leading `/`) are written as
`device <key>`; measurement-ID keys (leading `/`) are written as
`measurement <key>`. If the state has no `bad`-tagged keys, `bad_measurements`
is removed.
"""
module Tags

export TagDef, TagState, load, save, effective, dominant_color,
       tags_path, TAGS_FILENAME

const TAGS_FILENAME = "tags.txt"
const LEGACY_BAD_FILENAME = "bad_measurements"
const LEGACY_BAD_TAG = "bad"
const LEGACY_BAD_COLOR = (UInt8(0xff), UInt8(0x30), UInt8(0x30))
const LEGACY_BAD_PRIORITY = 100

tags_path(root::AbstractString) = joinpath(String(root), TAGS_FILENAME)
_legacy_bad_path(root::AbstractString) = joinpath(String(root), LEGACY_BAD_FILENAME)

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

Keys are arbitrary strings — either device-path keys (slash-joined segments,
e.g. `"RuO2test/A9/VI/D1"`) or measurement-ID strings (absolute paths with
optional `#cycle=N` / `#split=X` suffixes). The two namespaces never overlap.

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

# Merge entries from bad_measurements into a mutable catalog + assignments pair.
function _merge_legacy!(catalog::Vector{TagDef},
                        assignments::Dict{String,Set{String}},
                        root::AbstractString)
    path = _legacy_bad_path(root)
    isfile(path) || return
    if !any(t -> t.name == LEGACY_BAD_TAG, catalog)
        push!(catalog, TagDef(LEGACY_BAD_TAG, LEGACY_BAD_COLOR, LEGACY_BAD_PRIORITY))
    end
    for raw in eachline(path)
        line = strip(raw)
        isempty(line) && continue
        startswith(line, '#') && continue
        parts = split(line; limit=2)
        length(parts) == 2 || continue
        kind = parts[1]
        value = strip(parts[2])
        isempty(value) && continue
        if kind == "device" || kind == "measurement"
            set = get!(() -> Set{String}(), assignments, String(value))
            push!(set, LEGACY_BAD_TAG)
        end
    end
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

    _merge_legacy!(catalog, assignments, root)

    return TagState(catalog, assignments)
end

function _format_color(color::NTuple{3,UInt8})
    return string(
        lpad(string(color[1]; base=16), 2, '0'),
        lpad(string(color[2]; base=16), 2, '0'),
        lpad(string(color[3]; base=16), 2, '0'),
    )
end

function _write_legacy_bad!(root::AbstractString, state::TagState)
    path = _legacy_bad_path(root)
    bad_keys = String[]
    for (key, tags) in state.assignments
        LEGACY_BAD_TAG in tags && push!(bad_keys, key)
    end
    if isempty(bad_keys)
        rm(path; force=true)
        return
    end
    sort!(bad_keys)
    open(path, "w") do io
        for key in bad_keys
            kind = startswith(key, '/') ? "measurement" : "device"
            println(io, kind, ' ', key)
        end
    end
end

function save(root::AbstractString, state::TagState)
    path = tags_path(root)
    if isempty(state.catalog) && isempty(state.assignments)
        rm(path; force=true)
        _write_legacy_bad!(root, state)
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
    _write_legacy_bad!(root, state)
end

"""
    effective(state, path, ancestor_paths) -> Set{String}

Return tags that apply to `path`: its own assignments unioned with any
assignments on ancestor paths. The caller supplies `ancestor_paths`
(order doesn't affect the union).

To get the full applicable set for a measurement, call:

    effective(state, measurement_id, [device_path; device_ancestors...])

This works because all keys — device paths and measurement IDs — live in the
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
