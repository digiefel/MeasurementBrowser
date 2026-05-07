"""
Tags — catalog of named tags plus per-path and per-measurement-ID assignments stored at
`<root>/tags.txt`.

The file has two sections, each introduced by a bracketed header:

    [catalog]
    bad\\tff3030\\t100
    todo\\t30c0ff\\t50

    [assignments]
    device         RuO2test/A9/VI/D1    bad
    measurement    abc123hash           bad
    device         RuO2test/A10/VI      todo

Assignment rows carry an explicit kind prefix (`device` or `measurement`), followed by the key
and the tag name. Fields are tab-separated on write; whitespace-tolerant on read. Lines starting
with `#` are comments. Unknown kind tokens raise `TagsParseError`.

`Tags.load` reads both `tags.txt` and `bad_measurements` when both are present. Entries from
`bad_measurements` are merged as `bad` assignments: `device <path>` lines into `assignments`,
`measurement <id>` lines into `measurement_assignments`. If the catalog has no `bad` entry it
is added. This makes repeated `load → save` cycles idempotent as long as `bad_measurements`
has not grown.

`Tags.save` writes only `tags.txt`. After `load → save`, subsequent `load → save` cycles
produce byte-identical output as long as `bad_measurements` has not grown.
"""
module Tags

export TagDef, TagState, load, save, effective, dominant_color, assigned_to_measurement,
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

Catalog of `TagDef`s plus two assignment maps:

- `assignments` — device-path-keyed: explicit tags per path string.
- `measurement_assignments` — measurement-ID-keyed: explicit tags per ID.

Both maps record only explicitly attached tags; inheritance is applied at lookup time via
`effective` (device paths) or `assigned_to_measurement` (measurement IDs). Measurement-ID
lookup is not ancestor-walked; callers compose `effective` and `assigned_to_measurement`
themselves.
"""
struct TagState
    catalog::Vector{TagDef}
    assignments::Dict{String,Set{String}}
    measurement_assignments::Dict{String,Set{String}}
end

TagState() = TagState(TagDef[], Dict{String,Set{String}}(), Dict{String,Set{String}}())

struct TagsParseError <: Exception
    path::String
    line_number::Int
    message::String
end

function Base.showerror(io::IO, err::TagsParseError)
    print(io, "Invalid $(basename(err.path)) at line $(err.line_number): $(err.message)")
end

function _strip_comment(line::AbstractString)
    idx = findfirst('#', line)
    idx === nothing && return line
    return line[1:prevind(line, idx)]
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

# Merge entries from bad_measurements into a mutable catalog + assignment pair.
function _merge_legacy!(catalog::Vector{TagDef},
                        assignments::Dict{String,Set{String}},
                        measurement_assignments::Dict{String,Set{String}},
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
        if kind == "device"
            set = get!(() -> Set{String}(), assignments, String(value))
            push!(set, LEGACY_BAD_TAG)
        elseif kind == "measurement"
            set = get!(() -> Set{String}(), measurement_assignments, String(value))
            push!(set, LEGACY_BAD_TAG)
        end
    end
end

function load(root::AbstractString)::TagState
    path = tags_path(root)

    catalog = TagDef[]
    assignments = Dict{String,Set{String}}()
    measurement_assignments = Dict{String,Set{String}}()

    if isfile(path)
        section = :none
        for (line_number, raw) in enumerate(eachline(path))
            body = strip(_strip_comment(raw))
            isempty(body) && continue
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
                length(fields) == 3 || throw(TagsParseError(path, line_number,
                    "expected '<kind>\\t<key>\\t<tag_name>'"))
                kind_tok = String(fields[1])
                key = String(fields[2])
                tag = String(fields[3])
                (isempty(key) || isempty(tag)) &&
                    throw(TagsParseError(path, line_number, "empty key or tag name"))
                if kind_tok == "device"
                    set = get!(() -> Set{String}(), assignments, key)
                    push!(set, tag)
                elseif kind_tok == "measurement"
                    set = get!(() -> Set{String}(), measurement_assignments, key)
                    push!(set, tag)
                else
                    throw(TagsParseError(path, line_number,
                        "unknown assignment kind '$kind_tok'; expected 'device' or 'measurement'"))
                end
            else
                throw(TagsParseError(path, line_number,
                    "row outside of any section; expected '[catalog]' or '[assignments]' header first"))
            end
        end
    end

    _merge_legacy!(catalog, assignments, measurement_assignments, root)

    return TagState(catalog, assignments, measurement_assignments)
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
    if isempty(state.catalog) && isempty(state.assignments) && isempty(state.measurement_assignments)
        rm(path; force=true)
        return
    end

    open(path, "w") do io
        println(io, "[catalog]")
        for tag in state.catalog
            println(io, tag.name, '\t', _format_color(tag.color), '\t', tag.priority)
        end

        has_assignments = !isempty(state.assignments) || !isempty(state.measurement_assignments)
        if has_assignments
            println(io)
            println(io, "[assignments]")
            for key in sort!(collect(keys(state.assignments)))
                tags = sort!(collect(state.assignments[key]))
                for tag in tags
                    println(io, "device\t", key, '\t', tag)
                end
            end
            for key in sort!(collect(keys(state.measurement_assignments)))
                tags = sort!(collect(state.measurement_assignments[key]))
                for tag in tags
                    println(io, "measurement\t", key, '\t', tag)
                end
            end
        end
    end
end

"""
    effective(state, path, ancestor_paths) -> Set{String}

Return tags that apply to `path`: its own device-path assignments unioned with any
assignments on ancestor paths. The caller supplies `ancestor_paths`
(order doesn't affect the union).

Measurement-ID tags are not included here; callers compose `effective` and
`assigned_to_measurement` themselves.
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
    assigned_to_measurement(state, measurement_id) -> Set{String}

Return the tags explicitly assigned to `measurement_id`. Returns an empty set when
no tags are recorded for that ID.

Measurement-ID tags are not ancestor-walked. Callers combine this with `effective`
on the device path to get the full applicable set.
"""
function assigned_to_measurement(state::TagState, measurement_id::AbstractString)::Set{String}
    tags = get(state.measurement_assignments, String(measurement_id), nothing)
    tags === nothing && return Set{String}()
    return copy(tags)
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
