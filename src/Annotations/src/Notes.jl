"""
Notes — fenced per-path note sections stored at `<root>/notes.txt`.

Each section is `[<path>]` on its own line, then a triple-backtick
fence opening, then the body, then a closing triple-backtick fence:

    [ChipB]
    ```
    Oxygen flow: 5%
    ```
    [ChipB/SiteVI]
    ```
    dust particle visible
    [observed in second pass]
    ```

The parser keys on the **fences**, not the brackets, so users can
write `[anything]` freely inside bodies. A bracketed line is recognized
as a section header only when it's the bare `[<text>]` form and the
next non-blank line is exactly ```. Trailing newline before the
closing fence is trimmed; all other whitespace is preserved.
"""
module Notes

export read_section, merged_view, write_section!, notes_path,
       NOTES_FILENAME, NotesParseError

const NOTES_FILENAME = "notes.txt"
const FENCE = "```"

notes_path(root::AbstractString) = joinpath(String(root), NOTES_FILENAME)

struct NotesParseError <: Exception
    path::String
    line_number::Int
    message::String
end

function Base.showerror(io::IO, err::NotesParseError)
    print(io, "Invalid $(basename(err.path)) at line $(err.line_number): $(err.message)")
end

function _try_header(line::AbstractString)
    s = strip(line)
    (length(s) >= 2 && first(s) == '[' && last(s) == ']') || return nothing
    inner = strip(s[2:end-1])
    isempty(inner) && return nothing
    return String(inner)
end

function _read_all_sections(path::AbstractString)
    sections = Vector{Tuple{String,String}}()
    isfile(path) || return sections
    lines = readlines(path; keep=false)
    n = length(lines)
    i = 1
    while i <= n
        header = _try_header(lines[i])
        if header === nothing
            stripped = strip(lines[i])
            if !isempty(stripped)
                throw(NotesParseError(path, i,
                    "unexpected content outside of any section"))
            end
            i += 1
            continue
        end

        # Look ahead for the opening fence (skipping blanks).
        j = i + 1
        while j <= n && isempty(strip(lines[j]))
            j += 1
        end
        if j > n || strip(lines[j]) != FENCE
            # Not a real section header — but per spec, a header is only
            # recognized when followed by a fence. Treat as malformed.
            throw(NotesParseError(path, i,
                "expected opening ``` fence after [$header]"))
        end

        # Body lines until the closing fence.
        body_lines = String[]
        k = j + 1
        closed = false
        while k <= n
            if strip(lines[k]) == FENCE
                closed = true
                break
            end
            push!(body_lines, lines[k])
            k += 1
        end
        closed || throw(NotesParseError(path, j,
            "unterminated fence for section [$header]"))

        body = join(body_lines, '\n')
        push!(sections, (header, body))
        i = k + 1
    end
    return sections
end

"""
    read_section(root, path) -> String

Return the body for `path` from `notes.txt`, or `""` if absent.
"""
function read_section(root::AbstractString, path::AbstractString)::String
    sections = _read_all_sections(notes_path(root))
    target = String(path)
    for (key, body) in sections
        key == target && return body
    end
    return ""
end

"""
    merged_view(root, path, ancestor_paths) -> Vector{NamedTuple}

Return ancestor sections in the order supplied (each marked
`editable=false`), followed by the section for `path` itself
(`editable=true`). Each entry has fields `(path, body, editable)`. Own
section is included even when its body is empty.
"""
function merged_view(root::AbstractString, path::AbstractString, ancestor_paths)
    sections = _read_all_sections(notes_path(root))
    by_path = Dict{String,String}()
    for (key, body) in sections
        by_path[key] = body
    end

    out = NamedTuple{(:path, :body, :editable),Tuple{String,String,Bool}}[]
    for ancestor in ancestor_paths
        key = String(ancestor)
        haskey(by_path, key) || continue
        push!(out, (path=key, body=by_path[key], editable=false))
    end
    own_key = String(path)
    own_body = get(by_path, own_key, "")
    push!(out, (path=own_key, body=own_body, editable=true))
    return out
end

"""
    write_section!(root, path, body)

Write `body` as the section for `path` in `notes.txt`. Other sections
are preserved in their original order; the target section is replaced
in place if present, otherwise appended.
"""
function write_section!(root::AbstractString, path::AbstractString, body::AbstractString)
    file = notes_path(root)
    sections = _read_all_sections(file)
    target = String(path)
    new_body = String(body)

    found = false
    updated = Vector{Tuple{String,String}}(undef, length(sections))
    for (i, (key, existing)) in enumerate(sections)
        if key == target && !found
            updated[i] = (key, new_body)
            found = true
        else
            updated[i] = (key, existing)
        end
    end
    if !found
        push!(updated, (target, new_body))
    end

    open(file, "w") do io
        for (i, (key, sect_body)) in enumerate(updated)
            i > 1 && println(io)
            println(io, '[', key, ']')
            println(io, FENCE)
            if !isempty(sect_body)
                print(io, sect_body)
                endswith(sect_body, '\n') || println(io)
            end
            println(io, FENCE)
        end
    end
end

end # module
