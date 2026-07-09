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

using ..DataBrowserAnnotations: AnnotationKey, collection_annotation_key

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
    return collection_annotation_key(inner)
end

function _read_all_sections(path::AbstractString)
    sections = Vector{Tuple{AnnotationKey,String}}()
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

        j = i + 1
        while j <= n && isempty(strip(lines[j]))
            j += 1
        end
        if j > n || strip(lines[j]) != FENCE
            throw(NotesParseError(path, i,
                "expected opening ``` fence after [$header]"))
        end

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
    read_section(root, key) -> String

Return the body for `key` from `notes.txt`, or `""` if absent.
"""
function read_section(root::AbstractString, key::AnnotationKey)::String
    sections = _read_all_sections(notes_path(root))
    target = collection_annotation_key(key)
    for (section_key, body) in sections
        section_key == target && return body
    end
    return ""
end

"""
    merged_view(root, key, ancestor_keys) -> Vector{NamedTuple}

Return ancestor sections in the order supplied (each marked
`editable=false`), followed by the section for `key` itself
(`editable=true`). Each entry has fields `(path, body, editable)`. Own
section is included even when its body is empty.
"""
function merged_view(root::AbstractString, key::AnnotationKey, ancestor_keys)
    sections = _read_all_sections(notes_path(root))
    by_path = Dict{AnnotationKey,String}()
    for (section_key, body) in sections
        by_path[section_key] = body
    end

    out = NamedTuple{(:path, :body, :editable),Tuple{AnnotationKey,String,Bool}}[]
    for ancestor in ancestor_keys
        ancestor_key = collection_annotation_key(ancestor)
        haskey(by_path, ancestor_key) || continue
        push!(out, (path=ancestor_key, body=by_path[ancestor_key], editable=false))
    end
    own_key = collection_annotation_key(key)
    own_body = get(by_path, own_key, "")
    push!(out, (path=own_key, body=own_body, editable=true))
    return out
end

"""
    write_section!(root, key, body)

Write `body` as the section for `key` in `notes.txt`. Other sections
are preserved in their original order; the target section is replaced
in place if present, otherwise appended.
"""
function write_section!(root::AbstractString, key::AnnotationKey, body::AbstractString)
    file = notes_path(root)
    sections = _read_all_sections(file)
    target = collection_annotation_key(key)
    new_body = String(body)

    found = false
    updated = Vector{Tuple{AnnotationKey,String}}(undef, length(sections))
    for (i, (section_key, existing)) in enumerate(sections)
        if section_key == target && !found
            updated[i] = (section_key, new_body)
            found = true
        else
            updated[i] = (section_key, existing)
        end
    end
    if !found
        push!(updated, (target, new_body))
    end

    open(file, "w") do io
        for (i, (section_key, sect_body)) in enumerate(updated)
            i > 1 && println(io)
            println(io, '[', section_key, ']')
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
