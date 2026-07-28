"""
Premade recipes: complete registrations for formats common enough to ship.

Each one is written against `register_item!` and nothing else. A premade recipe must never reach
past the public dialect — if one needs something that is not there, that is a gap in the dialect to
close rather than a private hook to add.

Reading formats is what this package is for, so their parsers are ordinary dependencies and a
premade recipe simply works after `using DataBrowser`. A format heavy enough that projects not
using it should not pay for it belongs behind a package extension instead; that trades the
zero-ceremony story for load time, and only earns it when the dependency is genuinely large and not
already in the graph.
"""

"""Whether one source item's name carries any of the given lowercase extensions."""
function _has_extension(source_item::AbstractDataSourceItem, extensions::Vector{String})::Bool
    name = lowercase(label(source_item))
    return any(extension -> endswith(name, extension), extensions)
end

"""Parse one source item's file into a `DataFrame`."""
function _read_delimited(source_item::AbstractDataSourceItem, options::NamedTuple)::DataFrame
    path = source_item_path(source_item)
    path === nothing && error(
        "register_csv! needs a source item with a filesystem path, and " *
        "$(typeof(source_item)) has none; give that registration its own `read` callback",
    )
    return CSV.read(path, DataFrame; options...)
end

"""
    register_csv!(project, [kind]; extensions, read_options, callbacks...) -> project

Register delimited text files as items, one item per file.

`detect` matches by file extension (case-insensitively, `[".csv"]` by default) and `read` parses
with CSV.jl into a `DataFrame`. Every other keyword is an ordinary `register_item!` callback and is
forwarded untouched, so a premade recipe is a starting point rather than a different kind of thing:

```julia
register_csv!(project, :sweep;
    extensions=[".csv", ".txt"],
    read_options=(; delim='\\t', comment="#"),
    collection=(data, metadata) -> ["sweeps", metadata[:device]],
    analyze=(data, _metadata) -> Dict(:rows => nrow(data)),
)
```

Detection runs in registration order and the first match wins, so register narrower recipes before
this one when several could claim the same file. Because detection and reading go through the
generic source-item contract rather than any file type, this recipe works with any source whose
items expose a path.
"""
function register_csv!(
    project::Project,
    kind::Symbol;
    extensions::AbstractVector{<:AbstractString}=[".csv"],
    read_options::NamedTuple=NamedTuple(),
    kwargs...,
)::Project
    normalized = String[lowercase(extension) for extension in extensions]
    return register_item!(
        project,
        kind;
        detect=source_item -> _has_extension(source_item, normalized),
        read=source_item -> _read_delimited(source_item, read_options),
        kwargs...,
    )
end

register_csv!(project::Project; kwargs...)::Project = register_csv!(project, :csv; kwargs...)
