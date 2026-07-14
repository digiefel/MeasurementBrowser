"""Name of a callback, or "λ" for an anonymous function."""
_callback_name(f)::String = (n = string(nameof(f)); startswith(n, "#") ? "λ" : n)

"""Count with a pluralized noun, e.g. `1 plot`, `3 plots`."""
_plural(n::Integer, word::AbstractString)::String = "$n $word" * (n == 1 ? "" : "s")

Base.show(io::IO, project::Project) =
    print(io, "Project(\"$(project.name)\", ", _plural(length(project.recipes), "item kind"), ")")

function Base.show(io::IO, ::MIME"text/plain", project::Project)
    print(io, "Project \"$(project.name)\"")
    isempty(project.description) || print(io, " · ", project.description)

    print(io, "\n  ", _plural(length(project.recipes), "item kind"))
    isempty(project.recipes) || print(io, " (detection order):")
    pad = isempty(project.recipes) ? 0 : maximum(length(string(r.kind)) for r in project.recipes)
    for recipe in project.recipes
        optional = [name for (name, f) in
            (("entries", recipe.entries), ("process", recipe.process),
             ("analyze", recipe.analyze), ("label", recipe.label),
             ("collection", recipe.collection), ("id", recipe.id))
            if f !== nothing]
        steps = isempty(optional) ? "" : " · " * join(optional, " · ")
        print(io, "\n    ", rstrip(rpad(string(recipe.kind), pad) * steps))
    end

    print(io, "\n  ", _plural(length(project.collections), "collection recipe"))
    isempty(project.collections) || print(io, ":")
    for recipe in values(project.collections)
        stages = [name for (name, f) in
            (("process", recipe.process), ("analyze", recipe.analyze)) if f !== nothing]
        print(io, "\n    ", recipe.kind, " · ", join(stages, " · "))
    end
end

"""Create an empty project. Register item and collection recipes into it afterwards."""
function define_project(name::AbstractString; description::AbstractString="")::Project
    return Project(
        String(name),
        String(description),
        ItemRecipe[],
        Dict{Symbol,CollectionRecipe}(),
        Dict{String,SourceItemProfile}(),
        ReentrantLock(),
    )
end

const UNNAMED_ITEM_REGISTRATION = Symbol("#unnamed")

"""
    register_item!(project, [registration_name]; read, callbacks...) -> project

Register or replace one item pipeline. `read` is the only required callback. Registration
callbacks receive ordinary project data and metadata dictionaries; DataBrowser owns item identity
and storage records internally.

Calling the unnamed form again replaces the project's unnamed registration. Calling the named form
again with the same name replaces that registration in place.
"""
function register_item!(
    project::Project,
    kind::Symbol;
    read::Function,
    detect::Function=Returns(true),
    entries::Union{Nothing,Function}=nothing,
    process::Union{Nothing,Function}=nothing,
    analyze::Union{Nothing,Function}=nothing,
    label::Union{Nothing,Function}=nothing,
    collection::Union{Nothing,Function}=nothing,
    id::Union{Nothing,Function}=nothing,
)::Project
    recipe = ItemRecipe(
        kind, detect, read, entries, process, analyze, label, collection, id)
    index = findfirst(r -> r.kind === kind, project.recipes)
    index === nothing ? push!(project.recipes, recipe) : (project.recipes[index] = recipe)
    return project
end

function register_item!(project::Project; kwargs...)::Project
    return register_item!(project, UNNAMED_ITEM_REGISTRATION; kwargs...)
end

"""
Register (or replace) the collection recipe for one item kind.

`process(data, metadata)` receives aligned vectors and returns one data value per input.
`analyze(data, metadata)` receives the processed values and returns metadata for the collection.
Re-calling with the same registration name replaces the recipe.
"""
function register_collection_analysis!(
    project::Project,
    kind::Symbol;
    process::Union{Nothing,Function}=nothing,
    analyze::Union{Nothing,Function}=nothing,
)::Project
    project.collections[kind] = CollectionRecipe(kind, process, analyze)
    return project
end
