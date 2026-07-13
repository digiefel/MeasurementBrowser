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

    plot_count = sum(length, values(project.plots); init=0)
    print(io, "\n  ", _plural(plot_count, "plot"))
    isempty(project.plots) || print(io, ":")
    for kind in sort!(collect(keys(project.plots)); by=String)
        for label in sort!(collect(keys(project.plots[kind])))
            print(io, "\n    ", kind, " → \"", label, "\"")
        end
    end
end

"""Create an empty project. Register items/stats/plots into it afterwards."""
function define_project(name::AbstractString; description::AbstractString="")::Project
    return Project(
        String(name),
        String(description),
        ItemRecipe[],
        Dict{Symbol,CollectionRecipe}(),
        Dict{Symbol,Dict{String,PlotRecipe}}(),
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

"""
Register (or replace) one plot recipe for an item kind.

`setup(workspace, items)` builds and returns the figure; `draw(workspace, items, figure)` fills it
in. `items` are the loaded, data-bearing items for the selection (`Vector{<:AbstractDataItem}`);
each has already passed through `process(item)`, if registered, so neither callback resolves item
data itself. A kind may have multiple plots; re-registering the same `label` for the same kind
replaces that plot, which keeps REPL iteration stable.
"""
function register_plot!(
    project::Project,
    kind::Symbol;
    label::AbstractString,
    setup::Function,
    draw::Function,
)::Project
    label_string = String(label)
    recipes = get!(project.plots, kind) do
        Dict{String,PlotRecipe}()
    end
    recipes[label_string] = PlotRecipe(kind, label_string, setup, draw)
    return project
end
