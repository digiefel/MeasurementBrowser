"""One registered item pipeline."""
mutable struct ItemRecipe
    kind::Symbol
    detect::Function
    read::Function
    entries::Union{Nothing,Function}
    process::Union{Nothing,Function}
    analyze::Union{Nothing,Function}
    label::Union{Nothing,Function}
    collection::Union{Nothing,Function}
    id::Union{Nothing,Function}
end

"""
One registered collection recipe for a kind.

`process(data, metadata)` rewrites each member's data (one output per input).
`analyze(data, metadata)` folds the post-process members into a `Dict` attached to the collection.
"""
struct CollectionRecipe
    kind::Symbol
    process::Union{Nothing,Function}
    analyze::Union{Nothing,Function}
end

"""
A callback project assembled from registered recipes.

Source interpretation and data processing are defined by the registered callbacks.
Plot registration lives in `DataBrowserPlots`. Package-owned cache, job, and browser state does
not belong here.
"""
mutable struct Project
    name::String
    description::String
    recipes::Vector{ItemRecipe}
    collections::Dict{Symbol,CollectionRecipe}
end

const PROJECTS = Project[]
const DEFAULT_PROJECT = Ref{Union{Project,Nothing}}(nothing)
