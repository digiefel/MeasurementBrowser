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

"""Timing retained for one source item in the current scan."""
mutable struct SourceItemProfile
    source_item_id::String
    source_item_label::String
    source_item_path::Union{Nothing,String}
    kind::Symbol
    item_count::Int
    detect_seconds::Float64
    read_seconds::Float64
    entries_seconds::Float64
    process_seconds::Float64
    analyze_seconds::Float64
    total_seconds::Float64
    thread_ids::Set{Int}
end

SourceItemProfile(source_item_id::AbstractString)::SourceItemProfile = SourceItemProfile(
    String(source_item_id), String(source_item_id), nothing,
    :unmatched, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, Set{Int}())

"""One item-kind summary derived from the current source-item timings."""
mutable struct KindProfileRow
    kind::Symbol
    source_items::Int
    items::Int
    detect_seconds::Float64
    read_seconds::Float64
    entries_seconds::Float64
    process_seconds::Float64
    analyze_seconds::Float64
    total_seconds::Float64
end

KindProfileRow(kind::Symbol)::KindProfileRow =
    KindProfileRow(kind, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

"""Read-only source-item row returned to the Performance window."""
struct SourceProfileRow
    source_item_id::String
    source_item_label::String
    source_item_path::Union{Nothing,String}
    kind::Symbol
    items::Int
    detect_seconds::Float64
    read_seconds::Float64
    entries_seconds::Float64
    process_seconds::Float64
    analyze_seconds::Float64
    total_seconds::Float64
    thread_ids::Vector{Int}
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
    # Transient timing for the latest scan, surfaced in the performance window. One row per source
    # item keeps the profile useful without retaining every expanded data-item event.
    scan_profile::Dict{String,SourceItemProfile}
    profile_lock::ReentrantLock
end

const PROJECTS = Project[]
const DEFAULT_PROJECT = Ref{Union{Project,Nothing}}(nothing)
