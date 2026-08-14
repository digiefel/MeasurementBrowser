# Typed collections

Built. Kept as the record of why.

## The problem

Collection stages take the collection value, so a project can group the same items by device in one
level and by temperature in another and process each differently:

```julia
process(device::Device, items) -> Vector{<:AbstractDataItem}
analyze(temperature::Temperature, items) -> Dict
```

The engine could not produce that value. A `CollectionRecord` held key, id, parent, label, metadata,
and `registration_name` — a dialect-only string that let a `NamedCollection` be rebuilt and left
typed collections with nothing. So the engine reached the value by materializing a member and asking
it where it lived (`member_collection`), and answered "does this collection have a stage" with two
booleans the project declared:

```julia
_has_collection_process(::AbstractProject, ::Symbol)::Bool = false
```

Recipe projects derived them from their registration table. Typed projects got `false`, so a defined
collection stage silently never ran. That default was the one place the dialect had an engine
integration a typed project lacked.

## What was built

**Collections are rebuilt from their row.** `CollectionRecord` and the persisted `CollectionRow`
carry `kind::Symbol` — `nameof(typeof(value))` — in place of `registration_name`, plus `identity`:
`id(collection)` kept verbatim. Rebuilding goes through a new method with no default:

```julia
reconstruct(::Type{T}, id, metadata) where {T<:AbstractCollection}
```

No default, unlike the item method. An item that cannot be rebuilt is recreated by rerunning
`read` → `entries`; a collection has no such path, so a missing method errors naming the type. It
fires from the first interpretation, since collection values are discarded at interpretation and
never retained.

`id(::AbstractCollection)` is now required and must return a `String`. Identity is what `id` was
always for; it was previously only hashed into the occurrence digest and discarded, which is what
made rebuilding look impossible. Two dead ends came first and are recorded because they are tempting:
rebuilding from `label` (display text — two levels may legitimately share one) and from `metadata`
(forces identity into user-visible metadata; `NamedCollection` gained a `:name` field duplicating its
own label, which `test_hierarchy_edit` caught).

Narrowing `id` to `String` deletes ~150 lines of `collection_id.jl`: the canonical encoder existed to
hash arbitrary structs, dicts, and arrays into the digest, and served collections only. The digest is
now three length-prefixed strings.

**`NamedCollection` left `DataBrowserAPI`.** It had two consumers and belonged to neither layer above
them. `DataBrowserSources` now owns `DirectoryCollection` for directory levels; `DataBrowserRecipes`
owns `NamedCollection` for registered string paths. `annotate_collection_path` no longer type-checks
for a package type — it merges `metadata.txt` entries into any level through that level's own
`reconstruct`, so it now works for typed collections instead of silently skipping them.

**Both booleans are gone.** Collection process and analyze are scheduled without asking; the stage
defaults (`items` unchanged, empty `Dict`) make a project without stages a no-op, so "defined means
run" holds identically for typed and recipe projects. Nothing became eager — every scheduling site
is inside `if workspace.background_processing`, which defaults to `false`.

**`member_collection` is gone**, replaced by `collection_value(workspace, collections, key)`.

**`delivery_gate` always returns `ITEM_PROCESS`**, and `_delivered_payload` reads the deepest cached
stage: `:collection_processed` if a fold rewrote this member, else `:processed`. Delivery no longer
waits for a whole collection to fold before showing anything.

## What the id overload was hiding

`run_collection_process` paired outputs to inputs by `id(item)`. That worked for the dialect, whose
carrier answered `id` with the engine-built id, and broke for typed items, which answered with a
key the engine then wrapped — so members collided and the fold falsely reported rewrites. Pairing is
now positional, which the contract already implied by requiring one output per input.

Two more callers had the same defect: `select_items!` looked items up by `id(item)`, and
`item_annotation_key` keyed annotations by it, so typed items were unselectable by value and shared
one annotation key.

All three came from `id` naming two things — a user's key and the engine's built id — so the fix was
to stop having two. `id(value)::String` is now the identity for sources, source items, items, and
collections alike, used verbatim, with no default and no wrapping. `_mint_id` is deleted;
`DataBrowserRecipes` builds ids for its own items internally, so a `register_item!` user never sees
one. Item ids read as whatever a project returns — `cycle-2`, not `run.csv#cycles:cycle-2`.

## Cost

An identity fold per collection per open in eager mode: members are materialized and handed to a
function that returns them. Nothing is persisted — `run_collection_process` only writes outputs whose
data is not `===` its input. If it shows up in a benchmark, skip the fold when the previous run
rewrote nothing; the engine already computes that set.

## Follow-on

- Collection-level plots, dispatched on the collection type. This work is what makes them
  expressible: a plot for `Device` can now be resolved from a collection node with no member in hand.
- `request_processed_items` taking the requested stage as an argument, for the last-stage slider.
  Which stage is final belongs to the request, not the project.
- Whether `item_data`/`metadata`/`reconstruct` can be derived from a type declaration rather than
  hand-written, for both items and collections.
