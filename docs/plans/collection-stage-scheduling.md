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
carry `kind::Symbol` — `nameof(typeof(value))` — in place of `registration_name`. Rebuilding goes
through a new method with no default:

```julia
reconstruct(::Type{T}, label, metadata) where {T<:AbstractCollection}
```

No default, unlike the item method. An item that cannot be rebuilt is recreated by rerunning
`read` → `entries`; a collection has no such path, so a missing method errors naming the type. It
fires from the first interpretation, since collection values are discarded at interpretation and
never retained.

The occurrence ID is a digest of the parent ID, the type, and `id(collection)`, so it cannot be
inverted. Label and own metadata are all a rebuild has. **A collection's identity must therefore be
reproducible from those two** — `test_collection_id_persistence.jl` is the case that proves it, with
two collections sharing a label and distinguished by an integer key that now has to live in
metadata.

**Both booleans are gone.** Collection process and analyze are scheduled without asking; the stage
defaults (`items` unchanged, empty `Dict`) make a project without stages a no-op, so "defined means
run" holds identically for typed and recipe projects. Nothing became eager — every scheduling site
is inside `if workspace.background_processing`, which defaults to `false`.

**`member_collection` is gone**, replaced by `collection_value(workspace, collections, key)`.

**`delivery_gate` always returns `ITEM_PROCESS`**, and `_delivered_payload` reads the deepest cached
stage: `:collection_processed` if a fold rewrote this member, else `:processed`. Delivery no longer
waits for a whole collection to fold before showing anything.

## One bug found on the way

`run_collection_process` paired outputs to inputs by `id(item)`. Typed items carry no id of their own
— identity is minted by the engine and lives on the record — so every member collided on `""` and the
identity fold falsely reported rewrites. Pairing is now positional, which the contract already
implied by requiring one output per input.

It was unreachable before only because typed collection stages never ran.

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
