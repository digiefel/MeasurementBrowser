# First-class Collections

## Purpose

Make collections concrete hierarchical values without making basic projects define Julia types.
This replaces the source-specific `collection_metadata(source, path)` mechanism with one collection
model shared by directory sources, registered projects, and type-based projects.

This is a focused foundation change. It does not add new collection views or workflows.

## Current model

Both public project styles currently produce collection paths:

```julia
collection = (data, metadata) -> ["Wafer A", "Chip A", "Device A"]
collection(item::MyItem) = ["Wafer A", "Chip A", "Device A"]
```

DataBrowser immediately converts every segment to `String`. The resulting `Vector{String}` is the
item's hierarchy address, the hierarchy index key, and the value persisted with the item record.
Hierarchy nodes retain the segment name and metadata, but not a public collection value.

`AbstractCollection` is exported but has no runtime role. `DirectorySource` instead answers
`collection_metadata(source, path)`, and the hierarchy merges the returned dictionaries into its
nodes. The cache stores source collection metadata separately so it can detect changes and restore
the hierarchy.

This model displays the hierarchy efficiently, but it cannot preserve a project-defined collection
type or provide one collection metadata interface across different sources.

## Settled behavior

### Collections are hierarchical

An item's collection is the complete path from the broadest group to the most specific. Different
branches may contain collections with the same label:

```text
Wafer A / Chip A / Device A
Wafer A / Chip B / Device A
```

The complete path distinguishes the two devices. A segment needs to be unique only among siblings;
its displayed label does not need to be globally unique.
The two "Device A" segments are both distinct, in the sense that they can have different metadata
and contain different data items, and are both alike, in the sense that it should be possible
to reorganize the hierarchy view to group their items under a common "Device A" parent.

### Basic projects use ordinary strings

The registration API continues to accept a path of strings:

```julia
collection = (data::Measurement, metadata::Dict) ->
    [metadata[:wafer], metadata[:chip], metadata[:device]]
```

DataBrowser supplies the corresponding collection values. Defining a collection type is not
required for ordinary directory-based projects.

When no `collection` callback is provided for a `DirectorySource`, the relative directory hierarchy
continues to provide the default path.

### Type-based projects can return concrete collections

A domain project may define lightweight collection types with typed fields:

```julia
struct Wafer <: AbstractCollection
    name::String
    batch::String
end

struct Chip <: AbstractCollection
    name::String
    column::Int
    row::Int
end

struct Device <: AbstractCollection
    name::String
    area_um2::Float64
end

collection(item::Measurement) = [item.wafer, item.chip, item.device]

metadata(wafer::Wafer) = Dict(:batch => wafer.batch)
metadata(chip::Chip) = Dict(:column => chip.column, :row => chip.row)
metadata(device::Device) = Dict(:area_um2 => device.area_um2)
```

The hierarchy retains these concrete values. Selecting or inspecting a collection can therefore
return the original `Wafer`, `Chip`, or `Device`, rather than a replacement string or package-defined
wrapper.

As with a custom data item, typed fields are the domain representation. `metadata(collection)` is
the queryable dictionary view exposed to DataBrowser.

### Reopening preserves collection values

A cached workspace reconstructs the collection hierarchy with the same concrete collection types
and field values. Reopening must not reduce a custom collection to its label or require all source
items to be interpreted before the collection is available.

The cache may use a private stable key and private reconstruction data. Those details do not change
the value returned to project code.

### Collection metadata has one public source

`metadata(collection)` supplies a collection's own metadata. Directory metadata files populate the
collection values created by `DirectorySource`; they do not use a second public metadata channel.

Metadata inherited by an item follows its collection path from parent to child. A more specific
collection wins when two levels provide the same key, and item-local metadata wins over inherited
collection metadata.

Collection `analyze` output remains distinct. It is an expensive pipeline result computed from the
collection's members, while `metadata(collection)` describes the collection before processing.

## Public contract to finalize before implementation

The implementation must settle these names and signatures in public documentation before changing
runtime code:

- how DataBrowser obtains the displayed label of a concrete `AbstractCollection`;
- how a concrete collection supplies stable sibling identity when its label may change;
- the exact accepted path container type while retaining the current vector syntax;
- how custom collection reconstruction is expressed when a value cannot use the default cache
  representation.

These are public API choices. The implementation must not infer them from hierarchy storage or cache
schema convenience.

## Internal requirements

### Path normalization

Registration strings and concrete `AbstractCollection` values enter one normalization boundary.
After that boundary, the hierarchy works with collection values plus a stable internal path key.
Display labels are not assumed to be globally unique keys.

The internal key represents the complete parent path, so duplicate leaf labels under different
parents remain distinct. Its physical representation may be strings, tuples, integer surrogates, or
another measured choice; it is not part of the public collection contract.

### Hierarchy

Each hierarchy node retains its concrete collection value. Tree rendering reads its display label,
metadata inheritance reads `metadata(collection)`, and collection selection can return the value.

Incremental insertion, removal, sorting, and copy-on-write publication must remain proportional to
the affected path. The change must not rebuild the complete hierarchy for each published item.

### DirectorySource

`DirectorySource` creates collection values for directory levels and for collection paths described
by its metadata file. File parsing remains source-specific, but the result enters the same
`AbstractCollection` hierarchy used by typed projects.

After migration, the public hooks `collection_metadata(source, path)` and
`has_collection_metadata(source)` are removed.

### Cache and reopening

The cache stores enough information to restore each hierarchy node's:

- stable identity within its parent;
- displayed label;
- concrete collection type and required field values;
- own metadata;
- parent relationship;
- collection-analysis result and result state.

The cache schema remains private. Source collection metadata no longer needs a parallel public
meaning, but expensive collection-analysis results remain independently invalidatable and cacheable.

Changes to a collection's identity, metadata, or parent invalidate the affected descendants and
their dependent item and collection work. Unchanged branches remain reusable.

### GUI and annotations

The existing tree remains the initial collection view. It renders labels from collection values and
keeps using stable internal keys for selection and annotations. This migration does not redesign the
tree, tags, notes, or spatial browser.

## Implementation sequence

1. Finalize the small public collection contract and add contract tests for string paths, concrete
   paths, duplicate labels, metadata, and reconstruction.
2. Normalize registration strings and typed collection values without converting custom values to
   strings.
3. Retain collection values in hierarchy nodes and migrate hierarchy indexing, edits, selection,
   and collection work scheduling to stable internal keys.
4. Make `DirectorySource` construct collection values and remove the old source/path metadata
   hooks.
5. Persist and reconstruct collection values and migrate collection invalidation and analysis
   storage.
6. Update the tree, annotations, public documentation, and examples to use the completed contract.
7. Remove obsolete cache fields, helpers, and terminology only after all readers use the new model.

Each step must leave the package compiling and its tests passing. Cache schema changes may replace
the pre-alpha cache cleanly; no compatibility shim is required.

## Required tests

- Two identical leaf labels under different parents remain separate collections.
- A registration returning strings receives the same visible hierarchy as today.
- A typed item returning mixed concrete collection types gets those exact values back.
- Collection metadata is inherited parent-to-child and overridden at the expected level.
- Directory metadata and typed collection metadata use the same runtime interface.
- A warm reopen restores concrete collection types, fields, hierarchy, and metadata before source
  reinterpretation finishes.
- Changing one collection invalidates only that branch and its dependent work.
- Tags and selections remain attached across reopen when collection identity is unchanged.
- Incremental publication does not rebuild or linearly scan the complete hierarchy per item.

## Non-goals

- collection querying or filtering UI;
- new tagging or note behavior;
- collection visualizers;
- graph or multi-parent membership;
- spatial navigation changes;
- workflow changes;
- renaming item metadata cache tables;
- changing item processing or item analysis callbacks.
