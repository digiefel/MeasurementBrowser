const CollectionIdPath = Tuple{Vararg{String}}

function _write_id_length(io::IO, length::Integer)::Nothing
    value = UInt64(length)
    for shift in 56:-8:0
        write(io, UInt8((value >> shift) & 0xff))
    end
    return nothing
end

function _write_id_text(io::IO, text::AbstractString)::Nothing
    bytes = codeunits(String(text))
    _write_id_length(io, length(bytes))
    write(io, bytes)
    return nothing
end

"""
Return a collection occurrence ID.

Derived from the parent occurrence ID, the concrete collection type, and `id(collection)`. Lengths
are written before each part so no combination of contents can collide with another.

This is a one-way digest: a collection is rebuilt from the stored `id(collection)` and its metadata,
never from this.
"""
function collection_record_id(
    parent_id::AbstractString,
    collection::AbstractCollection,
)::String
    io = IOBuffer()
    _write_id_text(io, "DataBrowser.collection.id.v1")
    _write_id_text(io, parent_id)
    _write_id_text(io, string(typeof(collection)))
    _write_id_text(io, id(collection))
    return bytes2hex(SHA.sha256(take!(io)))
end

"""Return final deterministic occurrence IDs for one concrete collection value path."""
function collection_id_path(
    path::Union{Tuple,AbstractVector{<:AbstractCollection}},
)::CollectionIdPath
    ids = String[]
    parent = ""
    for collection in path
        collection_id = collection_record_id(parent, collection)
        push!(ids, collection_id)
        parent = collection_id
    end
    return Tuple(ids)
end
