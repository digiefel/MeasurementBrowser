const CollectionIdPath = Tuple{Vararg{String}}

function _write_id_length(io::IO, length::Integer)::Nothing
    value = UInt64(length)
    for shift in 56:-8:0
        write(io, UInt8((value >> shift) & 0xff))
    end
    return nothing
end

function _write_id_bytes(io::IO, bytes::AbstractVector{UInt8})::Nothing
    _write_id_length(io, length(bytes))
    write(io, bytes)
    return nothing
end

function _write_id_text(io::IO, text::AbstractString)::Nothing
    _write_id_bytes(io, codeunits(String(text)))
    return nothing
end

function _canonical_id_bytes(value)::Vector{UInt8}
    io = IOBuffer()
    _write_canonical_id!(io, value, IdDict{Any,Nothing}())
    return take!(io)
end

function _write_nested_id!(io::IO, value, seen::IdDict{Any,Nothing})::Nothing
    nested = IOBuffer()
    _write_canonical_id!(nested, value, seen)
    _write_id_bytes(io, take!(nested))
    return nothing
end

function _write_canonical_id!(io::IO, value, seen::IdDict{Any,Nothing})::Nothing
    type = typeof(value)
    write(io, UInt8('T'))
    _write_id_text(io, string(type))

    if value === nothing
        write(io, UInt8('n'))
    elseif value === missing
        write(io, UInt8('m'))
    elseif value isa Bool
        write(io, UInt8('b'), UInt8(value))
    elseif value isa Integer
        write(io, UInt8('i'))
        _write_id_text(io, string(value))
    elseif value isa BigFloat
        write(io, UInt8('F'))
        _write_id_text(io, string(precision(value)))
        _write_id_text(io, string(value))
    elseif value isa AbstractFloat
        isbitstype(type) || throw(ArgumentError(
            "collection id cannot canonically encode $(typeof(value)); " *
            "define id(collection) to return stable data",
        ))
        write(io, UInt8('f'))
        _write_id_text(io, bitstring(value))
    elseif value isa AbstractString
        write(io, UInt8('s'))
        _write_id_text(io, value)
    elseif value isa Symbol
        write(io, UInt8('y'))
        _write_id_text(io, String(value))
    elseif value isa Char
        write(io, UInt8('c'))
        _write_id_text(io, string(UInt32(value)))
    elseif value isa Dates.TimeType
        write(io, UInt8('d'))
        _write_id_text(io, string(Dates.value(value)))
    elseif value isa Type
        write(io, UInt8('t'))
        _write_id_text(io, string(value))
    elseif value isa Function || value isa Task || value isa Ptr
        throw(ArgumentError(
            "collection id cannot canonically encode $(typeof(value)); " *
            "define id(collection) to return stable data",
        ))
    else
        tracked = Base.ismutabletype(type)
        if tracked
            haskey(seen, value) && throw(ArgumentError(
                "collection id contains a reference cycle through $(typeof(value)); " *
                "define id(collection) to return stable acyclic data",
            ))
            seen[value] = nothing
        end
        try
            if value isa NamedTuple
                write(io, UInt8('N'))
                names = fieldnames(type)
                _write_id_length(io, length(names))
                for name in names
                    _write_id_text(io, String(name))
                    _write_nested_id!(io, getfield(value, name), seen)
                end
            elseif value isa Tuple
                write(io, UInt8('q'))
                _write_id_length(io, length(value))
                for element in value
                    _write_nested_id!(io, element, seen)
                end
            elseif value isa Pair
                write(io, UInt8('p'))
                _write_nested_id!(io, first(value), seen)
                _write_nested_id!(io, last(value), seen)
            elseif value isa AbstractDict
                write(io, UInt8('D'))
                entries = Vector{UInt8}[]
                for (key, item) in value
                    entry = IOBuffer()
                    _write_nested_id!(entry, key, seen)
                    _write_nested_id!(entry, item, seen)
                    push!(entries, take!(entry))
                end
                sort!(entries; by=bytes2hex)
                _write_id_length(io, length(entries))
                for entry in entries
                    _write_id_bytes(io, entry)
                end
            elseif value isa AbstractSet
                write(io, UInt8('S'))
                entries = Vector{UInt8}[]
                for item in value
                    entry = IOBuffer()
                    _write_nested_id!(entry, item, seen)
                    push!(entries, take!(entry))
                end
                sort!(entries; by=bytes2hex)
                _write_id_length(io, length(entries))
                for entry in entries
                    _write_id_bytes(io, entry)
                end
            elseif value isa AbstractArray
                write(io, UInt8('A'))
                _write_id_length(io, ndims(value))
                for dimension in size(value)
                    _write_id_length(io, dimension)
                end
                for element in value
                    _write_nested_id!(io, element, seen)
                end
            elseif isstructtype(type)
                write(io, UInt8('r'))
                names = fieldnames(type)
                _write_id_length(io, length(names))
                for name in names
                    _write_id_text(io, String(name))
                    _write_nested_id!(io, getfield(value, name), seen)
                end
            else
                throw(ArgumentError(
                    "collection id cannot canonically encode $(typeof(value)); " *
                    "define id(collection) to return stable data",
                ))
            end
        finally
            tracked && delete!(seen, value)
        end
    end
    return nothing
end

"""
Return a collection occurrence ID.

`id(collection)` supplies the user-defined value for one hierarchy level. The returned ID is the
final occurrence ID derived from the parent occurrence ID, concrete collection type, and that
user-defined value. Canonical bytes exist only while computing the digest and are not retained.
"""
function collection_record_id(
    parent_id::AbstractString,
    collection::AbstractCollection,
)::String
    io = IOBuffer()
    _write_id_text(io, "DataBrowser.collection.id.v1")
    _write_id_text(io, parent_id)
    _write_id_text(io, string(typeof(collection)))
    _write_id_bytes(io, _canonical_id_bytes(id(collection)))
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
