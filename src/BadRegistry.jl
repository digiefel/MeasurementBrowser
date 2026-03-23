"""
BadRegistry.jl - Persistent bad device / measurement registry
"""

struct BadRegistry
    devices::Set{String}
    measurements::Set{String}
end

BadRegistry() = BadRegistry(Set{String}(), Set{String}())

const BAD_REGISTRY_FILENAME = "bad_measurements"

struct BadRegistryParseError <: Exception
    path::String
    line_number::Int
    message::String
end

function Base.showerror(io::IO, err::BadRegistryParseError)
    print(io, "Invalid $(basename(err.path)) entry at line $(err.line_number): $(err.message)")
end

struct BadRegistryIOError <: Exception
    path::String
    operation::Symbol
    cause::Exception
end

function Base.showerror(io::IO, err::BadRegistryIOError)
    print(io, "Could not $(String(err.operation)) $(err.path): ")
    showerror(io, err.cause)
end

bad_registry_path(root_path::AbstractString) = joinpath(String(root_path), BAD_REGISTRY_FILENAME)

function _parse_bad_registry_line(path::String, line_number::Int, raw_line::AbstractString)
    line = strip(String(raw_line))
    isempty(line) && return nothing
    startswith(line, '#') && return nothing

    parts = split(line; limit=2)
    length(parts) == 2 || throw(BadRegistryParseError(path, line_number, "expected 'device <path>' or 'measurement <id>'"))

    kind = parts[1]
    value = strip(parts[2])
    isempty(value) && throw(BadRegistryParseError(path, line_number, "missing value for $kind entry"))

    if kind == "device"
        device_path_tuple(value)
        return (:device, value)
    elseif kind == "measurement"
        return (:measurement, value)
    end

    throw(BadRegistryParseError(path, line_number, "unknown entry type '$kind'"))
end

function load_bad_registry(root_path::AbstractString)
    path = bad_registry_path(root_path)
    isfile(path) || return BadRegistry()

    devices = Set{String}()
    measurements = Set{String}()

    try
        for (line_number, raw_line) in enumerate(eachline(path))
            parsed = _parse_bad_registry_line(path, line_number, raw_line)
            parsed === nothing && continue
            kind, value = parsed
            if kind == :device
                push!(devices, value)
            else
                push!(measurements, value)
            end
        end
    catch err
        err isa BadRegistryParseError && rethrow()
        if err isa IOError || err isa SystemError
            throw(BadRegistryIOError(path, :read, err))
        end
        rethrow()
    end

    return BadRegistry(devices, measurements)
end

function save_bad_registry(root_path::AbstractString, registry::BadRegistry)
    path = bad_registry_path(root_path)

    try
        if isempty(registry.devices) && isempty(registry.measurements)
            rm(path; force=true)
            return
        end

        open(path, "w") do io
            for device_key in sort!(collect(registry.devices))
                println(io, "device ", device_key)
            end
            for measurement_id in sort!(collect(registry.measurements))
                println(io, "measurement ", measurement_id)
            end
        end
    catch err
        if err isa IOError || err isa SystemError
            throw(BadRegistryIOError(path, :write, err))
        end
        rethrow()
    end
end
