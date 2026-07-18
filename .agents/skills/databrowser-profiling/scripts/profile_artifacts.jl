module ProfileArtifacts

using Profile
using Serialization

const FORMAT_VERSION = 1
const FrameRecord = NamedTuple{
    (:func, :file, :line, :from_c, :inlined, :pointer),
    Tuple{String,String,Int,Bool,Bool,UInt64},
}

"""
    profile_sample_count(data::AbstractVector{<:Unsigned}) -> Int

Return the number of sampled backtraces in Julia `Profile` data.
"""
function profile_sample_count(data::AbstractVector{<:Unsigned})::Int
    isempty(data) && return 0
    if isdefined(Profile, :has_meta) && Profile.has_meta(data)
        return count(i -> Profile.is_block_end(data, i), eachindex(data))
    end
    return count(iszero, data)
end

function frame_record(frame::Base.StackTraces.StackFrame)::FrameRecord
    return (
        func=string(frame.func),
        file=string(frame.file),
        line=Int(frame.line),
        from_c=frame.from_c,
        inlined=frame.inlined,
        pointer=UInt64(frame.pointer),
    )
end

"""
    save_profile(path, data, lidict; sampling_delay_ns) -> String

Save Julia profile samples and portable line records without serializing method instances.
"""
function save_profile(
    path::AbstractString,
    data::AbstractVector{<:Unsigned},
    lidict::AbstractDict;
    sampling_delay_ns::Integer=ccall(:jl_profile_delay_nsec, UInt64, ()),
)::String
    frames = Dict{UInt64,Vector{FrameRecord}}()
    for (instruction_pointer, line_info) in lidict
        entries = line_info isa AbstractVector ? line_info : (line_info,)
        frames[UInt64(instruction_pointer)] =
            FrameRecord[frame_record(frame) for frame in entries]
    end
    payload = (
        format_version=FORMAT_VERSION,
        sampling_delay_ns=UInt64(sampling_delay_ns),
        data=UInt64[data...],
        frames=frames,
    )
    open(path, "w") do io
        serialize(io, payload)
    end
    return String(path)
end

function stack_frame(record::FrameRecord)::Base.StackTraces.StackFrame
    return Base.StackTraces.StackFrame(
        Symbol(record.func),
        Symbol(record.file),
        record.line,
        nothing,
        record.from_c,
        record.inlined,
        record.pointer,
    )
end

"""
    load_profile(path) -> NamedTuple

Load a profile written by [`save_profile`](@ref) and reconstruct its `Profile.LineInfoDict`.
"""
function load_profile(path::AbstractString)
    payload = open(deserialize, path)
    payload.format_version == FORMAT_VERSION ||
        error("Unsupported profile artifact format $(payload.format_version)")
    lidict = Profile.LineInfoDict()
    for (instruction_pointer, records) in payload.frames
        lidict[instruction_pointer] =
            Base.StackTraces.StackFrame[stack_frame(record) for record in records]
    end
    return (
        data=payload.data,
        lidict=lidict,
        sampling_delay_ns=payload.sampling_delay_ns,
    )
end

end
