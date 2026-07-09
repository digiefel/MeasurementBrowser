using Printf
using Statistics: mean
import ModernGL as gl

"""Retain a bounded history of one performance sample."""
function _append_perf_sample!(
    state::BrowserState,
    key::Symbol,
    duration_ms::Float64,
    alloc_bytes::Int,
)::Nothing
    samples = get!(() -> Float64[], state.performance.timings, key)
    allocations = get!(() -> Int[], state.performance.allocations, key)
    push!(samples, duration_ms)
    length(samples) > 400 && popfirst!(samples)
    push!(allocations, alloc_bytes)
    length(allocations) > 400 && popfirst!(allocations)
    return nothing
end

"""Measure one browser operation and retain its latest timing and allocation samples."""
function _time!(
    f::Function,
    state::BrowserState,
    key::Symbol,
)::Nothing
    started_ns = time_ns()
    allocated_bytes = @allocated f()
    _append_perf_sample!(
        state,
        key,
        (time_ns() - started_ns) / 1e6,
        allocated_bytes,
    )
    return nothing
end

"""Format a byte count for the performance window."""
function _fmt_bytes(bytes::Integer)::String
    gib = bytes / (1024^3)
    gib >= 1 && return @sprintf("%.2f GiB", gib)
    mib = bytes / (1024^2)
    mib >= 1 && return @sprintf("%.0f MiB", mib)
    bytes >= 1024 && return @sprintf("%.1f KB", bytes / 1024)
    return "$bytes B"
end

"""Read identifying strings from the active OpenGL context."""
function _gl_info()::Dict{Symbol,String}
    try
        return Dict(
            :vendor => unsafe_string(gl.glGetString(gl.GL_VENDOR)),
            :renderer => unsafe_string(gl.glGetString(gl.GL_RENDERER)),
            :version => unsafe_string(gl.glGetString(gl.GL_VERSION)),
            :sl => unsafe_string(gl.glGetString(gl.GL_SHADING_LANGUAGE_VERSION)),
        )
    catch err
        @warn "GL info query failed" error = err
        return Dict{Symbol,String}()
    end
end

"""Write the collected performance samples to the debug log at shutdown."""
function _print_perf_summary(state::BrowserState)::Nothing
    @debug begin
        gl_info = state.performance.gl_info
        timings = state.performance.timings
        allocations = state.performance.allocations
        message = """\n
        ==== Performance Summary ====
        GL Vendor:   $(get(gl_info, :vendor, "?"))
        GL Renderer: $(get(gl_info, :renderer, "?"))
        GL Version:  $(get(gl_info, :version, "?"))
        """
        if !isempty(timings)
            message *= @sprintf(
                "%-12s %5s %9s %9s %9s %12s %12s\n",
                "Key",
                "n",
                "Mean(ms)",
                "Max(ms)",
                "Last(ms)",
                "AllocMean(KB)",
                "AllocLast(KB)",
            )
        end
        for key in sort(collect(keys(timings)))
            samples = timings[key]
            isempty(samples) && continue
            allocated = get(allocations, key, Int[])
            message *= @sprintf(
                "%-12s %5d %9.2f %9.2f %9.2f %12.1f %12.1f\n",
                String(key),
                length(samples),
                mean(samples),
                maximum(samples),
                samples[end],
                isempty(allocated) ? 0.0 : mean(allocated) / 1024,
                isempty(allocated) ? 0.0 : allocated[end] / 1024,
            )
        end
        message * "=============================="
    end
    return nothing
end
