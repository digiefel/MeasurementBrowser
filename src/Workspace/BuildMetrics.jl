"""Cheap always-on aggregate cache-write metrics for status and unbiased benchmarks."""
mutable struct BuildMetrics
    interpreted_write_ns::Base.Threads.Atomic{Int64}
    interpreted_writes::Base.Threads.Atomic{Int64}
    processed_write_ns::Base.Threads.Atomic{Int64}
    processed_writes::Base.Threads.Atomic{Int64}
    stats_write_ns::Base.Threads.Atomic{Int64}
    stats_writes::Base.Threads.Atomic{Int64}
end

"""Construct zeroed build metrics."""
function BuildMetrics()::BuildMetrics
    counter() = Base.Threads.Atomic{Int64}(0)
    return BuildMetrics(counter(), counter(), counter(), counter(), counter(), counter())
end

"""Reset all aggregate counters before one workspace build."""
function reset_build_metrics!(metrics::BuildMetrics)::Nothing
    for counter in (
        metrics.interpreted_write_ns,
        metrics.interpreted_writes,
        metrics.processed_write_ns,
        metrics.processed_writes,
        metrics.stats_write_ns,
        metrics.stats_writes,
    )
        counter[] = 0
    end
    return nothing
end

"""Add one completed cache phase to an aggregate counter pair."""
@inline function record_cache_phase!(
    total_ns::Base.Threads.Atomic{Int64},
    calls::Base.Threads.Atomic{Int64},
    t0_ns::UInt64,
)::Nothing
    Base.Threads.atomic_add!(total_ns, Int64(time_ns() - t0_ns))
    Base.Threads.atomic_add!(calls, Int64(1))
    return nothing
end
