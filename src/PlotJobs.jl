struct PlotJob
    job_key
    plot_key
    project::AbstractProject
    target::Symbol
    target_id::String
    measurements::Vector{MeasurementInfo}
    plot_kind::Symbol
    device_params::Vector{Dict{Symbol,Any}}
    cache_identity::Union{Nothing,ProjectCacheIdentity}
    cache_version
    debug::Bool
end

function _plot_parameters_key(params::Vector{Dict{Symbol,Any}})
    return Tuple(
        Tuple((key, repr(value)) for (key, value) in sort(collect(p); by=x -> String(first(x))))
        for p in params
    )
end

function _plot_job_key(
    project::AbstractProject,
    target::Symbol,
    target_id::String,
    measurements::Vector{MeasurementInfo},
    plot_kind::Symbol,
    device_params::Vector{Dict{Symbol,Any}},
    cache_version,
    debug::Bool,
)
    return (
        project_name(project),
        target,
        target_id,
        sort([measurement.id for measurement in measurements]),
        plot_kind,
        _plot_parameters_key(device_params),
        cache_version,
        debug,
    )
end

function _run_plot_job(job::PlotJob, should_cancel)
    if length(job.measurements) == 1
        measurement = only(job.measurements)
        params = only(job.device_params)
        if job.debug
            loaded = load_plot_for_file(
                job.project,
                measurement.filepath,
                measurement.measurement_kind;
                device_params=params,
                DEBUG=job.debug,
                should_cancel,
            )
            return analyze_plot_for_file(
                job.project,
                measurement.measurement_kind,
                loaded;
                device_params=params,
                DEBUG=job.debug,
                should_cancel,
            )
        end
        job.cache_identity isa ProjectCacheIdentity ||
            error("Plot job for '$(measurement.filepath)' is missing cache identity")
        return _measurement_group_for_cached_plot(job.cache_identity, measurement)
    end

    paths = [measurement.filepath for measurement in job.measurements]
    loaded = load_plot_for_files(
        job.project,
        paths,
        job.plot_kind;
        device_params_list=job.device_params,
        DEBUG=job.debug,
        should_cancel,
    )
    return analyze_plot_for_files(
        job.project,
        job.plot_kind,
        loaded;
        DEBUG=job.debug,
        should_cancel,
    )
end

function _draw_plot_job(job::PlotJob, data)
    if length(job.measurements) == 1
        return draw_plot_for_file(
            job.project,
            only(job.measurements).measurement_kind,
            data;
            device_params=only(job.device_params),
        )
    end
    return draw_plot_for_files(job.project, job.plot_kind, data; DEBUG=job.debug)
end
