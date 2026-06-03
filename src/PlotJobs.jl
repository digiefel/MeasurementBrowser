struct PlotJob
    job_key
    plot_key
    project::AbstractProject
    target::Symbol
    target_id::String
    measurements::Vector{MeasurementInfo}
    plot_kind::Symbol
    device_params::Vector{Dict{Symbol,Any}}
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
    debug::Bool,
)
    return (
        project_name(project),
        target,
        target_id,
        sort([measurement.unique_id for measurement in measurements]),
        plot_kind,
        _plot_parameters_key(device_params),
        debug,
    )
end

function _uses_new_plot_api(job::PlotJob)::Bool
    return job.project isa TASEProject ||
        (job.project isa RuO2Project && job.plot_kind === :cvsweep)
end

function _run_plot_job(job::PlotJob, should_cancel)
    if _uses_new_plot_api(job)
        return nothing
    end

    if job.debug
        if length(job.measurements) == 1
            measurement = only(job.measurements)
            params = only(job.device_params)
            return load_plot_for_file(
                job.project,
                measurement.filepath,
                measurement.measurement_kind;
                device_params=params,
                DEBUG=job.debug,
                should_cancel,
            )
        end

        paths = [measurement.filepath for measurement in job.measurements]
        return load_plot_for_files(
            job.project,
            paths,
            job.plot_kind;
            device_params_list=job.device_params,
            DEBUG=job.debug,
            should_cancel,
        )
    end

    if length(job.measurements) == 1
        measurement = only(job.measurements)
        params = only(job.device_params)
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
    if _uses_new_plot_api(job)
        job.debug && error("Debug plots are not implemented for $(project_name(job.project)) $(job.plot_kind)")
        fig = setup_plot(job.project, job.plot_kind, job.measurements)
        plot_data!(job.project, job.plot_kind, job.measurements, fig)
        return fig
    end

    if job.debug
        return debug_plot(
            job.project,
            job.measurements,
            data;
            device_params=job.device_params,
            plot_kind=job.plot_kind,
        )
    end

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
