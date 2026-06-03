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
            df = only(data_of_measurements(job.project, [measurement]; should_cancel))
            return _ruo2_plot_data(measurement, df; debug=job.debug)
        end

        data = data_of_measurements(job.project, job.measurements; should_cancel)
        return _ruo2_combined_plot_data(job.measurements, data, job.plot_kind)
    end

    if length(job.measurements) == 1
        measurement = only(job.measurements)
        df = only(data_of_measurements(job.project, [measurement]; should_cancel))
        loaded = _ruo2_plot_data(measurement, df; debug=job.debug)
        return _analyze_ruo2_file_plot(
            job.project,
            measurement.measurement_kind,
            loaded;
            DEBUG=job.debug,
            should_cancel,
        )
    end

    data = data_of_measurements(job.project, job.measurements; should_cancel)
    loaded = _ruo2_combined_plot_data(job.measurements, data, job.plot_kind)
    return _analyze_ruo2_files_plot(
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
        return _draw_ruo2_file_plot(
            job.project,
            only(job.measurements).measurement_kind,
            data;
            device_params=only(job.device_params),
        )
    end
    return _draw_ruo2_files_plot(job.project, job.plot_kind, data; DEBUG=job.debug)
end
