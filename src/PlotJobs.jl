struct PlotJob
    job_key
    plot_key
    project::AbstractProject
    target::Symbol
    target_id::String
    measurements::Vector{MeasurementInfo}
    plot_kind::Type{<:PlotKind}
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
    plot_kind::Type{<:PlotKind},
    device_params::Vector{Dict{Symbol,Any}},
    debug::Bool,
)
    return (
        project_name(project),
        target,
        target_id,
        sort([measurement.unique_id for measurement in measurements]),
        nameof(plot_kind),
        _plot_parameters_key(device_params),
        debug,
    )
end

function _run_plot_job(job::PlotJob, cancel_requested)
    return _with_cancel(cancel_requested) do
        return _plot_job_data(job.project, job.plot_kind, job.measurements; debug=job.debug)
    end
end

function _draw_plot_job(job::PlotJob, data)
    return _plot_job_figure(
        job.project,
        job.plot_kind,
        job.measurements,
        data;
        debug=job.debug,
        device_params=job.device_params,
    )
end
