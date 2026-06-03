using DataFrames

function load_source_data(
    ::AbstractProject,
    source_file::SourceFile;
    measurement::Union{Nothing,MeasurementInfo}=nothing,
    should_cancel::Union{Nothing,Function}=nothing,
)::DataFrame
    error("No source data loader for $(source_file.filepath)")
end

function data_of_measurements(
    project::AbstractProject,
    measurements::Vector{MeasurementInfo};
    should_cancel::Union{Nothing,Function}=nothing,
)::Vector{DataFrame}
    data = DataFrame[]
    sizehint!(data, length(measurements))
    for measurement in measurements
        _check_plot_cancel(should_cancel)
        source_file = index_source_file(measurement.filepath)
        push!(data, load_source_data(project, source_file; measurement, should_cancel))
    end
    return data
end
