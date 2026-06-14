const DEVICE_INFO_FILENAME = "device_info.txt"

"""
Send one structured progress update when a callback is present.
"""
function emit_progress(
    on_progress::Union{Nothing,Function};
    phase::Symbol,
    total_csv::Int,
    processed_csv::Int,
    loaded_measurements::Int,
    skipped_csv::Int,
    current_path::String="",
)::Nothing
    on_progress === nothing && return nothing
    on_progress((
        phase=phase,
        total_csv=total_csv,
        processed_csv=processed_csv,
        loaded_measurements=loaded_measurements,
        skipped_csv=skipped_csv,
        current_path=current_path,
    ))
    return nothing
end

"""
Count visible CSV files before a scan when exact progress totals are requested.
"""
function count_source_files(
    root_path::String;
    on_progress::Union{Nothing,Function}=nothing,
)::Int
    total = 0
    for (root, _, files) in walkdir(root_path)
        for file in files
            check_cancel()
            is_source_filename(file) || continue
            total += 1
            emit_progress(
                on_progress;
                phase=:counting,
                total_csv=total,
                processed_csv=total,
                loaded_measurements=0,
                skipped_csv=0,
                current_path=joinpath(root, file),
            )
        end
    end
    return total
end

"""
Parse one device-metadata cell into a primitive Julia value.
"""
function parse_metadata_value(text::AbstractString)::Any
    value = strip(text)
    isempty(value) && return nothing
    lowercase_value = lowercase(value)
    lowercase_value in ("true", "t", "yes", "y", "1") && return true
    lowercase_value in ("false", "f", "no", "n", "0") && return false
    for type in (Int, Float64)
        parsed = tryparse(type, value)
        parsed === nothing || return parsed
    end
    for format in (dateformat"yyyy-mm-dd", dateformat"yyyy/mm/dd")
        parsed = tryparse(Date, value, format)
        parsed === nothing || return parsed
    end
    for format in (
        dateformat"yyyy-mm-dd HH:MM:SS",
        dateformat"yyyy-mm-ddTHH:MM:SS",
    )
        parsed = tryparse(DateTime, value, format)
        parsed === nothing || return parsed
    end
    return value
end

"""
Read device metadata keyed by exact slash-separated path fragments.
"""
function load_device_metadata(
    path::AbstractString,
)::Dict{Tuple{Vararg{String}},Dict{Symbol,Any}}
    lines = readlines(path)
    isempty(lines) && return Dict{Tuple{Vararg{String}},Dict{Symbol,Any}}()
    header = strip.(split(first(lines), ','))
    length(header) >= 2 ||
        return Dict{Tuple{Vararg{String}},Dict{Symbol,Any}}()
    parameter_names = header[2:end]
    metadata = Dict{Tuple{Vararg{String}},Dict{Symbol,Any}}()
    for line in Iterators.drop(lines, 1)
        isempty(strip(line)) && continue
        cells = split(line, ',')
        path_text = strip(first(cells))
        isempty(path_text) && continue
        parameters = Dict{Symbol,Any}()
        for (offset, name) in enumerate(parameter_names)
            index = offset + 1
            index > length(cells) && continue
            value = parse_metadata_value(cells[index])
            value === nothing || (parameters[Symbol(strip(name))] = value)
        end
        metadata[Tuple(filter(!isempty, split(path_text, '/')))] = parameters
    end
    return metadata
end

"""
Merge every metadata fragment matching a parsed device path.

Longer path fragments are applied later and therefore override broader entries.
"""
function matching_device_parameters(
    metadata::Dict{Tuple{Vararg{String}},Dict{Symbol,Any}},
    location::Vector{String},
)::Union{Nothing,Dict{Symbol,Any}}
    merged = Dict{Symbol,Any}()
    for width in eachindex(location)
        for start in 1:(length(location) - width + 1)
            parameters = get(
                metadata,
                Tuple(location[start:(start + width - 1)]),
                nothing,
            )
            parameters === nothing || merge!(merged, parameters)
        end
    end
    return isempty(merged) ? nothing : merged
end

device_info_path(root_path::AbstractString)::String =
    joinpath(root_path, DEVICE_INFO_FILENAME)

has_device_metadata(root_path::AbstractString)::Bool =
    isfile(device_info_path(root_path))

"""
Read optional source-root device metadata for one scan.
"""
function load_scan_metadata(
    root_path::String,
)::Union{Nothing,Dict{Tuple{Vararg{String}},Dict{Symbol,Any}}}
    path = device_info_path(root_path)
    return isfile(path) ? load_device_metadata(path) : nothing
end

"""
Create the short title shown for a measurement in browser panels and plots.
"""
function build_clean_title(
    project::AbstractProject,
    filename::String,
    measurement_kind::Symbol,
    device_info::DeviceInfo,
)::String
    measurement_label = kind_label(project, measurement_kind)
    device_label = device_path_label(project, device_info)
    parts = filter(!isempty, (
        measurement_label == "Unknown" ? "" : measurement_label,
        device_label,
    ))
    return isempty(parts) ?
        strip(replace(filename, r"\.csv$" => "")) :
        join(parts, " ")
end

"""
Interpret one indexed file and apply matching source-root device metadata.
"""
function interpret_measurements(
    project::AbstractProject,
    file::SourceFile,
    metadata::Union{Nothing,Dict{Tuple{Vararg{String}},Dict{Symbol,Any}}},
)::Vector{MeasurementInfo}
    measurements = interpret_file(project, file)
    metadata === nothing && return measurements
    for measurement in measurements
        parameters = matching_device_parameters(metadata, measurement.device_info.location)
        parameters === nothing || merge!(measurement.device_info.parameters, parameters)
    end
    return measurements
end

"""
Interpret and analyze one physical file using the same project path as a full scan.
"""
function measurements_for_file(
    project::AbstractProject,
    filepath::AbstractString;
    meta::Union{Nothing,Dict{Tuple{Vararg{String}},Dict{Symbol,Any}}}=nothing,
)::Vector{MeasurementInfo}
    file = index_source_file(filepath)
    measurements = interpret_measurements(project, file, meta)
    compute_and_add_measurement_stats!(project, measurements, [file])
    return measurements
end

"""
Interpret indexed files concurrently and stream each successful measurement batch.

Failures remain attached to their physical files and do not stop unrelated files.
"""
function interpret_source_files(
    project::AbstractProject,
    files::Vector{SourceFile},
    metadata::Union{Nothing,Dict{Tuple{Vararg{String}},Dict{Symbol,Any}}};
    on_result::Function,
    on_measurements::Union{Nothing,Function}=nothing,
    on_progress::Union{Nothing,Function}=nothing,
)::NamedTuple
    processed_count = Base.Threads.Atomic{Int}(0)
    measurement_count = Base.Threads.Atomic{Int}(0)
    skipped_count = Base.Threads.Atomic{Int}(0)
    callback_lock = ReentrantLock()
    worker_limit = Base.Semaphore(max(1, Base.Threads.nthreads()))
    cancel_requested = get(task_local_storage(), CANCEL_CALLBACK_KEY, nothing)
    failures = MeasurementAnalysisFailure[]

    @sync for (index, file) in pairs(files)
        check_cancel()
        Base.acquire(worker_limit)
        Base.Threads.@spawn try
            with_cancel(cancel_requested) do
                check_cancel()
                measurements = try
                    interpret_measurements(project, file, metadata)
                catch error
                    is_job_cancelled(error) && rethrow()
                    backtrace = catch_backtrace()
                    @error(
                        "Source file interpretation failed",
                        project=project_name(project),
                        file=file.filepath,
                        exception=(error, backtrace),
                    )
                    lock(callback_lock) do
                        push!(failures, MeasurementAnalysisFailure(
                            file.filepath,
                            "",
                            sprint(showerror, error),
                        ))
                    end
                    MeasurementInfo[]
                end
                check_cancel()
                on_result(index, SourceFile(file, measurements))
                isempty(measurements) || on_measurements === nothing ||
                    lock(() -> on_measurements(measurements), callback_lock)
                isempty(measurements) ?
                    Base.Threads.atomic_add!(skipped_count, 1) :
                    Base.Threads.atomic_add!(measurement_count, length(measurements))
                processed = Base.Threads.atomic_add!(processed_count, 1) + 1
                on_progress === nothing || lock(callback_lock) do
                    on_progress((
                        total_csv=length(files),
                        processed_csv=processed,
                        loaded_measurements=measurement_count[],
                        skipped_csv=skipped_count[],
                        current_path=file.filepath,
                    ))
                end
            end
        finally
            Base.release(worker_limit)
        end
    end

    return (
        processed_csv=processed_count[],
        loaded_measurements=measurement_count[],
        skipped_csv=skipped_count[],
        total_csv=length(files),
        failures,
    )
end

"""Whether the freshly indexed files exactly match the cached set, by path and fingerprint."""
function _files_unchanged(
    files::Vector{SourceFile},
    cached::Dict{String,SourceFile},
)::Bool
    length(files) == length(cached) || return false
    for file in files
        cached_file = get(cached, file.filepath, nothing)
        (cached_file === nothing || cached_file.fingerprint != file.fingerprint) && return false
    end
    return true
end

"""
Scan one source root into its complete measurement hierarchy.

Measurements stream through `on_measurements` while interpretation is still running. Project
statistics run only after all files are known, allowing cross-file calculations.
"""
function scan_source(
    root_path::String;
    project::AbstractProject,
    cached_files::Union{Nothing,Dict{String,SourceFile}}=nothing,
    cached_source::Union{Nothing,SourceScan}=nothing,
    on_progress::Union{Nothing,Function}=nothing,
    on_measurements::Union{Nothing,Function}=nothing,
    count_first::Bool=false,
)::SourceScan
    root = normpath(abspath(expanduser(root_path)))
    metadata = load_scan_metadata(root)
    count_first && count_source_files(root; on_progress)
    emit_progress(
        on_progress;
        phase=:discovering,
        total_csv=0,
        processed_csv=0,
        loaded_measurements=0,
        skipped_csv=0,
    )
    files = collect_source_files(
        root;
        on_file=(file, count) -> emit_progress(
            on_progress;
            phase=:discovering,
            total_csv=0,
            processed_csv=count,
            loaded_measurements=0,
            skipped_csv=0,
            current_path=file.filepath,
        ),
    )
    check_cancel()
    # Fingerprinting is a cheap stat() per file; reading and analyzing is what's slow. When nothing
    # changed since the cache was written, the cache is authoritative — skip all reads and analysis.
    if cached_source !== nothing && cached_files !== nothing && _files_unchanged(files, cached_files)
        emit_progress(
            on_progress;
            phase=:analyzing,
            total_csv=1,
            processed_csv=1,
            loaded_measurements=length(cached_source.hierarchy.all_measurements),
            skipped_csv=cached_source.hierarchy.skipped_count,
            current_path=root,
        )
        return cached_source
    end
    scanned_files = Vector{SourceFile}(undef, length(files))
    emit_progress(
        on_progress;
        phase=:scanning,
        total_csv=length(files),
        processed_csv=0,
        loaded_measurements=0,
        skipped_csv=0,
    )
    summary = interpret_source_files(
        project,
        files,
        metadata;
        on_result=(index, file) -> (scanned_files[index] = file),
        on_measurements,
        on_progress=progress -> emit_progress(
            on_progress;
            phase=:scanning,
            total_csv=progress.total_csv,
            processed_csv=progress.processed_csv,
            loaded_measurements=progress.loaded_measurements,
            skipped_csv=progress.skipped_csv,
            current_path=progress.current_path,
        ),
    )
    check_cancel()
    measurements = reduce(
        append!,
        (file.measurements for file in scanned_files);
        init=MeasurementInfo[],
    )
    analysis_progress_seen = Ref(false)
    emit_progress(
        on_progress;
        phase=:analyzing,
        total_csv=0,
        processed_csv=0,
        loaded_measurements=length(measurements),
        skipped_csv=summary.skipped_csv,
        current_path=root,
    )
    analysis_failures = compute_and_add_measurement_stats!(
        project,
        measurements,
        scanned_files;
        on_progress=progress -> begin
            analysis_progress_seen[] = true
            emit_progress(
                on_progress;
                phase=:analyzing,
                total_csv=progress.total,
                processed_csv=progress.processed,
                loaded_measurements=length(measurements),
                skipped_csv=summary.skipped_csv,
                current_path=progress.current_path,
            )
        end,
    )
    analysis_progress_seen[] || emit_progress(
        on_progress;
        phase=:analyzing,
        total_csv=1,
        processed_csv=1,
        loaded_measurements=length(measurements),
        skipped_csv=summary.skipped_csv,
        current_path=root,
    )
    check_cancel()
    hierarchy = MeasurementHierarchy(
        measurements,
        root,
        metadata !== nothing,
        project,
        summary.skipped_csv,
    )
    append!(summary.failures, analysis_failures)
    return SourceScan(
        root,
        project,
        scanned_files,
        hierarchy,
        summary.failures,
    )
end
