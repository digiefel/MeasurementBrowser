"""Return source-specific failures accumulated outside `data_items` itself."""
source_analysis_failures(::Project, ::AbstractDataSource)::Vector{ItemFailure} = ItemFailure[]

"""
One completed source-item pass.

`records` are retained by the index. `interpreted_items` carry effective item data until the cache
writer accepts it and processing is queued.
"""
struct SourceItemInterpretation
    records::Vector{ItemRecord}
    interpreted_items::Vector{AbstractDataItem}
    failures::Vector{ItemFailure}
end

"""
Send one structured progress update when a callback is present.
"""
function emit_progress(
    on_progress::Union{Nothing,Function};
    phase::Symbol,
    total_source_items::Int,
    processed_source_items::Int,
    loaded_items::Int,
    skipped_source_items::Int,
    current_source_item::String="",
)::Nothing
    on_progress === nothing && return nothing
    on_progress((
        phase=phase,
        total_source_items=total_source_items,
        processed_source_items=processed_source_items,
        loaded_items=loaded_items,
        skipped_source_items=skipped_source_items,
        current_source_item=current_source_item,
    ))
    return nothing
end

"""Return current source-item fingerprints keyed by source-item id."""
function fingerprints_by_source_item(
    items::Vector{<:AbstractDataSourceItem},
)::Dict{String,Any}
    return Dict(source_item_id(item) => fingerprint(item) for item in items)
end

"""Whether every current source item has a non-null fingerprint matching the cached scan."""
function source_unchanged(
    source::AbstractDataSource,
    fingerprints::Dict{String,Any},
    cached::SourceScan,
)::Bool
    cached.source_id == source_id(source) || return false
    all(!isnothing, values(fingerprints)) || return false
    return cached.source_item_fingerprints == fingerprints
end

"""
Interpret every logical data item produced by one source item.

The source item is read only by `data_items`. Each result is normalized to a package `DataItem` with
effective parameters. Processing and statistics belong to the workspace processing queue.
"""
function interpret_source_item(
    project::Project,
    source::AbstractDataSource,
    source_item::AbstractDataSourceItem,
    profiler::Union{Nothing,Profiling.ProfileSession}=nothing,
)::SourceItemInterpretation
    source_item_id_value = source_item_id(source_item)
    source_item_path_value = source_item_path(source_item)
    source_item_label_value = if source_item_path_value !== nothing &&
                                 isabspath(source_id(source))
        relpath(source_item_path_value, source_id(source))
    else
        source_item_label(source_item)
    end
    source_started = time_ns()
    handles = try
        Profiling.@profile_span profiler :project :interpret_source_item Profiling.ProfileAttributes(
            source_id=source_item_id_value,
        ) begin
            data_items(project, source, source_item)::Vector{<:AbstractDataItem}
        end
    catch
        finish_source_profile!(
            project, source_item_id_value, source_item_label_value, source_item_path_value,
            :unmatched, 0,
            (time_ns() - source_started) / 1e9, Set([Base.Threads.threadid()]))
        rethrow()
    end
    item_count = length(handles)
    item_kinds = unique(kind(handle) for handle in handles)
    source_kind = length(item_kinds) == 1 ? only(item_kinds) :
        isempty(item_kinds) ? :unmatched : :mixed
    records = Vector{ItemRecord}(undef, item_count)
    interpreted_items = Vector{AbstractDataItem}(undef, item_count)
    for (index, handle) in pairs(handles)
        check_cancel()
        record = ItemRecord(handle; source_item,
            parameters=_effective_parameters(source, collection(handle), parameters(handle)))
        records[index] = record
        interpreted_items[index] = DataItem(record, item_data(handle))
    end
    finish_source_profile!(
        project, source_item_id_value, source_item_label_value, source_item_path_value,
        source_kind, item_count,
        (time_ns() - source_started) / 1e9, Set([Base.Threads.threadid()]))
    return SourceItemInterpretation(records, interpreted_items, ItemFailure[])
end

"""
Interpret source items concurrently and stream each successful record batch.
"""
function interpret_source_items(
    project::Project,
    source::AbstractDataSource,
    source_items::Vector{<:AbstractDataSourceItem};
    on_items::Union{Nothing,Function}=nothing,
    on_source_item::Union{Nothing,Function}=nothing,
    on_kept_source_item::Union{Nothing,Function}=nothing,
    on_progress::Union{Nothing,Function}=nothing,
    on_failure::Union{Nothing,Function}=nothing,
    unchanged_source_item_ids::Set{String}=Set{String}(),
    cached_records_by_source_item::Dict{String,Vector{ItemRecord}}=
        Dict{String,Vector{ItemRecord}}(),
    profiler::Union{Nothing,Profiling.ProfileSession}=nothing,
)::NamedTuple
    processed_count::Int = 0
    item_count::Int = 0
    skipped_count::Int = 0
    cancel_requested = get(task_local_storage(), CANCEL_CALLBACK_KEY, nothing)
    failures::Vector{ItemFailure} = ItemFailure[]
    records_by_position::Vector{Vector{ItemRecord}} =
        Vector{Vector{ItemRecord}}(undef, length(source_items))
    isempty(source_items) && return (
        processed_source_items=0,
        loaded_items=0,
        skipped_source_items=0,
        total_source_items=0,
        failures,
        records=ItemRecord[],
    )

    worker_count = min(max(1, Base.Threads.nthreads() ÷ 2), length(source_items))
    work = Channel{Tuple{Int,eltype(source_items)}}(length(source_items))
    results = Channel{NamedTuple}(worker_count)
    for (index, source_item) in pairs(source_items)
        put!(work, (index, source_item))
    end
    close(work)

    received = falses(length(source_items))
    @sync begin
        for _ in 1:worker_count
            Base.Threads.@spawn with_cancel(cancel_requested) do
                for (index, source_item) in work
                    source_item_id_value = source_item_id(source_item)
                    if source_item_id_value in unchanged_source_item_ids
                        # Fingerprint unchanged: reuse the cached records without reading the origin
                        # or recomputing process/stats. No loaded data is produced, so the cache
                        # writer never touches this source item's existing rows.
                        put!(results, (
                            kind=:unchanged,
                            index,
                            source_item_id=source_item_id_value,
                            interpretation=SourceItemInterpretation(
                                get(cached_records_by_source_item, source_item_id_value,
                                    ItemRecord[]),
                                AbstractDataItem[],
                                ItemFailure[],
                            ),
                            failure=nothing,
                        ))
                        continue
                    end
                    interpretation = try
                        check_cancel()
                        interpreted = interpret_source_item(
                            project, source, source_item, profiler)
                        check_cancel()
                        interpreted
                    catch error
                        if is_job_cancelled(error)
                            put!(results, (
                                kind=:cancelled,
                                index,
                                source_item_id=source_item_id_value,
                                interpretation=SourceItemInterpretation(
                                    ItemRecord[],
                                    AbstractDataItem[],
                                    ItemFailure[],
                                ),
                                failure=nothing,
                            ))
                            return nothing
                        end
                        backtrace = catch_backtrace()
                        @error(
                            "Source item interpretation failed",
                            source=source_label(source),
                            source_item=source_item_id_value,
                            exception=(error, backtrace),
                        )
                        put!(results, (
                            kind=:failure,
                            index,
                            source_item_id=source_item_id_value,
                            interpretation=SourceItemInterpretation(
                                ItemRecord[],
                                AbstractDataItem[],
                                ItemFailure[],
                            ),
                            failure=ItemFailure(
                                source_item_id_value,
                                "",
                                sprint(showerror, error),
                            ),
                        ))
                        continue
                    end
                    put!(results, (
                        kind=:ok,
                        index,
                        source_item_id=source_item_id_value,
                        interpretation,
                        failure=nothing,
                    ))
                end
            end
        end

        for _ in eachindex(source_items)
            result = take!(results)
            result.kind === :cancelled && throw(JobCancelled())
            index = result.index
            1 <= index <= length(source_items) || error(
                "Scan worker returned invalid source item index $index",
            )
            received[index] && error(
                "Scan worker returned duplicate result for source item " *
                "'$(result.source_item_id)' at index $index",
            )
            expected_source_item_id = source_item_id(source_items[index])
            result.source_item_id == expected_source_item_id || error(
                "Scan worker result for index $index belongs to source item " *
                "'$(result.source_item_id)', expected '$expected_source_item_id'",
            )
            received[index] = true
            if result.failure !== nothing
                push!(failures, result.failure)
                on_failure === nothing || on_failure(result.failure)
            end
            interpretation = result.interpretation
            append!(failures, interpretation.failures)
            on_failure === nothing || foreach(on_failure, interpretation.failures)
            records = interpretation.records
            records_by_position[index] = records
            isempty(records) || on_items === nothing || on_items(records)
            if result.kind === :unchanged
                # Preserve the existing cached rows for this source item; never re-cache it.
                on_kept_source_item === nothing ||
                    on_kept_source_item(result.source_item_id)
            else
                isempty(records) || on_source_item === nothing ||
                    on_source_item(interpretation)
            end
            isempty(records) ? (skipped_count += 1) :
                (item_count += length(records))
            processed_count += 1
            on_progress === nothing || on_progress((
                total_source_items=length(source_items),
                processed_source_items=processed_count,
                loaded_items=item_count,
                skipped_source_items=skipped_count,
                current_source_item=result.source_item_id,
            ))
        end
    end
    all(received) || error("Scan completed without receiving every worker result")

    records = reduce(append!, records_by_position; init=ItemRecord[])
    return (
        processed_source_items=processed_count,
        loaded_items=item_count,
        skipped_source_items=skipped_count,
        total_source_items=length(source_items),
        failures,
        records,
    )
end

"""
Decide which source items can be reused from a cached scan without re-reading them.

A source item is reusable when it is still discovered and its non-null fingerprint exactly matches
the cached one. Returns the cached records grouped by source item and the set of reusable source-item
ids (including reusable items that previously produced no records, so the scan keeps their bare/error
rows).
"""
function _incremental_reuse_plan(
    cached_source::SourceScan,
    fingerprints::Dict{String,Any},
)::Tuple{Dict{String,Vector{ItemRecord}},Set{String}}
    cached_fingerprints = cached_source.source_item_fingerprints
    records_by_source_item = Dict{String,Vector{ItemRecord}}()
    for record in cached_source.hierarchy.all_items
        push!(
            get!(() -> ItemRecord[], records_by_source_item, record.source_item_id),
            record,
        )
    end
    unchanged = Set{String}()
    for (source_item, fingerprint) in fingerprints
        fingerprint === nothing && continue
        haskey(cached_fingerprints, source_item) || continue
        isequal(cached_fingerprints[source_item], fingerprint) || continue
        push!(unchanged, source_item)
    end
    return records_by_source_item, unchanged
end

"""
Scan one source into its complete item hierarchy.
"""
function scan_source(
    project::Project,
    source::AbstractDataSource;
    cached_source::Union{Nothing,SourceScan}=nothing,
    on_progress::Union{Nothing,Function}=nothing,
    on_items::Union{Nothing,Function}=nothing,
    on_source_item::Union{Nothing,Function}=nothing,
    on_kept_source_item::Union{Nothing,Function}=nothing,
    on_failure::Union{Nothing,Function}=nothing,
    count_first::Bool=false,
    profiler::Union{Nothing,Profiling.ProfileSession}=nothing,
)::SourceScan
    reset_scan_profile!(project)
    count_first && emit_progress(
        on_progress;
        phase=:counting,
        total_source_items=0,
        processed_source_items=0,
        loaded_items=0,
        skipped_source_items=0,
    )
    emit_progress(
        on_progress;
        phase=:discovering,
        total_source_items=0,
        processed_source_items=0,
        loaded_items=0,
        skipped_source_items=0,
    )
    discovered = Profiling.@profile_span profiler :source :discover Profiling.ProfileAttributes(
        source_id=source_id(source),
    ) begin
        source_items(source)
    end
    fingerprints = fingerprints_by_source_item(discovered)
    check_cancel()
    cached_fingerprints_match = cached_source !== nothing &&
        source_unchanged(source, fingerprints, cached_source)
    if cached_fingerprints_match && !has_collection_parameters(source)
        emit_progress(
            on_progress;
            phase=:scanning,
            total_source_items=length(discovered),
            processed_source_items=length(discovered),
            loaded_items=length(cached_source.hierarchy.all_items),
            skipped_source_items=cached_source.hierarchy.skipped_count,
        )
        return cached_source
    end

    emit_progress(
        on_progress;
        phase=:scanning,
        total_source_items=length(discovered),
        processed_source_items=0,
        loaded_items=0,
        skipped_source_items=0,
    )
    # Reuse unchanged source items from the cache (read once, never again across sessions). Collection
    # parameters bake into a record's stored parameters, so when they are present we re-interpret
    # everything rather than risk serving stale effective parameters from a reused record.
    cached_records_by_source_item, unchanged_source_item_ids =
        if cached_source !== nothing && !has_collection_parameters(source)
            _incremental_reuse_plan(cached_source, fingerprints)
        else
            (Dict{String,Vector{ItemRecord}}(), Set{String}())
        end
    summary = interpret_source_items(
        project,
        source,
        discovered;
        on_items,
        on_source_item,
        on_kept_source_item,
        on_failure,
        unchanged_source_item_ids,
        cached_records_by_source_item,
        profiler,
        on_progress=progress -> emit_progress(
            on_progress;
            phase=:scanning,
            total_source_items=progress.total_source_items,
            processed_source_items=progress.processed_source_items,
            loaded_items=progress.loaded_items,
            skipped_source_items=progress.skipped_source_items,
            current_source_item=progress.current_source_item,
        ),
    )
    check_cancel()
    seen_ids = Set{String}()
    for record in summary.records
        record.id in seen_ids && error(
            "Duplicate item id '$(record.id)' produced while scanning source item " *
            "'$(record.source_item_id)'. Item ids must be unique within a source.",
        )
        push!(seen_ids, record.id)
    end
    hierarchy = Hierarchy(
        summary.records,
        source,
        summary.skipped_source_items,
    )
    if cached_source !== nothing
        append!(
            summary.failures,
            (
                failure
                for failure in cached_source.analysis_failures
                if failure.source_item_id in unchanged_source_item_ids
            ),
        )
    end
    append!(summary.failures, source_analysis_failures(project, source))
    source_scan = SourceScan(
        source_id(source),
        source_label(source),
        fingerprints,
        hierarchy,
        summary.failures,
    )
    if cached_fingerprints_match &&
       source_parameter_state(cached_source) == source_parameter_state(source_scan)
        return cached_source
    end
    return source_scan
end
