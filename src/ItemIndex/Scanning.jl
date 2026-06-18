"""Return source-specific failures accumulated outside `data_items` itself."""
function source_analysis_failures end

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

"""Return the metadata table a source wants applied during indexing."""
source_collection_metadata(::AbstractDataSource) = nothing

apply_collection_metadata!(::Vector{ItemRecord}, ::Nothing)::Nothing = nothing

"""Return the current source-item fingerprints keyed by source-item id."""
function source_item_fingerprints(
    items::Vector{<:AbstractDataSourceItem},
)::Dict{String,Any}
    return Dict(source_item_id(item) => source_item_fingerprint(item) for item in items)
end

"""Whether a cached scan still matches the source-level and source-item fingerprints."""
function source_unchanged(
    source::AbstractDataSource,
    fingerprints::Dict{String,Any},
    cached::SourceScan,
)::Bool
    cached.source_id == source_id(source) || return false
    cached.source_fingerprint == source_fingerprint(source) || return false
    return cached.source_item_fingerprints == fingerprints
end

"""Interpret one source item into records and keep the original handles for collection folds."""
function interpret_source_item(
    project::Project,
    source::AbstractDataSource,
    source_item::AbstractDataSourceItem,
    metadata::Union{Nothing,Dict{Tuple{Vararg{String}},Dict{Symbol,Any}}},
)::Tuple{Vector{ItemRecord},Dict{String,AbstractDataItem}}
    handles = data_items(project, source, source_item)::Vector{<:AbstractDataItem}
    records = ItemRecord[
        ItemRecord(handle; source_item)
        for handle in handles
    ]
    apply_collection_metadata!(records, metadata)
    by_key = Dict{String,AbstractDataItem}()
    for (record, handle) in zip(records, handles)
        by_key[item_record_key(record)] = handle
    end
    return records, by_key
end

"""
Interpret source items concurrently and stream each successful record batch.
"""
function interpret_source_items(
    project::Project,
    source::AbstractDataSource,
    source_items::Vector{<:AbstractDataSourceItem},
    metadata::Union{Nothing,Dict{Tuple{Vararg{String}},Dict{Symbol,Any}}};
    on_items::Union{Nothing,Function}=nothing,
    on_progress::Union{Nothing,Function}=nothing,
)::NamedTuple
    processed_count = Base.Threads.Atomic{Int}(0)
    item_count = Base.Threads.Atomic{Int}(0)
    skipped_count = Base.Threads.Atomic{Int}(0)
    callback_lock = ReentrantLock()
    worker_limit = Base.Semaphore(max(1, Base.Threads.nthreads()))
    cancel_requested = get(task_local_storage(), CANCEL_CALLBACK_KEY, nothing)
    failures = ItemFailure[]
    records_by_position = Vector{Vector{ItemRecord}}(undef, length(source_items))
    handles_by_key = Dict{String,AbstractDataItem}()

    @sync for (index, source_item) in pairs(source_items)
        check_cancel()
        Base.acquire(worker_limit)
        Base.Threads.@spawn try
            with_cancel(cancel_requested) do
                check_cancel()
                records, handles = try
                    interpret_source_item(project, source, source_item, metadata)
                catch error
                    is_job_cancelled(error) && rethrow()
                    backtrace = catch_backtrace()
                    @error(
                        "Source item interpretation failed",
                        source=source_label(source),
                        source_item=source_item_id(source_item),
                        exception=(error, backtrace),
                    )
                    lock(callback_lock) do
                        push!(failures, ItemFailure(
                            source_item_id(source_item),
                            "",
                            sprint(showerror, error),
                        ))
                    end
                    ItemRecord[], Dict{String,AbstractDataItem}()
                end
                check_cancel()
                records_by_position[index] = records
                lock(callback_lock) do
                    merge!(handles_by_key, handles)
                    isempty(records) || on_items === nothing || on_items(records)
                end
                isempty(records) ?
                    Base.Threads.atomic_add!(skipped_count, 1) :
                    Base.Threads.atomic_add!(item_count, length(records))
                processed = Base.Threads.atomic_add!(processed_count, 1) + 1
                on_progress === nothing || lock(callback_lock) do
                    on_progress((
                        total_source_items=length(source_items),
                        processed_source_items=processed,
                        loaded_items=item_count[],
                        skipped_source_items=skipped_count[],
                        current_source_item=source_item_id(source_item),
                    ))
                end
            end
        finally
            Base.release(worker_limit)
        end
    end

    records = reduce(append!, records_by_position; init=ItemRecord[])
    return (
        processed_source_items=processed_count[],
        loaded_items=item_count[],
        skipped_source_items=skipped_count[],
        total_source_items=length(source_items),
        failures,
        records,
        handles_by_key,
    )
end

"""Apply collection-node stats supplied by the source."""
function add_collection_stats!(
    project::Project,
    source::AbstractDataSource,
    hierarchy::Hierarchy,
    handles_by_key::Dict{String,AbstractDataItem},
)::Vector{ItemFailure}
    failures = ItemFailure[]
    for (path, node) in hierarchy.index
        isempty(node.items) && continue
        handles = AbstractDataItem[
            handles_by_key[item_record_key(record)]
            for record in node.items
        ]
        try
            merge!(node.stats, collection_stats(project, source, collect(path), handles))
        catch error
            first_item = first(node.items)
            push!(failures, ItemFailure(
                first_item.source_item_id,
                item_record_key(first_item),
                "collection_stats: " * sprint(showerror, error),
            ))
        end
    end
    return failures
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
    count_first::Bool=false,
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
    discovered = source_items(source)
    fingerprints = source_item_fingerprints(discovered)
    check_cancel()
    if cached_source !== nothing && source_unchanged(source, fingerprints, cached_source)
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
    summary = interpret_source_items(
        project,
        source,
        discovered,
        source_collection_metadata(source);
        on_items,
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
    hierarchy = Hierarchy(
        summary.records,
        source_id(source),
        source_collection_metadata(source) !== nothing,
        summary.skipped_source_items,
    )
    emit_progress(
        on_progress;
        phase=:analyzing,
        total_source_items=length(discovered),
        processed_source_items=length(discovered),
        loaded_items=length(summary.records),
        skipped_source_items=summary.skipped_source_items,
    )
    append!(summary.failures, add_collection_stats!(project, source, hierarchy, summary.handles_by_key))
    append!(summary.failures, source_analysis_failures(project, source))
    return SourceScan(
        source_id(source),
        source_label(source),
        source_fingerprint(source),
        fingerprints,
        hierarchy,
        summary.failures,
    )
end
