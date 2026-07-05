"""
One completed source-item pass.

`records` are retained by the index. `interpreted_items` carry effective item data only on direct
interpretation paths; workspace workers put that data in the memory cache before publishing a
lightweight completion.
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

"""
Interpret every logical data item produced by one source item.

The source item is read only by `data_items`. Each result is normalized to a package `DataItem` with
effective metadata. Processing and analysis belong to the workspace work graph.
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
        record = ItemRecord(handle; source_item)
        records[index] = record
        interpreted_items[index] = DataItem(
            ItemRecord(record; metadata=_effective_metadata(
                source, record.collection, record.metadata)),
            item_data(handle),
        )
    end
    finish_source_profile!(
        project, source_item_id_value, source_item_label_value, source_item_path_value,
        source_kind, item_count,
        (time_ns() - source_started) / 1e9, Set([Base.Threads.threadid()]))
    return SourceItemInterpretation(records, interpreted_items, ItemFailure[])
end

