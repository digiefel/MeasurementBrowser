module Cache

import ..ItemIndex:
    DataItem,
    ItemRecord,
    ItemFailure,
    Hierarchy,
    HierarchyNode,
    SourceScan,
    MetadataValue,
    MetadataDict,
    metadata_dict,
    insert_item!,
    collection_path_key,
    collection_path_tuple,
    check_cancel,
    emit_progress,
    source_parameter_state
import ..Projects:
    AbstractDataItem,
    cacheable,
    item_data,
    source_id,
    source_label
import ..Profiling
using ..Profiling: TIMER
using TimerOutputs: @timeit_debug

include("Cache/ProjectCache.jl")
include("Cache/ItemDataCache.jl")

Profiling.register_instrumented!(@__MODULE__)

end
