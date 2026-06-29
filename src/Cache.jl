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
using ..Profiling: @profile_span, ProfileAttributes

using DuckDB
using DBInterface
using DataFrames: AbstractDataFrame, DataFrame, names, nrow
using SHA
using Serialization
using Dates

# CacheBuffer is the lowest level (raw per-table DuckDB access); it names no higher-level type, so it
# is included first. ProjectCache wraps it: its CacheDB holds a CacheBuffer and drives it.
include("Cache/BuildMetrics.jl")
include("Cache/CacheBuffer.jl")
include("Cache/ProjectCache.jl")
include("Cache/ItemDataCache.jl")

end
