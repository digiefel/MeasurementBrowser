module Cache

import ..ItemIndex:
    DataItem,
    ItemRecord,
    ItemFailure,
    Hierarchy,
    SourceScan,
    MetadataValue,
    MetadataDict,
    metadata_dict,
    insert_item!,
    all_items,
    collection_path_tuple,
    emit_progress
import ..Projects
import ..Projects:
    AbstractDataItem,
    AbstractDataSourceItem,
    cacheable,
    fingerprint,
    item_data,
    metadata,
    source_id,
    source_item_id,
    source_item_path,
    source_item_timestamp,
    source_label
import ..Profiling
using ..Profiling: @profile_span, ProfileAttributes

using DuckDB
using DBInterface
using DataFrames: AbstractDataFrame, DataFrame, groupby, names, nrow
using SHA
using Serialization
using Dates

# CacheBuffer is the domain-free keyed mechanism and is included first. ProjectCache defines every
# concrete row, SQL operation, and the CacheDB that owns the typed buffers directly.
include("Cache/BuildMetrics.jl")
include("Cache/CacheBuffer.jl")
include("Cache/ProjectCache.jl")

end
