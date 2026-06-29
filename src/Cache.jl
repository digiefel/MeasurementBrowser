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

include("Cache/BuildMetrics.jl")
include("Cache/ProjectCache.jl")
include("Cache/ItemDataCache.jl")
# CacheBuffer is the lowest level (raw per-table DuckDB access); ProjectCache (above) wraps it. It is
# included after ProjectCache only because its owner struct still names CacheDB — step 2 inverts that.
include("Cache/CacheBuffer.jl")

end
