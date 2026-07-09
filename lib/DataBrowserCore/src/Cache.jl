"""DuckDB cache buffers, build metrics, and project cache domain."""
module Cache

import DataBrowserAPI
using DataBrowserAPI:
    AbstractDataItem,
    AbstractDataSourceItem,
    MetadataDict,
    cacheable,
    fingerprint,
    item_data,
    metadata,
    source_id,
    source_item_id,
    source_item_path,
    source_item_timestamp,
    source_label
import DataBrowserProfiling as Profiling
using DataBrowserProfiling: @profile_span, ProfileAttributes

using DuckDB
using DBInterface
using DataFrames: AbstractDataFrame, DataFrame, groupby, names, nrow
using SHA
using Serialization
using Dates

import ..ItemIndex:
    DataItem,
    Hierarchy,
    ItemFailure,
    ItemRecord,
    MetadataValue,
    SourceScan,
    all_items,
    collection_path_tuple,
    emit_progress,
    insert_item!,
    metadata_dict

const _CACHE_SRC = joinpath(@__DIR__, "..", "..", "DataBrowserCache", "src")

include(joinpath(_CACHE_SRC, "build_metrics.jl"))
include(joinpath(_CACHE_SRC, "cache_buffer.jl"))
include(joinpath(_CACHE_SRC, "project_cache_domain.jl"))

end
