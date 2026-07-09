"""
DuckDB cache buffers and build metrics (standalone slice).

The full project cache loads via [`cache_in_parent.jl`](@ref) inside the umbrella package.
"""
module DataBrowserCache

import DataBrowserAPI
using DataBrowserAPI: AbstractDataItem, AbstractDataSourceItem, MetadataDict, cacheable,
    fingerprint, item_data, metadata, source_id, source_item_id, source_item_path,
    source_item_timestamp, source_label
import DataBrowserProfiling as Profiling

using DuckDB
using DBInterface
using DataFrames: AbstractDataFrame, DataFrame, groupby, names, nrow
using SHA
using Serialization
using Dates

include("build_metrics.jl")
include("cache_buffer.jl")

end
