"""DuckDB cache buffers, build metrics, and project cache domain."""
module DataBrowserCache

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
using Serialization
using Dates

import DataBrowserAPI.ItemIndex:
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

include("build_metrics.jl")
include("cache_buffer.jl")
include("project_cache_domain.jl")

export AbstractCacheDB,
    BuildMetrics,
    CacheResultKey,
    CacheResultKind,
    CacheResultStatus,
    CacheDB,
    CacheStageSummary,
    COLLECTION_ANALYSIS_RESULT,
    COLLECTION_PROCESS_RESULT,
    ITEM_ANALYSIS_RESULT,
    PROCESSING_RESULT,
    ProjectCacheSchemaError,
    ProjectCacheIdentity,
    ProjectCacheStatus,
    RESULT_FAILED,
    RESULT_READY,
    cache_stage_summary,
    cache_built,
    cache_has_pending_writes,
    cache_pending_counts,
    clear_cache_index!,
    clear_cached_result_state!,
    clear_cached_source_state!,
    close_cache_db!,
    delete_collection_metadata!,
    delete_source_item!,
    load_cache_index,
    open_memory_cache_db,
    open_cache_db,
    project_cache_identity,
    query_items,
    read_item_data,
    record_cache_phase!,
    reset_build_metrics!,
    set_cache_memory_limit!,
    start_cache!,
    stop_cache!,
    store_collection_metadata!,
    store_collection_process_result!,
    store_interpreted!,
    store_item_metadata!,
    store_processed!,
    store_result_failure!,
    store_source_item_failure!,
    edit_source_collection_metadata!,
    edit_source_item_metadata!,
    wait_condition_deadline,
    write_meta_header!

end
