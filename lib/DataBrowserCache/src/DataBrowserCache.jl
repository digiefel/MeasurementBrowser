"""DuckDB cache buffers, build metrics, and project cache domain."""
module DataBrowserCache

import DataBrowserAPI
import DataBrowserAPI: label
using DataBrowserAPI:
    AbstractCollection,
    AbstractDataSourceItem,
    MetadataDict,
    PipelineStage,
    SOURCE_READ,
    SOURCE_INTERPRET,
    ITEM_PROCESS,
    ITEM_ANALYZE,
    COLLECTION_PROCESS,
    COLLECTION_ANALYZE,
    fingerprint,
    id,
    metadata,
    source_id,
    source_item_path,
    source_item_timestamp,
    source_label
using DataBrowserAPI: @timed_dbg

using DuckDB
using DBInterface
import Tables
using Serialization
using Dates

import DataBrowserAPI.ItemIndex:
    CollectionIndex,
    CollectionRecord,
    ItemFailure,
    ItemRecord,
    MetadataValue,
    SourceScan,
    append_item!,
    collection_path_keys,
    emit_progress,
    metadata_dict,
    register_collection!

include("build_metrics.jl")
include("cache_buffer.jl")
include("project_cache_domain.jl")
include("source_stages.jl")

export AbstractCacheDB,
    BuildMetrics,
    CacheResultKey,
    CacheResultStatus,
    CacheDB,
    CacheStageSummary,
    ProjectCacheError,
    ProjectCacheSchemaError,
    ProjectCacheDataError,
    ProjectCacheIdentity,
    ProjectCacheIndex,
    ProjectCacheStatus,
    RESULT_FAILED,
    RESULT_READY,
    cache_stage_summary,
    cache_has_pending_writes,
    cache_pending_counts,
    cached_result_state,
    cached_source_fingerprints,
    clear_cache_index!,
    clear_cached_result_state!,
    close_cache_db!,
    delete_collection_records!,
    delete_collection_metadata!,
    delete_source_item!,
    delete_source_output!,
    has_payload,
    load_cache_index,
    open_memory_cache_db,
    open_cache_db,
    project_cache_identity,
    query_items,
    read_payload,
    record_cache_phase!,
    reset_build_metrics!,
    set_cache_memory_limit!,
    store_collection_metadata!,
    store_collection_index!,
    store_collection_process_result!,
    store_interpreted!,
    store_interpreted_data!,
    store_item_metadata!,
    store_item_metadata_layer!,
    store_processed!,
    store_result_failure!,
    store_source_identity!,
    store_source_read!,
    store_source_result!,
    source_complete,
    cache_knows_source,
    edit_source_item_metadata!,
    source_item_key!,
    source_item_key,
    source_item_id

end
