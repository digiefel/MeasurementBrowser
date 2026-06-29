# --------------------------------------------------------------------------------------------------
# Item data
# --------------------------------------------------------------------------------------------------

# Cached data is keyed by item id and validated against both fingerprints. The source-item
# fingerprint is required; an absent item fingerprint intentionally means source-level invalidation.
_item_data_cacheable(item::ItemRecord)::Bool = item.source_item_fingerprint !== nothing

"""Quote an SQL identifier without treating user-provided column names as SQL."""
function _quote_identifier(identifier::AbstractString)::String
    return "\"" * replace(String(identifier), "\"" => "\"\"") * "\""
end

"""Return a stable id for one stage and a DataFrame's ordered column names and element types."""
function _dataframe_storage_id(data::AbstractDataFrame, stage::Symbol)::String
    schema = IOBuffer()
    print(schema, String(stage), ':')
    for name in names(data)
        type_name = string(eltype(data[!, name]))
        print(schema, ncodeunits(name), ':', name, ncodeunits(type_name), ':', type_name)
    end
    return bytes2hex(sha1(take!(schema)))
end

"""Return the physical table name for a validated DataFrame schema id."""
function _dataframe_table_name(storage_id::AbstractString)::String
    id = String(storage_id)
    occursin(r"^[0-9a-f]{40}$", id) ||
        error("Invalid DataFrame cache storage id '$id'; expected a SHA-1 digest")
    return "dataframe_$id"
end

"""Read the known DataFrame schemas into memory once for a cache operation."""
function _load_dataframe_schemas(connection)::Dict{String,Vector{String}}
    schemas = Dict{String,Vector{String}}()
    for row in DBInterface.execute(
        connection, "SELECT storage_id, column_names FROM dataframe_schemas")
        schemas[String(row.storage_id)] = String[name for name in row.column_names]
    end
    return schemas
end

"""Drop a physical DataFrame table after its final item-data entry has been removed."""
function _drop_unused_dataframe_storage!(
    connection,
    storage_id::AbstractString,
)::Nothing
    statement = DBInterface.prepare(
        connection, "SELECT count(*) AS count FROM item_data WHERE storage_id = ?")
    count = only(DBInterface.execute(statement, (String(storage_id),))).count
    count == 0 || return nothing
    DBInterface.execute(
        connection, "DROP TABLE $(_quote_identifier(_dataframe_table_name(storage_id)))")
    DBInterface.execute(
        DBInterface.prepare(connection, "DELETE FROM dataframe_schemas WHERE storage_id = ?"),
        (String(storage_id),),
    )
    return nothing
end

"""Delete all cached data derived from one source item."""
function _delete_cached_source_item_data!(
    connection,
    source_item_id::AbstractString,
)::Nothing
    storage_ids = String[
        String(row.storage_id)
        for row in DBInterface.execute(
            DBInterface.prepare(connection, """
                SELECT DISTINCT storage_id FROM item_data WHERE source_item_id = ?
            """),
            (String(source_item_id),),
        )
    ]
    for storage_id in storage_ids
        table = _quote_identifier(_dataframe_table_name(storage_id))
        DBInterface.execute(
            DBInterface.prepare(connection, """
                DELETE FROM $table
                WHERE __mb_seq IN (
                    SELECT seq FROM item_data
                    WHERE source_item_id = ? AND storage_id = ?
                )
            """),
            (String(source_item_id), storage_id),
        )
    end
    DBInterface.execute(
        DBInterface.prepare(connection, "DELETE FROM item_data WHERE source_item_id = ?"),
        (String(source_item_id),),
    )
    for storage_id in storage_ids
        _drop_unused_dataframe_storage!(connection, storage_id)
    end
    return nothing
end

"""Delete every cached DataFrame and its schema/index rows."""
function _delete_all_cached_item_data!(connection)::Nothing
    for storage_id in keys(_load_dataframe_schemas(connection))
        DBInterface.execute(
            connection, "DROP TABLE $(_quote_identifier(_dataframe_table_name(storage_id)))")
    end
    DBInterface.execute(connection, "DELETE FROM item_data")
    DBInterface.execute(connection, "DELETE FROM dataframe_schemas")
    return nothing
end

"""Fill the temporary request table used by the large metadata-only readiness query."""
function _fill_item_data_request!(connection, items::Vector{ItemRecord})::Nothing
    DBInterface.execute(connection, """
        CREATE TEMP TABLE IF NOT EXISTS requested_item_data(
            position BIGINT,
            item_id TEXT,
            sif_hash TEXT,
            if_hash TEXT)
    """)
    DBInterface.execute(connection, "DELETE FROM requested_item_data")
    appender = DuckDB.Appender(connection, "requested_item_data")
    try
        for (position, item) in pairs(items)
            check_cancel()
            _item_data_cacheable(item) || continue
            for value in (
                Int64(position),
                item.id,
                _fingerprint_hash(item.source_item_fingerprint),
                _fingerprint_hash(item.item_fingerprint),
            )
                DuckDB.append(appender, value)
            end
            DuckDB.end_row(appender)
        end
        DuckDB.flush(appender)
    finally
        DuckDB.close(appender)
    end
    return nothing
end

"""Return item ids whose requested cache stage exists with matching fingerprints."""
function cached_item_data_ids(
    cachedb::CacheDB,
    items::Vector{ItemRecord};
    stage::Symbol=:processed,
)::Set{String}
    isempty(items) && return Set{String}()
    stage in (:interpreted, :processed) ||
        throw(ArgumentError("unknown item-data cache stage '$stage'"))
    return with_reader(cachedb) do connection
        _fill_item_data_request!(connection, items)
        Set(String(row.item_id) for row in DBInterface.execute(connection, """
            SELECT r.item_id
            FROM requested_item_data r
            JOIN item_data d ON d.item_id = r.item_id
            WHERE d.stage = ?
              AND d.sif_hash = r.sif_hash
              AND (
                  (d.if_hash IS NULL AND r.if_hash IS NULL)
                  OR d.if_hash = r.if_hash
              )
        """, (String(stage),)))
    end
end
