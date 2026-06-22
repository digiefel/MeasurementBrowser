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

"""One standard cacheable item prepared for a native DataFrame write."""
struct DataFrameCacheEntry
    record::ItemRecord
    data::AbstractDataFrame
    storage_id::String
    column_names::Vector{String}
end

"""Choose a temporary column name absent from every DataFrame in a write group."""
function _temporary_column_name(
    entries::Vector{DataFrameCacheEntry},
    base::String,
)::String
    candidate = base
    suffix = 0
    while any(entry -> candidate in entry.column_names, entries)
        suffix += 1
        candidate = "$(base)_$suffix"
    end
    return candidate
end

"""Delete one item's physical rows and item-data index row for one cache stage."""
function _delete_cached_item_data!(
    connection,
    item_id::AbstractString,
    stage::Symbol,
)::Nothing
    statement = DBInterface.prepare(
        connection, "SELECT storage_id FROM item_data WHERE item_id = ? AND stage = ?")
    storage_id = nothing
    for row in DBInterface.execute(statement, (String(item_id), String(stage)))
        storage_id = String(row.storage_id)
    end
    if storage_id !== nothing
        table = _quote_identifier(_dataframe_table_name(storage_id))
        DBInterface.execute(
            DBInterface.prepare(connection, "DELETE FROM $table WHERE __mb_item_id = ?"),
            (String(item_id),),
        )
        DBInterface.execute(
            DBInterface.prepare(
                connection, "DELETE FROM item_data WHERE item_id = ? AND stage = ?"),
            (String(item_id), String(stage)),
        )
        _drop_unused_dataframe_storage!(connection, storage_id)
    end
    return nothing
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
                WHERE __mb_item_id IN (
                    SELECT item_id FROM item_data
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

"""Prepare one standard DataItem carrying a DataFrame or DataFrame view for native storage."""
function _dataframe_cache_entry(
    record::ItemRecord,
    item::DataItem,
    stage::Symbol,
)::DataFrameCacheEntry
    data = item_data(item)
    data isa AbstractDataFrame || throw(UnsupportedItemDataCacheError(
        record.id, typeof(item), typeof(data)))
    column_names = names(data)
    isempty(column_names) && throw(ArgumentError(
        "Cannot persist a DataFrame with no columns; set cacheable(item) = false for that item",
    ))
    return DataFrameCacheEntry(
        record,
        data,
        _dataframe_storage_id(data, stage),
        column_names,
    )
end

"""Fail clearly for custom cacheable item types until they provide a native cache implementation."""
function _dataframe_cache_entry(
    record::ItemRecord,
    item::AbstractDataItem,
    ::Symbol,
)::DataFrameCacheEntry
    data = item_data(item)
    throw(UnsupportedItemDataCacheError(record.id, typeof(item), typeof(data)))
end

"""Fill the temporary table identifying item rows replaced by one DataFrame write."""
function _fill_item_data_write_group!(
    connection,
    entries::Vector{DataFrameCacheEntry},
)::Nothing
    DBInterface.execute(connection, """
        CREATE TEMP TABLE IF NOT EXISTS writing_item_data(item_id TEXT)
    """)
    DBInterface.execute(connection, "DELETE FROM writing_item_data")
    appender = DuckDB.Appender(connection, "writing_item_data")
    try
        for entry in entries
            DuckDB.append(appender, entry.record.id)
            DuckDB.end_row(appender)
        end
        DuckDB.flush(appender)
    finally
        DuckDB.close(appender)
    end
    return nothing
end

"""
Write one worker-sized schema group through a single registered DataFrame.

A one-item group is a no-copy wrapper around the original columns. Larger groups are concatenated
once, bounded to at most one source item per worker, so DuckDB receives one contiguous columnar scan.
"""
function _write_dataframe_group!(
    connection,
    entries::Vector{DataFrameCacheEntry},
    schemas::Dict{String,Vector{String}},
)::Nothing
    first_entry = first(entries)
    storage_id = first_entry.storage_id
    column_names = first_entry.column_names
    all(entry -> entry.storage_id == storage_id, entries) ||
        error("DataFrame cache write group contains more than one storage id")
    all(entry -> entry.column_names == column_names, entries) ||
        error("DataFrame cache write group contains inconsistent column names")

    table = _quote_identifier(_dataframe_table_name(storage_id))
    exists = haskey(schemas, storage_id)
    if exists
        schemas[storage_id] == column_names ||
            error("DataFrame cache schema $storage_id has inconsistent column names")
        _fill_item_data_write_group!(connection, entries)
        DBInterface.execute(connection, """
            DELETE FROM $table
            WHERE __mb_item_id IN (SELECT item_id FROM writing_item_data)
        """)
        DBInterface.execute(
            DBInterface.prepare(connection, """
                DELETE FROM item_data
                WHERE storage_id = ?
                  AND item_id IN (SELECT item_id FROM writing_item_data)
            """),
            (storage_id,),
        )
    end

    item_column = _temporary_column_name(entries, "__mb_item_id")
    row_column = _temporary_column_name(entries, "__mb_row")
    staged = @timeit_debug TIMER "cache/stage_dataframe" begin
        if length(entries) == 1
            frame = DataFrame(first_entry.data; copycols=false)
            frame[!, Symbol(row_column)] = Base.OneTo(nrow(frame))
            frame
        else
            row_count = sum(entry -> nrow(entry.data), entries)
            staged_columns = Vector{AbstractVector}(undef, length(column_names))
            for (column_index, name) in pairs(column_names)
                column = Vector{eltype(first_entry.data[!, name])}(undef, row_count)
                offset = 0
                for entry in entries
                    source = entry.data[!, name]
                    copyto!(column, offset + 1, source, 1, length(source))
                    offset += length(source)
                end
                staged_columns[column_index] = column
            end
            frame = DataFrame(staged_columns, Symbol.(column_names); copycols=false)
            item_ids = Vector{String}(undef, row_count)
            row_numbers = Vector{Int64}(undef, row_count)
            offset = 0
            for entry in entries
                count = nrow(entry.data)
                rows = (offset + 1):(offset + count)
                fill!(view(item_ids, rows), entry.record.id)
                row_numbers[rows] .= 1:count
                offset += count
            end
            frame[!, Symbol(item_column)] = item_ids
            frame[!, Symbol(row_column)] = row_numbers
            frame
        end
    end

    source_columns = join(
        (
            "$(_quote_identifier(name)) AS $(_quote_identifier("c$index"))"
            for (index, name) in pairs(column_names)
        ),
        ", ",
    )
    item_expression = length(entries) == 1 ? "? AS __mb_item_id" :
        "$(_quote_identifier(item_column)) AS __mb_item_id"
    operation = exists ? "INSERT INTO $table" : "CREATE TABLE $table AS"
    registered = false
    @timeit_debug TIMER "cache/write_dataframe" begin
        try
            DuckDB.register_data_frame(connection, staged, ITEM_DATA_VIEW)
            registered = true
            statement = DBInterface.prepare(connection, """
                $operation
                SELECT
                    $item_expression,
                    $(_quote_identifier(row_column))::BIGINT AS __mb_row,
                    $source_columns
                FROM $(_quote_identifier(ITEM_DATA_VIEW))
            """)
            if length(entries) == 1
                DBInterface.execute(statement, (first_entry.record.id,))
            else
                DBInterface.execute(statement)
            end
        finally
            registered && DuckDB.unregister_data_frame(connection, ITEM_DATA_VIEW)
        end
    end
    if !exists
        DBInterface.execute(
            DBInterface.prepare(connection, "INSERT INTO dataframe_schemas VALUES (?, ?)"),
            (storage_id, column_names),
        )
        schemas[storage_id] = column_names
    end
    return nothing
end

"""Partition one schema by source item without ever splitting an expanded source."""
function _dataframe_write_batches(
    entries::Vector{DataFrameCacheEntry},
)::Vector{Vector{DataFrameCacheEntry}}
    by_source = Dict{String,Vector{DataFrameCacheEntry}}()
    for entry in entries
        push!(
            get!(() -> DataFrameCacheEntry[], by_source, entry.record.source_item_id),
            entry,
        )
    end
    source_groups = collect(values(by_source))
    groups_per_batch = max(1, Base.Threads.nthreads())
    batches = Vector{Vector{DataFrameCacheEntry}}()
    for group_range in Iterators.partition(eachindex(source_groups), groups_per_batch)
        batch = DataFrameCacheEntry[]
        for index in group_range
            append!(batch, source_groups[index])
        end
        push!(batches, batch)
    end
    return batches
end

"""Write aligned loaded items on one connection, removing entries that are no longer cacheable."""
function _write_cached_item_data!(
    connection,
    items::Vector{ItemRecord},
    data::Vector,
    stage::Symbol,
)::Nothing
    stage in (:interpreted, :processed) ||
        throw(ArgumentError("unknown item-data cache stage '$stage'"))
    schemas = _load_dataframe_schemas(connection)
    groups = Dict{String,Vector{DataFrameCacheEntry}}()
    @timeit_debug TIMER "cache/prepare_item_data" begin
        for (record, value) in zip(items, data)
            check_cancel()
            if !_item_data_cacheable(record) || value === nothing
                _delete_cached_item_data!(connection, record.id, stage)
                continue
            end
            value isa AbstractDataItem || throw(UnsupportedItemDataCacheError(
                record.id, typeof(value), typeof(value)))
            if !cacheable(value)
                _delete_cached_item_data!(connection, record.id, stage)
                continue
            end
            entry = _dataframe_cache_entry(record, value, stage)
            push!(get!(() -> DataFrameCacheEntry[], groups, entry.storage_id), entry)
        end
    end
    @timeit_debug TIMER "cache/write_dataframes" begin
        for entries in values(groups), batch in _dataframe_write_batches(entries)
            _write_dataframe_group!(connection, batch, schemas)
        end
    end

    @timeit_debug TIMER "cache/index_item_data" begin
        appender = DuckDB.Appender(connection, "item_data")
        try
            for entries in values(groups), entry in entries
                record = entry.record
                for value in (
                    record.id,
                    String(stage),
                    record.source_item_id,
                    _fingerprint_hash(record.source_item_fingerprint),
                    _fingerprint_hash(record.item_fingerprint),
                    entry.storage_id,
                    Int64(nrow(entry.data)),
                )
                    DuckDB.append(appender, value)
                end
                DuckDB.end_row(appender)
            end
            DuckDB.flush(appender)
        finally
            DuckDB.close(appender)
        end
    end
    return nothing
end

"""Whether a stored fingerprint hash matches the requested one, treating SQL NULL as `nothing`."""
function _hash_matches(stored, requested)::Bool
    stored = _null_to_nothing(stored)
    stored === nothing && return requested === nothing
    return requested !== nothing && String(stored) == requested
end

"""
A small selection — an interactive plot of a handful of items — is served by the temp-table-free
id-list path; a bulk materialization (collection analysis, a whole-tree reload) is served by the
position-ordered temp-table join, which sorts by an integer instead of by item id and slices rows
contiguously. Each path wins on its own regime.
"""
const SMALL_ITEM_DATA_READ = 16

"""Read valid native item data, choosing the read strategy that fits the selection size."""
function _read_cached_item_data(
    connection,
    items::Vector{ItemRecord},
    stage::Symbol,
)::Vector{Any}
    stage in (:interpreted, :processed) ||
        throw(ArgumentError("unknown item-data cache stage '$stage'"))
    return count(_item_data_cacheable, items) <= SMALL_ITEM_DATA_READ ?
        _read_cached_item_data_small(connection, items, stage) :
        _read_cached_item_data_bulk(connection, items, stage)
end

"""Fill the temporary request table with item-data cache keys and their original positions."""
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

"""Bulk read: join a temp request table so rows come back in integer-position order for slicing."""
function _read_cached_item_data_bulk(
    connection,
    items::Vector{ItemRecord},
    stage::Symbol,
)::Vector{Any}
    loaded = Any[nothing for _ in items]
    _fill_item_data_request!(connection, items)
    grouped = Dict{String,Vector{Tuple{Int,Int}}}()
    metadata_statement = DBInterface.prepare(connection, """
        SELECT r.position, d.storage_id, d.row_count
        FROM requested_item_data r
        JOIN item_data d ON d.item_id = r.item_id
        WHERE d.stage = ?
          AND d.sif_hash = r.sif_hash
          AND (
              (d.if_hash IS NULL AND r.if_hash IS NULL)
              OR d.if_hash = r.if_hash
          )
        ORDER BY r.position
    """)
    for row in DBInterface.execute(metadata_statement, (String(stage),))
        position = Int(row.position)
        1 <= position <= length(items) || error(
            "Invalid item-data cache read position $position for $(length(items)) requested items",
        )
        push!(
            get!(() -> Tuple{Int,Int}[], grouped, String(row.storage_id)),
            (position, Int(row.row_count)),
        )
    end

    schemas = _load_dataframe_schemas(connection)
    for (storage_id, entries) in grouped
        column_names = get(schemas, storage_id, nothing)
        column_names === nothing && error(
            "Cached item data refers to missing DataFrame schema '$storage_id'",
        )
        table = _quote_identifier(_dataframe_table_name(storage_id))
        selected_columns = join(
            (
                "d.$(_quote_identifier("c$index"))"
                for index in eachindex(column_names)
            ),
            ", ",
        )
        statement = DBInterface.prepare(connection, """
            SELECT r.position, $selected_columns
            FROM requested_item_data r
            JOIN item_data i ON i.item_id = r.item_id
            JOIN $table d ON d.__mb_item_id = r.item_id
            WHERE i.storage_id = ?
              AND i.stage = ?
              AND i.sif_hash = r.sif_hash
              AND (
                  (i.if_hash IS NULL AND r.if_hash IS NULL)
                  OR i.if_hash = r.if_hash
              )
            ORDER BY r.position, d.__mb_row
        """)
        raw = DataFrame(DBInterface.execute(statement, (storage_id, String(stage))))
        row_start = 1
        for (position, expected_rows) in entries
            check_cancel()
            if expected_rows > 0
                row_start <= nrow(raw) || error(
                    "Cached DataFrame for item '$(items[position].id)' is missing its rows",
                )
                Int(raw.position[row_start]) == position || error(
                    "Cached DataFrame rows are out of order for item '$(items[position].id)'",
                )
            end
            row_stop = row_start + expected_rows - 1
            row_stop <= nrow(raw) || error(
                "Cached DataFrame for item '$(items[position].id)' has fewer than " *
                "$expected_rows rows",
            )
            columns = [
                Symbol(name) => raw[row_start:row_stop, Symbol("c$index")]
                for (index, name) in pairs(column_names)
            ]
            loaded[position] = DataItem(
                items[position],
                DataFrame(columns; copycols=false),
            )
            row_start = row_stop + 1
        end
        row_start == nrow(raw) + 1 || error(
            "Cached DataFrame read returned $(nrow(raw) - row_start + 1) unexpected rows",
        )
    end
    return loaded
end

"""
Small read: look up ids with bound id lists — no temporary request table — and validate fingerprint
freshness in Julia, so an interactive selection pays two lightweight queries per schema instead of
building, filling, and joining a temp table on every read.
"""
function _read_cached_item_data_small(
    connection,
    items::Vector{ItemRecord},
    stage::Symbol,
)::Vector{Any}
    loaded = Any[nothing for _ in items]
    requested = Dict{String,Tuple{Int,Union{Nothing,String},Union{Nothing,String}}}()
    ids = String[]
    for (position, item) in pairs(items)
        check_cancel()
        _item_data_cacheable(item) || continue
        haskey(requested, item.id) && continue
        requested[item.id] = (
            position,
            _fingerprint_hash(item.source_item_fingerprint),
            _fingerprint_hash(item.item_fingerprint),
        )
        push!(ids, item.id)
    end
    isempty(ids) && return loaded

    # Map each fresh requested item to its storage table; the freshness check is done here in Julia.
    by_storage = Dict{String,Vector{Tuple{String,Int,Int}}}()
    meta_statement = DBInterface.prepare(connection, """
        SELECT item_id, sif_hash, if_hash, storage_id, row_count
        FROM item_data WHERE stage = ? AND item_id IN ?
    """)
    for row in DBInterface.execute(meta_statement, (String(stage), ids))
        id = String(row.item_id)
        entry = get(requested, id, nothing)
        entry === nothing && continue
        position, sif_hash, if_hash = entry
        (_hash_matches(row.sif_hash, sif_hash) && _hash_matches(row.if_hash, if_hash)) || continue
        push!(
            get!(() -> Tuple{String,Int,Int}[], by_storage, String(row.storage_id)),
            (id, position, Int(row.row_count)),
        )
    end
    isempty(by_storage) && return loaded

    schemas = _load_dataframe_schemas(connection)
    for (storage_id, group) in by_storage
        column_names = get(schemas, storage_id, nothing)
        column_names === nothing && error(
            "Cached item data refers to missing DataFrame schema '$storage_id'",
        )
        table = _quote_identifier(_dataframe_table_name(storage_id))
        selected_columns = join(
            (_quote_identifier("c$index") for index in eachindex(column_names)), ", ")
        statement = DBInterface.prepare(connection, """
            SELECT __mb_item_id, $selected_columns
            FROM $table
            WHERE __mb_item_id IN ?
            ORDER BY __mb_item_id, __mb_row
        """)
        raw = DataFrame(DBInterface.execute(statement, (String[entry[1] for entry in group],)))
        expected = Dict(entry[1] => (entry[2], entry[3]) for entry in group)

        reconstruct(range) = DataFrame(
            [Symbol(name) => raw[range, Symbol("c$index")]
             for (index, name) in pairs(column_names)];
            copycols=false,
        )

        identifiers = raw[!, :__mb_item_id]
        seen = Set{String}()
        first_row = 1
        total = nrow(raw)
        while first_row <= total
            check_cancel()
            id = String(identifiers[first_row])
            last_row = first_row
            while last_row < total && String(identifiers[last_row + 1]) == id
                last_row += 1
            end
            position, expected_rows = expected[id]
            (last_row - first_row + 1) == expected_rows || error(
                "Cached DataFrame for item '$id' has $(last_row - first_row + 1) rows; " *
                "expected $expected_rows",
            )
            loaded[position] = DataItem(items[position], reconstruct(first_row:last_row))
            push!(seen, id)
            first_row = last_row + 1
        end
        # Zero-row items never appear in the data table; rebuild them from the typed empty columns.
        for (id, (position, expected_rows)) in expected
            id in seen && continue
            expected_rows == 0 || error(
                "Cached DataFrame for item '$id' is missing its $expected_rows rows",
            )
            loaded[position] = DataItem(items[position], reconstruct(1:0))
        end
    end
    return loaded
end

"""Read valid cached item data for one stage without source access."""
function read_cached_item_data(
    cachedb::CacheDB,
    items::Vector{ItemRecord};
    stage::Symbol=:processed,
)::Vector{Any}
    isempty(items) && return Any[]
    return @timeit_debug TIMER "cache/read_item_data" with_persistent_reader(cachedb) do connection
        _read_cached_item_data(connection, items, stage)
    end
end

"""Persist aligned loaded items for one stage through a workspace's persistent writer."""
function write_cached_item_data!(
    cachedb::CacheDB,
    items::Vector{ItemRecord},
    data::Vector;
    stage::Symbol=:processed,
)::Nothing
    length(items) == length(data) ||
        throw(DimensionMismatch("items and data must have equal lengths"))
    isempty(items) && return nothing
    @timeit_debug TIMER "cache/write_item_data" with_writer(cachedb) do connection
        DBInterface.execute(connection, "BEGIN TRANSACTION")
        try
            _write_cached_item_data!(connection, items, data, stage)
            DBInterface.execute(connection, "COMMIT")
        catch
            DBInterface.execute(connection, "ROLLBACK")
            rethrow()
        end
    end
    return nothing
end
