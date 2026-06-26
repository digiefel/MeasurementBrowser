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

"""Return `count` scalar SQL placeholders for a dynamically sized `IN` clause."""
function _sql_placeholders(count::Integer)::String
    count > 0 || throw(ArgumentError("SQL placeholder count must be positive, got $count"))
    return join(fill("?", Int(count)), ", ")
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

"""
Mint one fresh integer surrogate per item id from the cache sequence.

The payload tables and catalog both store this integer (`__mb_seq` / `seq`) in place of the long item
id. A single `range` query draws the whole batch in one round trip.
"""
function _mint_item_seqs!(connection, ids::Vector{String})::Dict{String,Int64}
    seqs = Dict{String,Int64}()
    isempty(ids) && return seqs
    sizehint!(seqs, length(ids))
    index = 0
    for row in DBInterface.execute(
        connection, "SELECT nextval('item_data_seq') AS seq FROM range($(length(ids)))")
        index += 1
        seqs[ids[index]] = Int64(row.seq)
    end
    index == length(ids) || error(
        "Minted $index surrogate keys for $(length(ids)) items",
    )
    return seqs
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
        connection, "SELECT storage_id, seq FROM item_data WHERE item_id = ? AND stage = ?")
    storage_id = nothing
    seq = nothing
    for row in DBInterface.execute(statement, (String(item_id), String(stage)))
        storage_id = String(row.storage_id)
        seq = Int64(row.seq)
    end
    if storage_id !== nothing
        table = _quote_identifier(_dataframe_table_name(storage_id))
        DBInterface.execute(
            DBInterface.prepare(connection, "DELETE FROM $table WHERE __mb_seq = ?"),
            (seq,),
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

"""
Write one worker-sized schema group through a single registered DataFrame.

A one-item group is a no-copy wrapper around the original columns. Larger groups are concatenated
once, bounded to at most one source item per worker, so DuckDB receives one contiguous columnar scan.
"""
function _write_dataframe_group!(
    connection,
    profiler::Profiling.ProfileSession,
    entries::Vector{DataFrameCacheEntry},
    schemas::Dict{String,Vector{String}},
    seqs::Dict{String,Int64},
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
        ids = unique!(String[entry.record.id for entry in entries])
        placeholders = _sql_placeholders(length(ids))
        parameters = tuple(storage_id, ids...)
        old_seqs = Int64[
            Int64(row.seq)
            for row in DBInterface.execute(
                DBInterface.prepare(connection, """
                    SELECT seq FROM item_data
                    WHERE storage_id = ? AND item_id IN ($placeholders)
                """),
                parameters,
            )
        ]
        if !isempty(old_seqs)
            seq_placeholders = _sql_placeholders(length(old_seqs))
            @profile_span profiler :cache :delete_payload_rows ProfileAttributes(
                rows=Int64(sum(entry -> nrow(entry.data), entries)),
                batch_size=Int64(length(entries)),
            ) begin
                DBInterface.execute(
                    DBInterface.prepare(
                        connection, "DELETE FROM $table WHERE __mb_seq IN ($seq_placeholders)"),
                    Tuple(old_seqs),
                )
            end
            DBInterface.execute(
                DBInterface.prepare(connection, """
                    DELETE FROM item_data
                    WHERE storage_id = ? AND item_id IN ($placeholders)
                """),
                parameters,
            )
        end
    end

    seq_column = _temporary_column_name(entries, "__mb_seq")
    row_column = _temporary_column_name(entries, "__mb_row")
    staged = @profile_span profiler :cache :stage_dataframe ProfileAttributes(
        items=Int64(length(entries)),
        rows=Int64(sum(entry -> nrow(entry.data), entries)),
        batch_size=Int64(length(entries)),
    ) begin
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
            seq_values = Vector{Int64}(undef, row_count)
            row_numbers = Vector{Int64}(undef, row_count)
            offset = 0
            for entry in entries
                count = nrow(entry.data)
                rows = (offset + 1):(offset + count)
                fill!(view(seq_values, rows), seqs[entry.record.id])
                row_numbers[rows] .= 1:count
                offset += count
            end
            frame[!, Symbol(seq_column)] = seq_values
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
    seq_expression = length(entries) == 1 ? "?::BIGINT AS __mb_seq" :
        "$(_quote_identifier(seq_column))::BIGINT AS __mb_seq"
    operation = exists ? "INSERT INTO $table" : "CREATE TABLE $table AS"
    view_name = ITEM_DATA_VIEW
    registered = false
    statement = nothing
    @profile_span profiler :cache :write_dataframe ProfileAttributes(
        items=Int64(length(entries)),
        rows=Int64(sum(entry -> nrow(entry.data), entries)),
        batch_size=Int64(length(entries)),
    ) begin
        try
            DuckDB.register_data_frame(connection, staged, view_name)
            registered = true
            statement = DBInterface.prepare(connection, """
                $operation
                SELECT
                    $seq_expression,
                    $(_quote_identifier(row_column))::BIGINT AS __mb_row,
                    $source_columns
                FROM $(_quote_identifier(view_name))
            """)
            if length(entries) == 1
                DBInterface.execute(statement, (seqs[first_entry.record.id],))
            else
                DBInterface.execute(statement)
            end
        finally
            statement === nothing || DBInterface.close!(statement)
            # An aborted transaction may reject cleanup; preserve the original failure.
            if registered
                try
                    DuckDB.unregister_data_frame(connection, view_name)
                catch
                end
            end
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
    profiler::Profiling.ProfileSession,
    items::Vector{ItemRecord},
    data::Vector,
    stage::Symbol,
)::Nothing
    stage in (:interpreted, :processed) ||
        throw(ArgumentError("unknown item-data cache stage '$stage'"))
    schemas = _load_dataframe_schemas(connection)
    groups = Dict{String,Vector{DataFrameCacheEntry}}()
    @profile_span profiler :cache :prepare_item_data ProfileAttributes(
        stage=stage,
        items=Int64(length(items)),
    ) begin
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
    seqs = _mint_item_seqs!(
        connection, String[entry.record.id for entries in values(groups) for entry in entries])
    @profile_span profiler :cache :write_dataframes ProfileAttributes(
        stage=stage,
        items=Int64(length(items)),
        batch_size=Int64(length(items)),
    ) begin
        for entries in values(groups), batch in _dataframe_write_batches(entries)
            _write_dataframe_group!(connection, profiler, batch, schemas, seqs)
        end
    end

    @profile_span profiler :cache :index_item_data ProfileAttributes(
        stage=stage,
        items=Int64(length(items)),
    ) begin
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
                    seqs[record.id],
                )
                    DuckDB.append(appender, value)
                end
                DuckDB.end_row(appender)
            end
            DuckDB.flush(appender)
        finally
            # Closing flushes; on an already-aborted transaction that itself throws and would mask the
            # real error (see _write_dataframe_group!). The rollback discards the appended rows anyway.
            try
                DuckDB.close(appender)
            catch
            end
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

"""Read valid cached DataFrames with scalar predicates so DuckDB can prune clustered payloads."""
function _read_cached_item_data(
    connection,
    profiler::Profiling.ProfileSession,
    items::Vector{ItemRecord},
    stage::Symbol,
)::Vector{Any}
    stage in (:interpreted, :processed) ||
        throw(ArgumentError("unknown item-data cache stage '$stage'"))
    loaded = Any[nothing for _ in items]
    requested = Dict{String,Vector{Int}}()
    ids = String[]
    for (position, item) in pairs(items)
        check_cancel()
        _item_data_cacheable(item) || continue
        positions = get!(requested, item.id) do
            push!(ids, item.id)
            Int[]
        end
        push!(positions, position)
    end
    isempty(ids) && return loaded

    entries_type = Tuple{Int64,Vector{Int},Int,String}
    by_storage = Dict{String,Tuple{Vector{String},Vector{entries_type}}}()
    placeholders = _sql_placeholders(length(ids))
    meta_statement = DBInterface.prepare(connection, """
        SELECT d.item_id, d.sif_hash, d.if_hash, d.storage_id, d.row_count, d.seq,
               s.column_names
        FROM item_data d
        JOIN dataframe_schemas s ON s.storage_id = d.storage_id
        WHERE d.stage = ? AND d.item_id IN ($placeholders)
    """)
    for row in DBInterface.execute(meta_statement, tuple(String(stage), ids...))
        id = String(row.item_id)
        positions = Int[
            position
            for position in get(requested, id, Int[])
            if _hash_matches(
                row.sif_hash,
                _fingerprint_hash(items[position].source_item_fingerprint),
            ) && _hash_matches(
                row.if_hash,
                _fingerprint_hash(items[position].item_fingerprint),
            )
        ]
        isempty(positions) && continue
        seq = Int64(row.seq)
        storage_id = String(row.storage_id)
        column_names = String[String(name) for name in row.column_names]
        schema, entries = get!(by_storage, storage_id) do
            (column_names, entries_type[])
        end
        schema == column_names || error(
            "Cached DataFrame schema '$storage_id' has inconsistent column names",
        )
        push!(entries, (seq, positions, Int(row.row_count), id))
    end
    isempty(by_storage) && return loaded

    for (storage_id, (column_names, entries)) in by_storage
        table = _quote_identifier(_dataframe_table_name(storage_id))
        selected_columns = join(
            (_quote_identifier("c$index") for index in eachindex(column_names)), ", ")
        placeholders = _sql_placeholders(length(entries))
        statement = DBInterface.prepare(connection, """
            SELECT __mb_seq, $selected_columns
            FROM $table
            WHERE __mb_seq IN ($placeholders)
            ORDER BY __mb_seq, __mb_row
        """)
        seqs = Tuple(entry[1] for entry in entries)
        raw = @profile_span profiler :cache :payload_query ProfileAttributes(
            stage=stage,
            items=Int64(length(entries)),
        ) begin
            DataFrame(DBInterface.execute(statement, seqs))
        end
        expected = Dict(entry[1] => (entry[2], entry[3], entry[4]) for entry in entries)

        function reconstruct(range)::DataFrame
            columns = [
                Symbol(name) => view(raw[!, Symbol("c$index")], range)
                for (index, name) in pairs(column_names)
            ]
            return DataFrame(columns; copycols=false)
        end

        seq_column = raw[!, :__mb_seq]
        seen = Set{Int64}()
        first_row = 1
        total = nrow(raw)
        while first_row <= total
            check_cancel()
            seq = Int64(seq_column[first_row])
            last_row = first_row
            while last_row < total && Int64(seq_column[last_row + 1]) == seq
                last_row += 1
            end
            entry = get(expected, seq, nothing)
            entry === nothing && error(
                "Cached DataFrame read returned unexpected sequence $seq from '$storage_id'",
            )
            positions, expected_rows, id = entry
            (last_row - first_row + 1) == expected_rows || error(
                "Cached DataFrame for item '$id' has " *
                "$(last_row - first_row + 1) rows; expected $expected_rows",
            )
            data = reconstruct(first_row:last_row)
            for position in positions
                loaded[position] = DataItem(items[position], data)
            end
            push!(seen, seq)
            first_row = last_row + 1
        end
        # Zero-row items never appear in the data table; rebuild them from the typed empty columns.
        for (seq, (positions, expected_rows, id)) in expected
            seq in seen && continue
            expected_rows == 0 || error(
                "Cached DataFrame for item '$id' is missing its $expected_rows rows",
            )
            data = reconstruct(1:0)
            for position in positions
                loaded[position] = DataItem(items[position], data)
            end
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
    return @profile_span cachedb.profiler :cache :read_item_data ProfileAttributes(
        stage=stage,
        items=Int64(length(items)),
    ) begin
        with_reader_snapshot(cachedb) do connection
            _read_cached_item_data(connection, cachedb.profiler, items, stage)
        end
    end
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
    @profile_span cachedb.profiler :cache :write_item_data ProfileAttributes(
        stage=stage,
        items=Int64(length(items)),
        batch_size=Int64(length(items)),
    ) begin
        with_writer_transaction(cachedb) do connection
            _write_cached_item_data!(connection, cachedb.profiler, items, data, stage)
        end
    end
    return nothing
end
