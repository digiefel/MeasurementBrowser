# Does DuckDB actually let multiple connections write the same table in parallel, or does it
# serialize internally? This mimics our cache write pattern: one shared, unindexed table; each "item"
# runs BEGIN / DELETE its own rows / INSERT K rows / COMMIT on its own connection. We compare one
# connection doing all items serially against N connections doing them concurrently.

using DuckDB, DBInterface, Printf

const N_ITEMS = 400
const ROWS_PER_ITEM = 2_000
const PREFILL_ITEMS = 400          # start with a populated table so DELETE scans real data

function make_db(path)
    isfile(path) && rm(path)
    db = DBInterface.connect(DuckDB.DB, path)
    setup = DBInterface.connect(db)
    DBInterface.execute(setup, "CREATE TABLE shared(item_id TEXT, row_no BIGINT, val DOUBLE)")
    # Prefill so the table already holds data and the per-item DELETE scans a real table.
    DBInterface.execute(setup, "BEGIN TRANSACTION")
    app = DuckDB.Appender(db, "shared")
    for i in 1:PREFILL_ITEMS, r in 1:ROWS_PER_ITEM
        DuckDB.append(app, "prefill_$i"); DuckDB.append(app, Int64(r)); DuckDB.append(app, rand())
        DuckDB.end_row(app)
    end
    DuckDB.close(app)
    DBInterface.execute(setup, "COMMIT")
    DBInterface.close!(setup)
    return db
end

"""One item's write: replace its rows in a single transaction, exactly like the cache does."""
function write_item!(conn, item_id)
    DBInterface.execute(conn, "BEGIN TRANSACTION")
    try
        DBInterface.execute(
            DBInterface.prepare(conn, "DELETE FROM shared WHERE item_id = ?"), (item_id,))
        app = DuckDB.Appender(conn, "shared")
        for r in 1:ROWS_PER_ITEM
            DuckDB.append(app, item_id); DuckDB.append(app, Int64(r)); DuckDB.append(app, rand())
            DuckDB.end_row(app)
        end
        DuckDB.close(app)
        DBInterface.execute(conn, "COMMIT")
    catch err
        DBInterface.execute(conn, "ROLLBACK")
        rethrow()
    end
end

function run_serial(db)
    conn = DBInterface.connect(db)
    items = ["item_$i" for i in 1:N_ITEMS]
    elapsed = @elapsed for id in items
        write_item!(conn, id)
    end
    DBInterface.close!(conn)
    return elapsed
end

function run_concurrent(db, n_conns)
    items = ["item_$i" for i in 1:N_ITEMS]
    work = Channel{String}(length(items))
    for id in items; put!(work, id); end
    close(work)
    conflicts = Threads.Atomic{Int}(0)
    elapsed = @elapsed begin
        @sync for _ in 1:n_conns
            Threads.@spawn begin
                conn = DBInterface.connect(db)
                for id in work
                    while true
                        try
                            write_item!(conn, id); break
                        catch err
                            occursin("onflict", sprint(showerror, err)) || rethrow()
                            Threads.atomic_add!(conflicts, 1)
                        end
                    end
                end
                DBInterface.close!(conn)
            end
        end
    end
    return elapsed, conflicts[]
end

println("threads=$(Threads.nthreads())  items=$N_ITEMS  rows/item=$ROWS_PER_ITEM  prefill=$PREFILL_ITEMS")

db = make_db(tempname() * ".duckdb")
serial = run_serial(db)
@printf("serial (1 connection):        %7.0f ms  (%.0f items/s)\n", 1e3serial, N_ITEMS / serial)

for n in (2, 4, 8)
    n > Threads.nthreads() && continue
    db_n = make_db(tempname() * ".duckdb")
    el, conflicts = run_concurrent(db_n, n)
    @printf("concurrent (%d connections):   %7.0f ms  (%.0f items/s)  speedup %.2fx  conflicts=%d\n",
        n, 1e3el, N_ITEMS / el, serial / el, conflicts)
end
