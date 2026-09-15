using Printf

"""Write the latest measured metrics with their input identities and measurement dates."""
function write_report(records)
    path = joinpath(@__DIR__, "status.txt")
    open(path, "w") do io
        println(io, "# Julia $(VERSION); $(Threads.nthreads()) threads; $(Sys.CPU_NAME); $(Sys.KERNEL)/$(Sys.ARCH)")
        for name in ("precompile", "engine", "browser")
            haskey(records, name) || continue
            record = records[name]
            println(io, "# $name: $(record["measured_at"]); inputs $(record["inputs"])")
            values = name == "precompile" ? Dict("precompile_s" => record["seconds"]) : record["result"]
            for key in sort(collect(keys(values)))
                @printf(io, "%-30s %12.4f\n", name == "browser" && key == "peak_rss_mib" ? "browser_peak_rss_mib" : key, values[key])
            end
        end
    end
    print(read(path, String))
end
