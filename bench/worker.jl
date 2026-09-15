using TOML, Statistics

workload, output = ARGS[1:2]
import_s = @elapsed @eval using DataBrowser

if workload == "engine"
    include("workloads.jl")
    repeats = parse(Int, get(ENV, "MB_BENCH_REPEATS", "3"))
    repeats > 0 || error("MB_BENCH_REPEATS must be positive")
    samples = run_benchmark(repeats)
    result = Dict(String(name) => median(getproperty.(samples, name)) for name in keys(first(samples)))
else
    include("browser.jl")
    root, depot, opening, started_ns = ARGS[3:6]
    # Change only application-cache storage after packages have loaded.
    pushfirst!(DEPOT_PATH, depot)
    result = try
        browser_workload(root, opening, parse(UInt64, started_ns))
    finally
        popfirst!(DEPOT_PATH)
    end
    result["import_s"] = import_s
end

open(io -> TOML.print(io, result; sorted=true), output, "w")
