using TOML, Statistics

const RESULTS_DIR = joinpath(@__DIR__, "results")

"""Run a worker and sample its resident memory, including native allocations, every 100 ms."""
function run_worker(arguments)
    command = `$(Base.julia_cmd()) --startup-file=no --project=$(@__DIR__) --threads=$(Threads.nthreads()) $(@__FILE__) --worker $arguments`
    process = run(pipeline(command; stdout=stdout, stderr=stderr); wait=false)
    peak_kib = 0
    try
        while process_running(process)
            rss = try
                read(pipeline(`ps -o rss= -p $(getpid(process))`, stderr=devnull), String)
            catch
                process_running(process) && rethrow()
                ""
            end
            isempty(strip(rss)) || (peak_kib = max(peak_kib, parse(Int, strip(rss))))
            sleep(0.1)
        end
        success(process) || error("Benchmark worker failed")
    finally
        if process_running(process)
            kill(process)
            wait(process)
        end
    end
    return peak_kib / 1024
end

function benchmark_main(workload)
    mkpath(RESULTS_DIR)
    output = joinpath(RESULTS_DIR, workload * ".toml")
    rm(output; force=true)
    if workload == "engine"
        rss = run_worker(["engine", output])
        result = TOML.parsefile(output)
        result["peak_rss_mib"] = rss
    elseif workload == "browser"
        result = mktempdir() do dir
            root = joinpath(dir, "project")
            cp(joinpath(@__DIR__, "..", "test", "fixtures", "public_api"), root)
            depot = joinpath(dir, "data-cache")
            mkpath(depot)
            peaks = Float64[]
            measurements = Dict{String,Any}()
            for opening in ("open", "reopen")
                part = joinpath(dir, opening * ".toml")
                push!(peaks, run_worker(["browser", part, root, depot, opening, string(time_ns())]))
                merge!(measurements, TOML.parsefile(part))
            end
            measurements["peak_rss_mib"] = maximum(peaks)
            measurements
        end
    else
        error("Unknown benchmark: $workload")
    end
    open(io -> TOML.print(io, result; sorted=true), output, "w")
end

function benchmark_worker(args)
    workload, output = args[1:2]
    import_s = @elapsed @eval using DataBrowser
    if workload == "engine"
        include("workloads.jl")
        result = Base.invokelatest() do
            repeats = parse(Int, get(ENV, "MB_BENCH_REPEATS", "3"))
            repeats > 0 || error("MB_BENCH_REPEATS must be positive")
            samples = run_benchmark(repeats)
            Dict(String(name) => median(getproperty.(samples, name)) for name in keys(first(samples)))
        end
    else
        include("browser.jl")
        root, depot, opening, started_ns = args[3:6]
        # Packages are loaded before the temporary application-cache depot is installed.
        pushfirst!(DEPOT_PATH, depot)
        result = try
            Base.invokelatest() do
                browser_workload(root, opening, parse(UInt64, started_ns))
            end
        finally
            popfirst!(DEPOT_PATH)
        end
        result["import_s"] = import_s
    end
    open(io -> TOML.print(io, result; sorted=true), output, "w")
end

if abspath(PROGRAM_FILE) == @__FILE__
    if !isempty(ARGS) && first(ARGS) == "--worker"
        benchmark_worker(ARGS[2:end])
    else
        foreach(benchmark_main, isempty(ARGS) ? ["engine", "browser"] : ARGS)
    end
end
