using TOML

"""Run a worker and sample its resident memory, including native allocations, every 100 ms."""
function run_worker(command::Cmd)
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

"""Measure one workload in fresh Julia processes using the supplied Julia command."""
function measure_benchmark(workload, julia::Cmd)
    worker = joinpath(@__DIR__, "worker.jl")
    return mktempdir() do dir
        output = joinpath(dir, "result.toml")
        if workload == "engine"
            rss = run_worker(`$julia $worker engine $output`)
            result = TOML.parsefile(output)
        elseif workload == "browser"
            root = joinpath(dir, "project")
            cp(joinpath(@__DIR__, "..", "test", "fixtures", "public_api"), root)
            depot = mkpath(joinpath(dir, "data-cache"))
            rss = 0.0
            result = Dict{String,Any}()
            for opening in ("open", "reopen")
                peak = run_worker(`$julia $worker browser $output $root $depot $opening $(time_ns())`)
                rss = max(rss, peak)
                measurements = TOML.parsefile(output)
                # Import and frame time refer to the initial opening; reopen has its own latency.
                opening == "reopen" && delete!(measurements, "import_s")
                opening == "reopen" && delete!(measurements, "frame_ui_ms")
                merge!(result, measurements)
            end
        else
            error("Unknown benchmark: $workload")
        end
        result["peak_rss_mib"] = rss
        return result
    end
end
