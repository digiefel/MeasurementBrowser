# Package isolation and resolution belong to Pkg. This runner only chooses workloads.
using Pkg, SHA, TOML, Dates

const ROOT = dirname(@__DIR__)
const BENCH = joinpath(ROOT, "bench")
const RESULTS = joinpath(BENCH, "results")
const PACKAGES = sort(["DataBrowser" * name for name in
    ("API", "Annotations", "Cache", "Core", "GUI", "Plots", "Profiling", "Recipes", "Sources")])
const RECORD_FILE = joinpath(RESULTS, "checks.toml")

function hash_files(io, path)
    ispath(path) || return
    files = isfile(path) ? [path] : sort([joinpath(dir, file)
        for (dir, _, names) in walkdir(path) for file in names])
    for file in files
        println(io, relpath(file, ROOT))
        write(io, read(file))
    end
end

"""Fingerprint resolved dependencies and local source contents, including uncommitted edits."""
function fingerprint(names; tests=false, workload=nothing)
    dependencies = Pkg.dependencies()
    byname = Dict(info.name => uuid for (uuid, info) in dependencies)
    visited = Set()
    io = IOBuffer()
    println(io, VERSION, Sys.MACHINE, Sys.CPU_NAME, Threads.nthreads())
    hash_files(io, joinpath(BENCH, "LocalPreferences.toml"))
    hash_files(io, joinpath(ROOT, "LocalPreferences.toml"))
    function visit(uuid)
        uuid in visited && return
        push!(visited, uuid)
        info = dependencies[uuid]
        println(io, uuid, info.version, info.tree_hash)
        if info.is_tracking_path
            for name in ("Project.toml", "src", "ext", "deps", "LocalPreferences.toml")
                hash_files(io, joinpath(info.source, name))
            end
        end
        for dependency in sort(collect(values(info.dependencies)); by=string)
            visit(dependency)
        end
    end
    for name in sort(names)
        visit(byname[name])
        tests && hash_files(io, joinpath(dependencies[byname[name]].source, "test"))
    end
    (tests || workload !== nothing) && hash_files(io, @__FILE__)
    if workload !== nothing
        for file in workload
            hash_files(io, joinpath(ROOT, file))
        end
        println(io, get(ENV, "MB_BENCH_REPEATS", "3"))
    end
    return bytes2hex(sha256(take!(io)))
end

"""Run one command only when its inputs differ from the last successful run."""
function checked(f, records, name, key; force=false)
    if !force && get(get(records, name, Dict()), "inputs", nothing) == key
        println(name, ": unchanged (previous success reused)")
        return records[name]
    end
    println(name, ": running")
    result = try
        f()
    catch
        delete!(records, name)
        mkpath(RESULTS)
        open(io -> TOML.print(io, records; sorted=true), RECORD_FILE, "w")
        rethrow()
    end
    record = Dict{String,Any}("inputs" => key, "measured_at" => string(Dates.now()))
    result === nothing || merge!(record, result)
    records[name] = record
    mkpath(RESULTS)
    open(io -> TOML.print(io, records; sorted=true), RECORD_FILE, "w")
    return record
end

"""An isolated compiled cache, sharing installed package sources and binary artifacts."""
function measurement_depot()
    depot = joinpath(RESULTS, "depot")
    mkpath(depot)
    for name in ("packages", "artifacts", "registries")
        target = joinpath(depot, name)
        if !ispath(target)
            source = findfirst(d -> isdir(joinpath(d, name)), DEPOT_PATH)
            source === nothing || symlink(joinpath(DEPOT_PATH[source], name), target)
        end
    end
    return depot
end

function julia_command(arguments; depot=nothing)
    command = `$(Base.julia_cmd()) --startup-file=no --project=$BENCH --threads=$(Threads.nthreads()) $arguments`
    return depot === nothing ? command : addenv(command, "JULIA_DEPOT_PATH" => depot * string(Sys.iswindows() ? ';' : ':'))
end

"""Run selected package tests; the full selection also runs clean compilation and smoke benchmarks."""
function check(args=ARGS)
    force = "--force" in args
    selected = filter(!=("--force"), args)
    full = isempty(selected)
    package = full ? nothing : startswith(first(selected), "DataBrowser") ? first(selected) : "DataBrowser" * first(selected)
    benchmark = !full && first(selected) == "bench"
    precompile_only = !full && first(selected) == "precompile"
    full || benchmark || precompile_only || package in PACKAGES || error("Unknown package: $(first(selected))")
    records = isfile(RECORD_FILE) ? TOML.parsefile(RECORD_FILE) : Dict{String,Any}()
    depot = measurement_depot()
    code = fingerprint(vcat(PACKAGES, ["DataBrowser"]))
    if full || benchmark || precompile_only
        checked(records, "precompile", code; force=force && precompile_only) do
            # Only this runner-owned compiled directory is removed. Sources/artifacts are retained.
            rm(joinpath(depot, "compiled"); recursive=true, force=true)
            command = julia_command(["-e", "using Pkg; Pkg.precompile($(repr(vcat(PACKAGES, ["DataBrowser"]))); strict=true)"]; depot)
            seconds = @elapsed run(command)
            fingerprint(vcat(PACKAGES, ["DataBrowser"])) == code || error("Package code changed during precompilation")
            Dict("seconds" => seconds)
        end
    end
    if full || (!benchmark && !precompile_only)
        for name in (full ? PACKAGES : [package])
            files = full ? String[] : selected[2:end]
            key = fingerprint([name]; tests=true)
            # A selected file is an explicit diagnostic run, not a pass for the whole package.
            test_code = "using Pkg; Pkg.test($(repr(name)); allow_reresolve=false, julia_args=[\"--check-bounds=auto\"], test_args=$(repr(files)))"
            command = julia_command(["-e", test_code]; depot=isdir(joinpath(depot, "compiled")) ? depot : nothing)
            if isempty(files)
                checked(records, name, key; force) do
                    run(command)
                    fingerprint([name]; tests=true) == key || error("$name changed while its tests ran")
                    nothing
                end
            else
                run(command)
            end
        end
    end
    if full || benchmark
        workloads = full || length(selected) == 1 ? ["engine", "browser"] : selected[2:end]
        for workload in workloads
            workload in ("engine", "browser") || error("Unknown benchmark: $workload")
            files = workload == "engine" ? ["bench/run.jl", "bench/workloads.jl"] :
                ["bench/run.jl", "bench/browser.jl", "bench/smoke_project.jl", "test/fixtures/public_api", "test/fixtures/public_api_variants"]
            key = fingerprint(["DataBrowser"]; workload=files)
            checked(records, workload, key; force) do
                run(julia_command([joinpath(BENCH, "run.jl"), workload]; depot))
                fingerprint(["DataBrowser"]; workload=files) == key || error("Benchmark inputs changed during the run")
                Dict("result" => TOML.parsefile(joinpath(RESULTS, workload * ".toml")))
            end
        end
        if workloads == ["engine", "browser"]
            include(joinpath(BENCH, "report.jl"))
            Base.invokelatest() do
                write_report(records)
            end
        end
    end
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    check()
end
