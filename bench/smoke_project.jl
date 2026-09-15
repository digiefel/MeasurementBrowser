using DataBrowser
using GLMakie: Figure, Axis, lines!
using TOML

const PUBLIC_FIXTURE = joinpath(@__DIR__, "..", "test", "fixtures", "public_api")
const PUBLIC_VARIANTS = joinpath(@__DIR__, "..", "test", "fixtures", "public_api_variants")

struct ToyCounters
    reads::Dict{String,Threads.Atomic{Int}}
    entries::Threads.Atomic{Int}
    processes::Threads.Atomic{Int}
    analyses::Threads.Atomic{Int}
    collection_processes::Threads.Atomic{Int}
    collection_analyses::Threads.Atomic{Int}
    draws::Threads.Atomic{Int}
end

ToyCounters() = ToyCounters(
    Dict(name => Threads.Atomic{Int}(0) for name in
        ("a.dbitem", "b.dbitem", "extra.dbitem")),
    Threads.Atomic{Int}(0), Threads.Atomic{Int}(0), Threads.Atomic{Int}(0),
    Threads.Atomic{Int}(0), Threads.Atomic{Int}(0), Threads.Atomic{Int}(0),
)

function toy_project(name::AbstractString, counters::ToyCounters)::Project
    project = define_project(name)
    register_item!(project, :trace;
        detect=file -> endswith(file.filename, ".dbitem"),
        read=function (file)
            Threads.atomic_add!(counters.reads[file.filename], 1)
            return (data=TOML.parsefile(file.filepath), metadata=Dict(:instrument => "toy"))
        end,
        entries=function (source, _metadata)
            Threads.atomic_add!(counters.entries, 1)
            values = Float64.(source["values"])
            return [(
                data=(x=Float64.(eachindex(values)), y=values .+ (index - 1)),
                metadata=Dict{Symbol,Any}(
                    :device => source["device"], :channel => channel),
            ) for (index, channel) in pairs(source["channels"])]
        end,
        id=(_data, metadata) -> metadata[:channel],
        label=(_data, metadata) -> "$(metadata[:filename])/$(metadata[:channel])",
        collection=(_data, metadata) -> ["runs", metadata[:device]],
        process=function (data, metadata)
            Threads.atomic_add!(counters.processes, 1)
            return (x=data.x, y=data.y .* metadata[:scale])
        end,
        analyze=function (data, _metadata)
            Threads.atomic_add!(counters.analyses, 1)
            return Dict{Symbol,Any}(:peak => maximum(data.y))
        end,
    )
    register_collection_analysis!(project, :trace;
        process=function (items, _metadata)
            Threads.atomic_add!(counters.collection_processes, 1)
            return [(; item..., members=fill(length(items), length(item.x))) for item in items]
        end,
        analyze=function (items, _metadata)
            Threads.atomic_add!(counters.collection_analyses, 1)
            return Dict{Symbol,Any}(:members => length(items))
        end,
    )
    register_plot!(project, :trace;
        label="Toy",
        setup=(_workspace, _items) -> begin
            figure = Figure()
            Axis(figure[1, 1])
            figure
        end,
        draw=function (_workspace, items, _figure)
            all(item -> label(typeof(item)) === :trace && haskey(metadata(item), :peak), items) ||
                error("plot did not receive analyzed trace items")
            axis = _figure[1, 1][]
            empty!(axis)
            for item in items
                lines!(axis, item_data(item).x, item_data(item).y)
            end
            Threads.atomic_add!(counters.draws, 1)
            return nothing
        end,
    )
    return project
end

function copied_public_fixture(work::Function)
    mktempdir() do dir
        root = joinpath(dir, "project")
        cp(PUBLIC_FIXTURE, root)
        work(root)
    end
end
