using DataBrowser
using GLMakie: Figure
using Test
using TOML

const PUBLIC_FIXTURE = joinpath(@__DIR__, "fixtures", "public_api")
const PUBLIC_VARIANTS = joinpath(@__DIR__, "fixtures", "public_api_variants")

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
        setup=(_workspace, _items) -> Figure(),
        draw=function (_workspace, items, _figure)
            all(item -> kind(item) === :trace && haskey(metadata(item), :peak), items) ||
                error("plot did not receive analyzed trace items")
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

@testset "public registration rejects an invalid collection path" begin
    project = define_project("InvalidCollection")
    register_item!(project;
        read=file -> TOML.parsefile(file.filepath),
        collection=(_data, _metadata) -> "flat",
    )
    @test_throws ArgumentError items_for_file(project, joinpath(PUBLIC_FIXTURE, "a.dbitem"))
end

@testset "public API drives the complete registered pipeline" begin
    copied_public_fixture() do root
        counters = ToyCounters()
        project = toy_project("PublicAPI_$(basename(dirname(root)))", counters)
        workspace = open_workspace(project, root; background_processing=true)
        ids = String[]
        reads_after_build = Dict{String,Int}()
        try
            wait_workspace_idle!(workspace; timeout=30)
            status = workspace_status(workspace)
            @test !status.busy && isempty(status.errors)
            @test status.counts.cache.interpreted_items == 3
            @test status.counts.cache.processed == 3
            @test status.counts.cache.analyzed == 3
            @test status.counts.cache.collection_processed == 2
            @test status.counts.cache.collection_analyzed == 2

            ids = query_items(workspace)
            @test length(ids) == 3
            select_items!(workspace, ids)
            items = materialize_items(workspace)
            @test read_item_data(workspace) == item_data.(items)
            @test Set(item_label.(items)) ==
                Set(["a.dbitem/up", "a.dbitem/down", "b.dbitem/only"])
            @test all(item -> collection(item)[1] == "runs", items)
            @test sort([metadata(item)[:peak] for item in items]) == [4.0, 6.0, 9.0]
            @test sort([first(item_data(item).members) for item in items]) == [1, 2, 2]

            plot_kind = only(registered_plot_kinds(project, :trace))
            figure = setup_plot(workspace, plot_kind, items)
            @test plot_data!(workspace, plot_kind, items, figure) === nothing
            @test counters.draws[] == 1

            profile = only(row for row in scan_profile_summary(project) if row.kind === :trace)
            @test profile.source_items == 2
            @test profile.items == 3
            @test profile.detect_seconds >= 0
            @test profile.read_seconds >= 0
            @test profile.entries_seconds >= 0
            @test profile.process_seconds >= 0
            @test profile.analyze_seconds >= 0
            @test profile.total_seconds >= profile.read_seconds

            source_profile = filter(row -> row.kind === :trace, scan_source_profile(project))
            @test length(source_profile) == 2
            @test all(row -> row.items > 0 && !isempty(row.thread_ids), source_profile)
            @test Set(row.source_item_label for row in source_profile) ==
                Set(["a.dbitem", joinpath("nested", "b.dbitem")])
            @test Set(row.source_item_path for row in source_profile) == Set([
                joinpath(root, "a.dbitem"),
                joinpath(root, "nested", "b.dbitem"),
            ])
            @test all(row -> !isabspath(row.source_item_label), source_profile)
            @test issorted(source_profile; by=row -> row.total_seconds, rev=true)
            @test counters.collection_processes[] == 2
            @test counters.collection_analyses[] == 2
            reads_after_build = Dict(name => count[] for (name, count) in counters.reads)
        finally
            close_workspace!(workspace)
        end

        reopened = open_workspace(project, root; background_processing=true)
        try
            wait_workspace_idle!(reopened; timeout=30)
            @test query_items(reopened) == ids
            @test query_items(reopened, "peak >= 9") == [
                only(id for id in ids if occursin("b.dbitem", id))]
            select_items!(reopened, ids)
            @test length(materialize_items(reopened)) == 3
            @test Dict(name => count[] for (name, count) in counters.reads) == reads_after_build
        finally
            close_workspace!(reopened)
        end

        cp(joinpath(PUBLIC_VARIANTS, "a_changed.dbitem"), joinpath(root, "a.dbitem"); force=true)
        changed = open_workspace(project, root; background_processing=true)
        try
            wait_workspace_idle!(changed; timeout=30)
            select_items!(changed, query_items(changed))
            @test sort([metadata(item)[:peak] for item in materialize_items(changed)]) ==
                [9.0, 40.0, 42.0]
            @test counters.reads["a.dbitem"][] == reads_after_build["a.dbitem"] + 1
            @test counters.reads["b.dbitem"][] == reads_after_build["b.dbitem"]

            reads_before_metadata = Dict(name => count[] for (name, count) in counters.reads)
            processes_before_metadata = counters.processes[]
            cp(joinpath(PUBLIC_VARIANTS, "metadata_changed.txt"),
                joinpath(root, "metadata.txt"); force=true)
            @test Base.timedwait(
                () -> counters.processes[] >= processes_before_metadata + 2, 5) === :ok
            wait_workspace_idle!(changed; timeout=30)
            select_items!(changed, query_items(changed))
            @test sort([metadata(item)[:peak] for item in materialize_items(changed)]) ==
                [9.0, 80.0, 84.0]
            @test Dict(name => count[] for (name, count) in counters.reads) ==
                reads_before_metadata
        finally
            close_workspace!(changed)
        end
    end
end

@testset "public API observes recursive and live source changes" begin
    copied_public_fixture() do root
        counters = ToyCounters()
        project = toy_project("PublicLive_$(basename(dirname(root)))", counters)

        shallow = open_workspace(project, root; recursive=false, cache=false)
        try
            wait_workspace_idle!(shallow; timeout=30)
            @test length(query_items(shallow)) == 2
        finally
            close_workspace!(shallow)
        end

        workspace = open_workspace(project, root; cache=false)
        try
            wait_workspace_idle!(workspace; timeout=30)
            @test length(query_items(workspace)) == 3

            extra = joinpath(root, "extra.dbitem")
            cp(joinpath(PUBLIC_VARIANTS, "extra.dbitem"), extra)
            @test Base.timedwait(() -> length(query_items(workspace)) == 4, 5) === :ok

            rm(extra)
            @test Base.timedwait(() -> length(query_items(workspace)) == 3, 5) === :ok

            cp(joinpath(PUBLIC_VARIANTS, "metadata_malformed.txt"),
                joinpath(root, "metadata.txt"); force=true)
            @test Base.timedwait(
                () -> workspace_status(workspace).label == "Source Error", 5) === :ok

            cp(joinpath(PUBLIC_FIXTURE, "metadata.txt"),
                joinpath(root, "metadata.txt"); force=true)
            @test Base.timedwait(
                () -> workspace_status(workspace).label != "Source Error", 5) === :ok
        finally
            close_workspace!(workspace)
        end
    end
end
