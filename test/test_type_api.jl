using MeasurementBrowser
using Test
using GLMakie: Figure, Axis

const MB = MeasurementBrowser

# A custom item via the type API: subtype AbstractDataItem, carry metadata as typed fields and the
# payload directly. The engine indexes it through the same contract the package's DataItem answers.
struct Photo <: MB.AbstractDataItem
    id::String
    collection::Vector{String}
    exposure::Float64
    data::Matrix{Float64}
end
MB.item_id(p::Photo) = p.id
MB.item_label(p::Photo) = "exp=$(p.exposure)"
MB.kind(::Photo) = :photo
MB.collection(p::Photo) = p.collection
MB.parameters(p::Photo) = Dict{Symbol,Any}(:exposure => p.exposure)
MB.stats(p::Photo) = Dict{Symbol,Any}()
MB.item_data(p::Photo) = p.data

@testset "type API: custom AbstractDataItem end to end" begin
    dir = mktempdir()
    write(joinpath(dir, "a.photo"), "2.0\n")
    write(joinpath(dir, "b.photo"), "4.0\n")

    project = MB.define_project("Photos")
    MB.register_item!(
        project,
        :photo;
        detect=file -> endswith(file.filename, ".photo"),
        read=file -> parse(Float64, strip(read(file.filepath, String))),
        entries=function (file, exposure)
            name = splitext(file.filename)[1]
            return [Photo(file.filepath, [name], exposure, fill(exposure, 2, 2))]
        end,
    )

    drawn_pixels = Ref(0)
    MB.register_plot!(
        project,
        :photo;
        label="Image",
        setup=(ws, items) -> Figure(),
        draw=function (ws, items, figure)
            Axis(figure[1, 1])
            for it in items
                # The bridge hands back the project's own subtype, carrying its data.
                @test it isa Photo
                @test MB.item_data(it) === it.data
                drawn_pixels[] += length(it.data)
            end
            nothing
        end,
    )

    # Scan: entries returns Photos; the engine derives records via the contract (no filepath needed
    # on the item — it comes from the SourceFile) and frees the data-bearing items.
    scan = scan_test_source(project, dir)
    records = scan.hierarchy.all_items
    @test length(records) == 2
    @test all(r -> r.kind == :photo, records)
    @test Set(r.parameters[:exposure] for r in records) == Set([2.0, 4.0])
    @test all(r -> !isempty(r.collection), records)

    # Plot: the bridge re-runs read+entries for the type API and matches items to records by id,
    # so draw receives the Photos with item.data populated.
    workspace = MB.Workspace.Workspace(project, test_source(project, dir))
    plot_kind = MB.RegisteredPlot{:photo,Symbol("Image")}
    figure = MB.setup_plot(workspace, plot_kind, records)
    @test figure isa Figure
    @test MB.plot_data!(workspace, plot_kind, records, figure) === nothing
    @test drawn_pixels[] == 8   # two 2x2 matrices
end
