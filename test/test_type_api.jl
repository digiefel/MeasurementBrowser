using DataBrowser
using Test
using GLMakie: Figure, Axis

const MB = DataBrowser

struct PhotoSource <: DataBrowserAPI.AbstractDataSource
    name::String
end

struct PhotoSourceItem <: DataBrowserAPI.AbstractDataSourceItem
    key::String
    exposure::Float64
    camera::String
    gain::Int
end

DataBrowserAPI.source_id(source::PhotoSource) = source.name
DataBrowserAPI.source_label(source::PhotoSource) = source.name
DataBrowserAPI.source_items(::PhotoSource; kwargs...) = [
    PhotoSourceItem("a.photo", 2.0, "typed camera", 1),
    PhotoSourceItem("b.photo", 4.0, "typed camera", 1),
]
DataBrowserAPI.source_item_id(item::PhotoSourceItem) = item.key
DataBrowserAPI.source_item_label(item::PhotoSourceItem) = item.key
DataBrowserAPI.metadata(item::PhotoSourceItem) =
    Dict(:camera => item.camera, :gain => item.gain)

struct Photo <: MB.AbstractDataItem
    exposure::Float64
    pixels::Matrix{Float64}
    camera::String
    gain::Int
    collection::Vector{String}
end

MB.collection(photo::Photo) = photo.collection
MB.metadata(photo::Photo) = Dict(
    :exposure => photo.exposure,
    :camera => photo.camera,
    :gain => photo.gain,
)

function MB.process(photo::Photo)
    return Photo(
        photo.exposure,
        photo.pixels .* photo.gain,
        photo.camera,
        photo.gain,
        photo.collection,
    )
end

MB.analyze(photo::Photo) =
    Dict(:mean_intensity => sum(photo.pixels) / length(photo.pixels))

function DataBrowserAPI.data_items(
    ::DataBrowserAPI.Project,
    ::PhotoSource,
    source_item::PhotoSourceItem,
)
    return [Photo(
        source_item.exposure,
        fill(source_item.exposure, 2, 2),
        source_item.camera,
        source_item.gain,
        ["micrographs"],
    )]
end

@testset "type API preserves custom AbstractDataItem values" begin
    project = MB.define_project("Photos")
    drawn_pixels = Ref(0)
    MB.register_plot!(project, :Photo;
        label="Image",
        setup=(workspace, items) -> Figure(),
        draw=function (workspace, items, figure)
            Axis(figure[1, 1])
            @test all(item -> item isa Photo, items)
            drawn_pixels[] += sum(length(item.pixels) for item in items)
            nothing
        end,
    )

    workspace = MB.open_workspace(
        project,
        PhotoSource("typed photos");
        cache=true,
        background_processing=true,
    )
    plot_kind = RegisteredPlot{:Photo,Symbol("Image")}
    try
        DataBrowserCore.Workspace.wait_workspace_idle!(workspace)
        records = DataBrowserAPI.ItemIndex.all_items(workspace.index.hierarchy)
        @test length(records) == 2
        @test all(record -> record.kind == :Photo, records)
        @test all(record -> !isempty(record.id), records)
        @test all(record -> record.item_label in ("a.photo", "b.photo"), records)
        @test all(record -> record.collection == ["micrographs"], records)
        @test Set(record.metadata[:exposure] for record in records) == Set([2.0, 4.0])
        @test all(record -> record.metadata[:camera] == "typed camera", records)
        @test all(record -> record.metadata[:gain] == 1, records)

        loaded = DataBrowserCore.Workspace.materialize_items(workspace, records)
        @test all(item -> item isa Photo, loaded)
        @test Set(item.exposure for item in loaded) == Set([2.0, 4.0])
        @test all(item -> item.camera == "typed camera", loaded)
        @test all(item -> item.gain == 1, loaded)

        figure = MB.setup_plot(workspace, plot_kind, records)
        @test figure isa Figure
        @test MB.plot_data!(workspace, plot_kind, records, figure) === nothing
        @test drawn_pixels[] == 8
    finally
        MB.close_workspace!(workspace)
    end
end
