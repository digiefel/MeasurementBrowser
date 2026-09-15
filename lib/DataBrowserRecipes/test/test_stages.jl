using DataBrowserAPI: AbstractDataSource, AbstractDataSourceItem, entries, process, analyze, reconstruct, id, label, item_data, metadata
using Test, DataBrowserRecipes, DataBrowserAPI

struct RecipeSource <: AbstractDataSource end
struct RecipeFile <: AbstractDataSourceItem
    path::String
end
DataBrowserAPI.id(f::RecipeFile) = basename(f.path)
DataBrowserAPI.label(f::RecipeFile) = basename(f.path)
DataBrowserAPI.source_item_path(f::RecipeFile) = f.path

@testset "registered callbacks implement the type API" begin
    project = define_project("recipe")
    register_item!(project, :trace;
        read=_ -> (data=[1, 2], metadata=Dict(:scale => 3)),
        entries=(data, _) -> [(data=(x=data,), metadata=Dict(:channel => n)) for n in 1:2],
        id=(_, meta) -> meta[:channel],
        process=(data, meta) -> (x=data.x .* meta[:scale],),
        analyze=(data, _) -> Dict(:total => sum(data.x)))
    file = RecipeFile("trace.dat")
    items = entries(project, file, read(project, RecipeSource(), file))
    @test length(unique(id.(items))) == 2
    for item in items
        processed = process(project, item)
        @test id(processed) == id(item)
        @test item_data(processed).x == [3, 6]
        @test analyze(project, processed)[:total] == 9
        restored = reconstruct(typeof(processed), id(processed), item_data(processed), metadata(processed))
        @test metadata(restored) == metadata(processed)
    end
end

@testset "CSV recipe uses the source-item contract" begin
    mktempdir() do dir
        file = RecipeFile(joinpath(dir, "sample.csv"))
        write(file.path, "x,y\n1,2\n3,4\n")
        project = register_csv!(define_project("csv"))
        item = only(entries(project, file, read(project, RecipeSource(), file)))
        @test item_data(item).x == [1, 3]
        @test item_data(item).y == [2, 4]
    end
end
