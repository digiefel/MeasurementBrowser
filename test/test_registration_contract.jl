const REGISTRATION_API = DataBrowserAPI
const REGISTRATION_CORE = DataBrowserCore
const REGISTRATION_SOURCES = DataBrowserSources

@testset "registration callbacks receive data and metadata" begin
    directory = mktempdir()
    filepath = joinpath(directory, "cycles.dat")
    write(filepath, "2,4")
    source = REGISTRATION_SOURCES.DirectorySource(directory; metadata_file=nothing)
    source_file = REGISTRATION_SOURCES.index_source_file(filepath)

    seen_entries_metadata = Ref{Dict}()
    project = REGISTRATION_API.define_project("Registration contract")
    REGISTRATION_API.register_item!(project, :cycles;
        detect=(file::REGISTRATION_SOURCES.SourceFile) -> endswith(file.filename, ".dat"),
        read=(file::REGISTRATION_SOURCES.SourceFile) -> (
            data=parse.(Int, split(read(file.filepath, String), ',')),
            metadata=Dict(:instrument => "tester"),
        ),
        entries=(values::Vector{Int}, source_metadata::Dict) -> begin
            seen_entries_metadata[] = copy(source_metadata)
            [
                (data=value, metadata=Dict(:position => position))
                for (position, value) in pairs(values)
            ]
        end,
        label=(value::Int, item_metadata::Dict) ->
            "$(item_metadata[:filename]) #$(item_metadata[:position])",
        collection=(value::Int, item_metadata::Dict) ->
            ["instrument $(item_metadata[:instrument])"],
        process=(value::Int, item_metadata::Dict) -> value * item_metadata[:scale],
        analyze=(value::Int, item_metadata::Dict) -> Dict(:squared => value^2),
    )

    items = REGISTRATION_API.data_items(project, source, source_file)
    @test length(items) == 2
    @test seen_entries_metadata[][:filename] == "cycles.dat"
    @test seen_entries_metadata[][:instrument] == "tester"
    @test REGISTRATION_API.item_data.(items) == [2, 4]
    @test REGISTRATION_API.metadata(items[1]) == Dict(
        :filename => "cycles.dat",
        :instrument => "tester",
        :position => 1,
    )
    @test REGISTRATION_API.label.(items) == ["cycles.dat #1", "cycles.dat #2"]
    @test [REGISTRATION_API.label.(REGISTRATION_API.collection(item)) for item in items] ==
        [["instrument tester"], ["instrument tester"]]
    # `data_items` adapts entries without minting identity; ids are minted at interpretation.
    @test REGISTRATION_API.id(items[1]) == REGISTRATION_API.id(items[2]) == ""
    interpretation = REGISTRATION_CORE.interpret_source_item(project, source, source_file)
    @test length(interpretation.records) == 2
    @test allunique(record.id for record in interpretation.records)
    @test [REGISTRATION_API.id(item) for item in interpretation.interpreted_items] ==
        [record.id for record in interpretation.records]

    effective = merge(REGISTRATION_API.metadata(items[1]), Dict(:scale => 10))
    input = REGISTRATION_API.ItemIndex.RegisteredDataItem(
        REGISTRATION_API.ItemIndex.ItemRecord(
            interpretation.records[1]; metadata=effective),
        REGISTRATION_API.item_data(items[1]),
    )
    processed = REGISTRATION_API.process(project, source, input)
    @test REGISTRATION_API.item_data(processed) == 20
    @test REGISTRATION_API._analyze_item(project, source, processed) ==
        Dict(:squared => 400)
end

@testset "registration collection callback requires a string vector" begin
    directory = mktempdir()
    filepath = joinpath(directory, "item.dat")
    write(filepath, "data")
    source = REGISTRATION_SOURCES.DirectorySource(directory; metadata_file=nothing)
    source_file = REGISTRATION_SOURCES.index_source_file(filepath)
    project = REGISTRATION_API.define_project("Invalid collection")
    REGISTRATION_API.register_item!(
        project;
        read=_ -> :data,
        collection=(_data, _metadata) -> "flat",
    )
    @test_throws ArgumentError REGISTRATION_API.data_items(project, source, source_file)
end

@testset "unnamed registration is optional and replaceable" begin
    project = REGISTRATION_API.define_project("Unnamed")
    REGISTRATION_API.register_item!(project; read=(file) -> :first)
    @test length(project.recipes) == 1
    @test project.recipes[1].read(nothing) == :first

    REGISTRATION_API.register_item!(project; read=(file) -> :second)
    @test length(project.recipes) == 1
    @test project.recipes[1].read(nothing) == :second
    @test project.recipes[1].detect(nothing)
    @test project.recipes[1].entries === nothing
end
