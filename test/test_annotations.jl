using Annotations
using Test

const ANNOT_FIXTURES = joinpath(@__DIR__, "fixtures", "annotations")

@testset "Annotations" begin

    @testset "Coords" begin
        info = Annotations.Coords.DevicesInfoTable(
            ("ChipA",) => Dict{Symbol,Any}(:notes => "root only"),
            ("ChipA", "Site1") => Dict{Symbol,Any}(:x_um => 0.0, :y_um => 0.0,
                                                    :w_um => 100.0, :h_um => 50.0),
            ("ChipA", "Site2") => Dict{Symbol,Any}(:x_um => "200.0", :y_um => "150.0"),
            ("ChipA", "Site3") => Dict{Symbol,Any}(:x_um => 50.0),
        )
        positions = Annotations.Coords.read_positions(info)
        @test positions["ChipA/Site1"] == (0.0, 0.0)
        @test positions["ChipA/Site2"] == (200.0, 150.0)
        @test !haskey(positions, "ChipA/Site3")
        @test !haskey(positions, "ChipA")

        overrides = Annotations.Coords.read_overrides(info)
        @test overrides["ChipA/Site1"] == (100.0, 50.0)
        @test !haskey(overrides, "ChipA/Site2")

        bbox = Annotations.Coords.bounding_box(positions, ["ChipA/Site1", "ChipA/Site2"])
        @test bbox === Annotations.Coords.Rect(0.0, 0.0, 200.0, 150.0)

        bbox_override = Annotations.Coords.bounding_box(
            positions, ["ChipA/Site1", "ChipA/Site2"]; override=(500.0, 400.0))
        @test bbox_override.x == 0.0 && bbox_override.y == 0.0
        @test bbox_override.w == 500.0 && bbox_override.h == 400.0

        @test Annotations.Coords.bounding_box(positions, ["Nonexistent/Path"]) === nothing
        @test Annotations.Coords.read_positions(nothing) == Dict{String,Tuple{Float64,Float64}}()
    end

    @testset "Layout" begin
        mktempdir() do dir
            positions = Annotations.Layout.PositionMap(
                "ChipA" => (10.0, 20.0),
                "ChipB" => (-5.5, 0.0),
            )
            Annotations.Layout.save(dir, positions)
            path = Annotations.Layout.layout_path(dir)
            @test isfile(path)
            roundtrip = Annotations.Layout.load(dir)
            @test roundtrip == positions

            empty_map = Annotations.Layout.PositionMap()
            Annotations.Layout.save(dir, empty_map)
            @test !isfile(path)
        end

        @testset "reset! grid shapes" begin
            for (n, expected_cols) in ((1, 1), (4, 2), (7, 3))
                positions = Annotations.Layout.PositionMap()
                paths = ["P$i" for i in 1:n]
                Annotations.Layout.reset!(positions, paths; spacing_um=10.0)
                @test length(positions) == n
                xs = sort(unique(p[1] for p in values(positions)))
                @test length(xs) == min(expected_cols, n)
                # row-major fill: P1 at origin
                @test positions["P1"] == (0.0, 0.0)
            end
        end

        mktempdir() do dir
            write(joinpath(dir, "layout.txt"), "bad row only one field\n")
            @test_throws Annotations.Layout.LayoutParseError Annotations.Layout.load(dir)
        end
    end

    @testset "Tags" begin
        @testset "missing files" begin
            mktempdir() do dir
                state = Annotations.Tags.load(dir)
                @test isempty(state.catalog)
                @test isempty(state.assignments)
                @test isempty(state.measurement_assignments)
            end
        end

        @testset "parse and roundtrip" begin
            fixture = joinpath(ANNOT_FIXTURES, "tags_full")
            state = Annotations.Tags.load(fixture)
            names = [t.name for t in state.catalog]
            @test "bad" in names && "todo" in names
            bad = state.catalog[findfirst(t -> t.name == "bad", state.catalog)]
            @test bad.color == (0xff, 0x30, 0x30)
            @test bad.priority == 100
            @test "bad" in state.assignments["RuO2test/A9/VI/D1"]
            @test "todo" in state.assignments["RuO2test/A10/VI"]
            @test isempty(state.measurement_assignments)

            mktempdir() do tmp
                Annotations.Tags.save(tmp, state)
                roundtripped = Annotations.Tags.load(tmp)
                @test Set(t.name for t in roundtripped.catalog) ==
                    Set(t.name for t in state.catalog)
                @test roundtripped.assignments == state.assignments
                @test roundtripped.measurement_assignments == state.measurement_assignments
            end
        end

        @testset "legacy bad_measurements migration — device only" begin
            mktempdir() do dir
                write(joinpath(dir, "bad_measurements"),
                    """
                    # comment
                    device RuO2test/A9/VI/D1
                    device RuO2test/A9/VI/D2
                    """)
                state = Annotations.Tags.load(dir)
                @test length(state.catalog) == 1
                @test state.catalog[1].name == "bad"
                @test state.catalog[1].color == (0xff, 0x30, 0x30)
                @test state.catalog[1].priority == 100
                @test state.assignments["RuO2test/A9/VI/D1"] == Set(["bad"])
                @test state.assignments["RuO2test/A9/VI/D2"] == Set(["bad"])
                @test isempty(state.measurement_assignments)
                @test !isfile(Annotations.Tags.tags_path(dir))
            end
        end

        @testset "legacy bad_measurements migration — measurement lines" begin
            mktempdir() do dir
                write(joinpath(dir, "bad_measurements"),
                    """
                    # comment
                    device RuO2test/A9/VI/D1
                    device RuO2test/A9/VI/D2
                    measurement abc123hash
                    """)
                state = Annotations.Tags.load(dir)
                @test state.assignments["RuO2test/A9/VI/D1"] == Set(["bad"])
                @test state.assignments["RuO2test/A9/VI/D2"] == Set(["bad"])
                @test Annotations.Tags.assigned_to_measurement(state, "abc123hash") ==
                    Set(["bad"])
                @test length(state.catalog) == 1
                @test state.catalog[1].name == "bad"
                @test !isfile(Annotations.Tags.tags_path(dir))
            end
        end

        @testset "combined tags.txt + bad_measurements" begin
            mktempdir() do dir
                write(joinpath(dir, "tags.txt"),
                    "[catalog]\ntodo\t30c0ff\t50\n\n[assignments]\ndevice\tRuO2test/A10/VI\ttodo\n")
                write(joinpath(dir, "bad_measurements"),
                    "device RuO2test/A9/VI/D1\nmeasurement abc123hash\n")
                state = Annotations.Tags.load(dir)

                names = Set(t.name for t in state.catalog)
                @test "todo" in names
                @test "bad" in names

                @test "todo" in state.assignments["RuO2test/A10/VI"]
                @test "bad" in state.assignments["RuO2test/A9/VI/D1"]
                @test Annotations.Tags.assigned_to_measurement(state, "abc123hash") ==
                    Set(["bad"])
            end
        end

        @testset "idempotency: load → save → load → save" begin
            mktempdir() do dir
                write(joinpath(dir, "tags.txt"),
                    "[catalog]\ntodo\t30c0ff\t50\n\n[assignments]\ndevice\tRuO2test/A10/VI\ttodo\n")
                write(joinpath(dir, "bad_measurements"),
                    "device RuO2test/A9/VI/D1\nmeasurement abc123hash\n")

                state1 = Annotations.Tags.load(dir)
                Annotations.Tags.save(dir, state1)
                content1 = read(Annotations.Tags.tags_path(dir), String)

                state2 = Annotations.Tags.load(dir)
                Annotations.Tags.save(dir, state2)
                content2 = read(Annotations.Tags.tags_path(dir), String)

                @test content1 == content2
            end
        end

        @testset "assigned_to_measurement" begin
            mktempdir() do dir
                write(joinpath(dir, "bad_measurements"), "measurement abc123\n")
                state = Annotations.Tags.load(dir)
                @test Annotations.Tags.assigned_to_measurement(state, "abc123") ==
                    Set(["bad"])
                @test isempty(Annotations.Tags.assigned_to_measurement(state, "nothere"))
                @test state.catalog[1].name == "bad"
            end
        end

        @testset "measurement_assignments in TagState constructor" begin
            state = Annotations.Tags.TagState(
                Annotations.Tags.TagDef[],
                Dict{String,Set{String}}("ChipA" => Set(["bad"])),
                Dict{String,Set{String}}("meas42" => Set(["bad"])),
            )
            @test Annotations.Tags.assigned_to_measurement(state, "meas42") == Set(["bad"])
            @test isempty(Annotations.Tags.assigned_to_measurement(state, "other"))
        end

        @testset "unknown assignment kind raises error" begin
            mktempdir() do dir
                write(joinpath(dir, "tags.txt"),
                    "[catalog]\nbad\tff3030\t100\n\n[assignments]\nunknown\tsome/path\tbad\n")
                @test_throws Annotations.Tags.TagsParseError Annotations.Tags.load(dir)
            end
        end

        @testset "effective inheritance" begin
            state = Annotations.Tags.TagState(
                Annotations.Tags.TagDef[],
                Dict{String,Set{String}}(
                    "ChipA" => Set(["bad"]),
                    "ChipA/SiteVI/D1" => Set(["todo"]),
                ),
                Dict{String,Set{String}}(),
            )
            eff = Annotations.Tags.effective(state, "ChipA/SiteVI/D1",
                ["ChipA", "ChipA/SiteVI"])
            @test eff == Set(["bad", "todo"])

            eff_no_anc = Annotations.Tags.effective(state, "ChipA/SiteVI/D1", String[])
            @test eff_no_anc == Set(["todo"])
        end

        @testset "dominant_color priority" begin
            catalog = [
                Annotations.Tags.TagDef("low", (0x10, 0x10, 0x10), 1),
                Annotations.Tags.TagDef("hi", (0xff, 0x00, 0x00), 100),
                Annotations.Tags.TagDef("mid", (0x80, 0x80, 0x80), 50),
            ]
            state = Annotations.Tags.TagState(catalog, Dict{String,Set{String}}(),
                Dict{String,Set{String}}())
            @test Annotations.Tags.dominant_color(state, Set(["low", "hi", "mid"])) ==
                (0xff, 0x00, 0x00)
            @test Annotations.Tags.dominant_color(state, Set(["low", "mid"])) ==
                (0x80, 0x80, 0x80)
            @test Annotations.Tags.dominant_color(state, Set(["unknown"])) === nothing
            @test Annotations.Tags.dominant_color(state, Set{String}()) === nothing
        end
    end

    @testset "Notes" begin
        @testset "round-trip with bracketed body content" begin
            fixture = joinpath(ANNOT_FIXTURES, "notes_basic")
            chipB = Annotations.Notes.read_section(fixture, "ChipB")
            @test occursin("Oxygen flow", chipB)
            site = Annotations.Notes.read_section(fixture, "ChipB/SiteVI")
            @test occursin("[observed in second pass]", site)
            @test occursin("dust particle", site)
            @test Annotations.Notes.read_section(fixture, "Missing/Path") == ""
        end

        @testset "merged_view ordering and editable flags" begin
            fixture = joinpath(ANNOT_FIXTURES, "notes_basic")
            view = Annotations.Notes.merged_view(fixture, "ChipB/SiteVI", ["ChipB"])
            @test length(view) == 2
            @test view[1].path == "ChipB"
            @test view[1].editable == false
            @test view[2].path == "ChipB/SiteVI"
            @test view[2].editable == true
        end

        @testset "missing notes file" begin
            mktempdir() do dir
                @test Annotations.Notes.read_section(dir, "Anything") == ""
                view = Annotations.Notes.merged_view(dir, "ChipA", ["Root"])
                @test length(view) == 1
                @test view[1].path == "ChipA"
                @test view[1].body == ""
                @test view[1].editable == true
            end
        end

        @testset "write_section! preserves and inserts" begin
            mktempdir() do dir
                src = joinpath(ANNOT_FIXTURES, "notes_basic", "notes.txt")
                cp(src, joinpath(dir, "notes.txt"))

                Annotations.Notes.write_section!(dir, "ChipB", "Updated body\n")
                @test Annotations.Notes.read_section(dir, "ChipB") == "Updated body"
                # other section preserved
                @test occursin("dust particle",
                    Annotations.Notes.read_section(dir, "ChipB/SiteVI"))

                # New section appended
                Annotations.Notes.write_section!(dir, "ChipC", "fresh note")
                @test Annotations.Notes.read_section(dir, "ChipC") == "fresh note"
                @test Annotations.Notes.read_section(dir, "ChipB") == "Updated body"

                # Order: original two then ChipC at end
                contents = read(joinpath(dir, "notes.txt"), String)
                @test findfirst("[ChipC]", contents).start >
                    findfirst("[ChipB/SiteVI]", contents).start
            end
        end
    end

end
