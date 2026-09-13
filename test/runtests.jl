using Test
using DataBrowserAPI
using DataBrowserCache
using DataBrowserCore
using DataBrowserGUI
using DataBrowserSources

const TEST_CACHE_DEPOT = mktempdir()
pushfirst!(DEPOT_PATH, TEST_CACHE_DEPOT)
atexit(() -> rm(TEST_CACHE_DEPOT; force=true, recursive=true))

include(joinpath(@__DIR__, "test_project.jl"))

const UNIT_TEST_FILES = (
    "test_project_view_state.jl",
    "test_work_graph.jl",
    "test_collection_cache.jl",
    "test_payload_reconstruction.jl",
    "test_hierarchy_edit.jl",
    "test_open_options.jl",
    "test_label_cache.jl",
    "test_public_api.jl",
    "test_wide_cache.jl",
    "test_collection_identity.jl",
    "test_collection_index.jl",
    "test_collection_id_persistence.jl",
    "test_type_api.jl",
    "test_table_inspector.jl",
    "test_debug_tools.jl",
    "test_timing.jl",
    "test_gui_extensions.jl",
    "test_plots_extension_boundary.jl",
    "test_annotations.jl",
)

@testset verbose = true "package tests" begin
    for file in UNIT_TEST_FILES
        @testset "$file" begin
            include(joinpath(@__DIR__, file))
        end
    end
end

include(joinpath(@__DIR__, "..", "bench", "run.jl"))

@testset "performance snapshot" begin
    redirect_stdout(run_performance, devnull)
    @test isfile(STATUS_PATH) && filesize(STATUS_PATH) > 0
    println("Wrote ", STATUS_PATH)
end
