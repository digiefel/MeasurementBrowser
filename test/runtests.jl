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

@testset "DataBrowser Tests" begin
    include(joinpath(@__DIR__, "test_project_view_state.jl"))
    include(joinpath(@__DIR__, "test_work_graph.jl"))
    include(joinpath(@__DIR__, "test_collection_cache.jl"))
    include(joinpath(@__DIR__, "test_hierarchy_edit.jl"))
    include(joinpath(@__DIR__, "test_open_options.jl"))
    include(joinpath(@__DIR__, "test_label_cache.jl"))
    include(joinpath(@__DIR__, "test_public_api.jl"))
    include(joinpath(@__DIR__, "test_wide_cache.jl"))
    include(joinpath(@__DIR__, "test_collection_identity.jl"))
    include(joinpath(@__DIR__, "test_collection_index.jl"))
    include(joinpath(@__DIR__, "test_collection_id_persistence.jl"))
    include(joinpath(@__DIR__, "test_type_api.jl"))
    include(joinpath(@__DIR__, "test_table_inspector.jl"))
    include(joinpath(@__DIR__, "test_debug_tools.jl"))
    include(joinpath(@__DIR__, "test_timing.jl"))
    include(joinpath(@__DIR__, "test_gui_extensions.jl"))
    include(joinpath(@__DIR__, "test_plots_extension_boundary.jl"))
    include(joinpath(@__DIR__, "test_annotations.jl"))
end

include(joinpath(@__DIR__, "..", "bench", "run.jl"))
run_performance()
