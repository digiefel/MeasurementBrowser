const PERFORMANCE_BROWSER = DataBrowserGUI.Browser

import CImGui

@testset "performance window sparklines" begin
    live = PERFORMANCE_BROWSER.LivePlotsState(capacity=3)
    timings = Dict{Symbol,Vector{Float64}}(
        :plot_load => [1.0, 2.0],
        :plot_setup => [3.0, 4.0],
        :plot_data => [5.0, 6.0],
        :plot_draw => [7.0, 8.0, 9.0],
    )

    PERFORMANCE_BROWSER._update_live_timings!(live, timings)
    @test live.load_buf == Float32[1, 2]
    @test live.total_buf == Float32[7, 8, 9]

    for value in (1.0f0, 2.0f0, 3.0f0, 4.0f0)
        PERFORMANCE_BROWSER._ring_push!(live.elapsed_buf, value, live.capacity)
    end
    @test live.elapsed_buf == Float32[2, 3, 4]

    # Pin the positional PlotLines signature `_render_sparkline!` relies on; CImGui's
    # PlotLines accepts no keyword arguments, so a kwarg call would MethodError at render.
    @test hasmethod(
        CImGui.PlotLines,
        Tuple{String,Vector{Float32},Int,Int,Ptr{Nothing},Float32,Float32,
              Tuple{Float32,Float32}},
    )

end
