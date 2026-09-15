using Test

files = isempty(ARGS) ? sort(filter(name -> startswith(name, "test_") && endswith(name, ".jl"), readdir(@__DIR__))) : ARGS
@testset "DataBrowserCore" begin
    for file in files
        include(joinpath(@__DIR__, file))
    end
end
