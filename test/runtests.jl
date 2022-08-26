using Pkg
Pkg.add(path="S:/CMSM-Moa")
using Moa
using Test

@testset "Moa.jl" begin
    # Is the test set connected to Moa?
    @test Moa.HelloWorld()

    state = Moa.parse_input("test/twopoint/twopoint.toml")
    @test state ≠ Nothing
    @test length(state.nodes) == 2
    @test length(state.bonds) == 2


end
