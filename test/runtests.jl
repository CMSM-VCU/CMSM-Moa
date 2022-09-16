# using Pkg
# Pkg.add(path="D:/CMSM-Moa")
using Revise
using Moa
using Test
using Logging
println("\n\n\n\n\n")

Logging.disable_logging(LogLevel(0))
# Logging.disable_logging(LogLevel(-1000))

@testset "Connectivity" begin
    # Is the test set connected to Moa?
    @test typeof(Moa.version()) <: String
end

# RACE CONDITION IN PROXIMITY SEARCH!!!!!!!!!!!!!!!!!!
@testset "ProximitySearch.jl" begin
    # find_bounds
    @test Moa.ProximitySearch.find_bounds([[0.01, 1, 0], [-1, 1, 0], [-1, 0.5, -1], [-0.2, -0.5, 0.6], [0.7, -1, -0.2], [1,0,0], [0,0,1], [0,0,0]]) == ([-1,-1,-1],[1,1,1])
    @test Moa.ProximitySearch.find_bounds([[0.01, 1, 1.1], [2,1,0]]) == ([0.01,1,0],[2,1,1.1])

    state = Moa.parse_input("test/neighbortest/neighbortest.toml")
    @test length(state.nodes[1].family) == 1
    @test length(state.nodes[2].family) == 2
    @test length(state.nodes[3].family) == 1

    @test state.nodes[1] ∈ [bond.to for bond in state.nodes[2].family]
    @test state.nodes[2] ∈ [bond.to for bond in state.nodes[1].family]
    @test state.nodes[3] ∈ [bond.to for bond in state.nodes[2].family]
    @test state.nodes[2] ∈ [bond.to for bond in state.nodes[3].family]

    @test state.nodes[2] ∉ [bond.to for bond in state.nodes[2].family]
    @test state.nodes[1] ∉ [bond.to for bond in state.nodes[1].family]
    @test state.nodes[3] ∉ [bond.to for bond in state.nodes[3].family]

    @test state.nodes[1] ∉ [bond.to for bond in state.nodes[3].family]
    @test state.nodes[3] ∉ [bond.to for bond in state.nodes[1].family]

    # @test state.nodes[2] ∉ state.nodes[2].family
    # @test state.nodes[1] ∉ state.nodes[1].family
    # @test state.nodes[3] ∉ state.nodes[3].family

    # @test state.nodes[1] ∉ state.nodes[3].family
    # @test state.nodes[3] ∉ state.nodes[1].family

    # state = Moa.parse_input("test/twopoint/twopoint.toml")
    # print(Moa.ProximitySearch.find_bounds([Vector{Float64}(node.position) for node in state.nodes]))

end;



@testset "Moa.jl" begin
    state = Moa.parse_input("test/twopoint/twopoint.toml")
    @test state ≠ Nothing
    @test length(state.nodes) == 2
    @test length(state.bonds) == 2
end

;