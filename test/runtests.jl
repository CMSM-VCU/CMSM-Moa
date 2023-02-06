using Revise
using Moa
using Test
using Logging
using Plots
println("\n\n\n\n\n")

Logging.disable_logging(LogLevel(Logging.Info))
# Logging.disable_logging(LogLevel(Logging.Debug))



@testset "BundleForceResponse" begin
    state = Moa.parse_input("test/bundle/twopointbundle.toml");
    push!(state.materials[1].stronglyConnected, 1)

    forces = []
    strains = []
    bondstrains = []

    push!(forces, Moa.ForceProbes.measure_force(state.forceProbes[1]))
    push!(strains, 0.0)

    for i in 1:20
        m_strain = i / 20 * 0.2
        for j in 1:1000
            strain = j/ 1000 * m_strain
            for node in state.nodes
                node.displacement[1] = node.position[1] * strain
            end
            for bond in state.bonds
                if !bond.isBroken && Moa.Bonds.shouldbreak(bond)
                    Moa.Bonds.break!(bond)
                end
            end
            push!(forces, Moa.ForceProbes.measure_force(state.forceProbes[1]))
            push!(strains, strain)
            
        end
        for j in 1:1000
            strain = (1000-j)/ 1000 * m_strain
            for node in state.nodes
                node.displacement[1] = node.position[1] * strain
            end
            for bond in state.bonds
                if !bond.isBroken && Moa.Bonds.shouldbreak(bond)
                    Moa.Bonds.break!(bond)
                end
            end
            push!(forces, Moa.ForceProbes.measure_force(state.forceProbes[1]))
            push!(strains, strain)
        end
    end

    @test isapprox(maximum(forces), 1422.912142068498 * tanh(20.059494696928038 * 0.038885466353955414) / 6, atol=0.5)
    # Bring these outside the testset
    plot(strains, forces, legend=false, markershape=:circle)
end
##




@testset "Connectivity" begin
    # Is the test set connected to Moa?
    @test typeof(Moa.version()) <: String
end

@testset "ProximitySearch.jl" begin
    # find_bounds
    @test Moa.ProximitySearch.find_bounds([[0.01, 1, 0], [-1, 1, 0], [-1, 0.5, -1], [-0.2, -0.5, 0.6], [0.7, -1, -0.2], [1,0,0], [0,0,1], [0,0,0]]) == ([-1,-1,-1],[1,1,1])
    @test Moa.ProximitySearch.find_bounds([[0.01, 1, 1.1], [2,1,0]]) == ([0.01,1,0],[2,1,1.1])
    for i in 1:100
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
    end

    # state = Moa.parse_input("test/twopoint/twopoint.toml")
    # print(Moa.ProximitySearch.find_bounds([Vector{Float64}(node.position) for node in state.nodes]))

end;

@testset "Dynamic Integration" begin
    state = Moa.parse_input("test/twopointdynamic/twopointdynamic.toml")
    state.dt = 10^-2
    x = []
    y = []
    for i in 1:250
        Moa.dynamic_integration(state, 1.0, false, false)
        push!(x,i)
        push!(y,state.nodes[2].displacement[1])
    end
    # plot(x,y)
    @test isapprox(minimum(y),-0.05, rtol=0.0001) 
    @test isapprox(maximum(y),0.05, rtol=0.0001) 
end

@testset "Adaptive Dynamic Relaxation" begin
    state = Moa.parse_input("test/twopointdynamic/twopointdynamic.toml")
    x = []
    y = []
    push!(x,0)
    push!(y,state.nodes[2].displacement[1])
    Moa.adr(state, true, false, 1.2e4)
    for i in 1:50
        Moa.adr(state, false, false, 1.2e4)
        push!(x,i)
        push!(y,state.nodes[2].displacement[1])
    end
    # plot(x,y)
    @test isapprox(minimum(y), 0.0, atol=10^-5)
    @test isapprox(y[end], 0.0, atol=10^-5)
    @test isapprox(maximum(y), 0.05, atol=0.001)
end

@testset "NodeMaterialReference" begin
    # Nodes should referernce a material in state.materials instead of having their own copy
    state = Moa.parse_input("test/twopoint/twopoint.toml")
    @test state.nodes[1].material ∈ values(state.materials)
end;

@testset "CNT" begin
    state = Moa.parse_input("test/tanhelastic/tanhelastic.toml")

    push!(state.materials[1].stronglyConnected, 1)
    forces = []
    strains = []
    bondstrains = []


    push!(forces, Moa.ForceProbes.measure_force(state.forceProbes[1]))
    push!(strains, 0.0)
    for i in 1:600
        strain = i / 2000
        Threads.@threads for node in state.nodes
            node.displacement[1] = node.position[1] * strain
        end
        Threads.@threads for bond in state.bonds
            if !bond.isBroken && Moa.Bonds.shouldbreak(bond)
                Moa.Bonds.break!(bond)
            end
        end

        push!(forces, Moa.ForceProbes.measure_force(state.forceProbes[1]))
        
        if strain < state.materials[1].critical_strain
            @test forces[end] ≈ state.materials[1].a * tanh(strain * state.materials[1].b) * 6
        else
            @test forces[end] ≈ 0.0
        end
        push!(strains, strain)
    end
    
    # println(Moa.ForceProbes.measure_force(state.forceProbes[1]))
end;
# plot(strains, forces, legend=false, markershape=:circle)
# plot(forces)

##
# state = Moa.parse_input("test/tanhelastic/twopointinterfacetest.toml");
# # Connect materials to themselves
# for (id, mat) in state.materials
#     append!(mat.stronglyConnected, id)
# end
# tensile = []
# shear = []
# append!(tensile, Moa.ForceProbes.measure_force(state.forceProbes[1]))
# append!(shear, Moa.ForceProbes.measure_force(state.forceProbes[2]))

# for ts in 1:2500
#     for bc in state.boundaryConditions
#         if bc isa Moa.BoundaryConditions.StagedLoadingBC
#             Moa.BoundaryConditions.increment_staged_loading(bc)
#         end
#     end

#     Threads.@threads for bond in state.bonds
#         if !bond.isBroken && Moa.Bonds.shouldbreak(bond)
#             Moa.Bonds.break!(bond)
#         end
#     end

#     Threads.@threads for bond in state.bonds
#         # If interface bond
#         if Moa.Bonds.getstrain(bond) < 0 && bond.from.material.id ∉ bond.to.material.stronglyConnected
#             Moa.Bonds.break!(bond)
#         end
#     end

#     append!(tensile, Moa.ForceProbes.measure_force(state.forceProbes[1]))
#     append!(shear, Moa.ForceProbes.measure_force(state.forceProbes[2]))
# end
# plot(tensile)
# plot!(shear, legend=false)
# plot(tensile[1:300])
# plot!(shear[1:300], legend=false)

##
@testset "Break Compressed Interface Bonds" begin
    forces = []
    displacements = []
    strains = []
    state = Moa.parse_input("test/tanhelastic/twopointtanh.toml")

    push!(forces, Moa.ForceProbes.measure_force(state.forceProbes[1]))
    push!(displacements, 0)
    push!(strains, 0)
    
    for x in 0:0.01:3*π
        state.nodes[2].displacement[1] = 0.9*sin(x)
        Threads.@threads for bond in state.bonds
            if Moa.Bonds.getstrain(bond) < 0 && bond.from.material.id ∉ bond.to.material.stronglyConnected
                Moa.Bonds.break!(bond)
            end
        end

        push!(displacements, 0.9*sin(x))
        push!(forces, Moa.ForceProbes.measure_force(state.forceProbes[1]))
        push!(strains, Moa.Bonds.getstrain(state.bonds[1]))
        if any(map(x->x<0, strains))
            @test forces[end] ≈ 0.0
        end
    end
    # plot(forces*0.1, label="Force")
    # plot!(displacements, label="Displacement")
    # plot(displacements, legend=false)
    # plot(strains)
    # plotDisplacement(state, 1.0, 1)
end;

@testset "Contact" begin
    state = Moa.parse_input("test/twopoint/twopointcontact.toml");
    state.nodes[2].displacement[1] = -0.1
    Moa.apply_contact_force(state)
    @test state.nodes[1].force[1] ≈ 0.5
    @test state.nodes[2].force[1] ≈ -0.5
end;


@testset "Moa.jl" begin
    state = Moa.parse_input("test/twopoint/twopoint.toml")
    @test state ≠ Nothing
    
    @test length(state.nodes) == 2

    # Bond connection (should change to 1 bond in the future)
    @test length(state.bonds) == 2
end;

