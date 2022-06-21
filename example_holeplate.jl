# START
include("Moa.jl")
Moa.parse_input("example_holeplate.toml")
dt = 9999.
for node in Moa.nodes
    node.stableMass = Moa.Nodes.stableMass(node, dt, Moa.horizon, Moa.gridspacing)
end

# Increment SL stage:

for bc in Moa.boundaryConditions
    if bc isa Moa.BoundaryConditions.StagedLoadingBC
        Moa.BoundaryConditions.apply_bc(bc)
    end
end

#  Relax setup
using Profile
Profile.clear()
# kehis = Vector{Float64}()
Moa.TimeIntegration.adr(Moa.nodes, Moa.bonds, Moa.boundaryConditions, 9999., true)

# Relax


for i in 1:100
    # @profview Moa.TimeIntegration.adr(Moa.nodes, Moa.bonds, Moa.boundaryConditions, 9999., false)
    @profile Moa.TimeIntegration.adr(Moa.nodes, Moa.bonds, Moa.boundaryConditions, 1., false)
    println(i)
    # push!(kehis, Moa.MoaUtil.KineticEnergy(Moa.nodes))
end
view_profile()

# count = 1
# print(count)
# while Moa.MoaUtil.KineticEnergy(Moa.nodes) > 0.00001
#     Moa.TimeIntegration.adr(Moa.nodes, Moa.bonds, Moa.boundaryConditions, 9999., false)
#     push!(kehis, Moa.MoaUtil.KineticEnergy(Moa.nodes))
#     print("\r", count)
# end




## PYVISTA VISUALIZE

include("Visualize.jl")
plot(Moa.nodes)

##
using Plots
plot(kehis)


## WRITE OUTPUT

Moa.write_output("output.csv", Moa.nodes)


## Increment

for timestep in 1:300
    println("Displacing to ", newdisp)
    # Move tabs
    for bc in Moa.boundaryConditions
        if bc isa Moa.BoundaryConditions.StagedLoadingBC
            Moa.BoundaryConditions.apply_bc(bc)
        end
    end

    # Relax
    Moa.TimeIntegration.adr(Moa.nodes, Moa.bonds, Moa.boundaryConditions, 9999., true)
    for i in 1:10
        Moa.TimeIntegration.adr(Moa.nodes, Moa.bonds, Moa.boundaryConditions, 9999., false)
    end
    
    while Moa.MoaUtil.KineticEnergy(Moa.nodes) > 0.00001
        Moa.TimeIntegration.adr(Moa.nodes, Moa.bonds, Moa.boundaryConditions, 9999., false)
    end

    anyBrokenBonds = false
    while !anyBrokenBonds
        Threads.@threads for bond in Moa.bonds
            if !bond.isBroken && Moa.Bonds.should_break(bond)
                Moa.Bonds.break!(bond)
                # vvv THIS MAY BE A MISTAKE HERE (ACCESSING AND WRITING TO SHARED VARIABLE)
                anyBrokenBonds = true
            end
        end

        # Relax
        Moa.TimeIntegration.adr(Moa.nodes, Moa.bonds, Moa.boundaryConditions, 9999., true)
        for i in 1:10
            Moa.TimeIntegration.adr(Moa.nodes, Moa.bonds, Moa.boundaryConditions, 9999., false)
        end
        while Moa.MoaUtil.KineticEnergy(Moa.nodes) > 0.00001
            Moa.TimeIntegration.adr(Moa.nodes, Moa.bonds, Moa.boundaryConditions, 9999., false)
        end
    end
    # Repeat
end


