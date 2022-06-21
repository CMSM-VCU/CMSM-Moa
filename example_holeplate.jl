# START
include("Moa.jl")
include("Visualize.jl")
Moa.parse_input("example_holeplate.toml")
dt = 9999.
for node in Moa.nodes
    node.stableMass = Moa.Nodes.stableMass(node, dt, Moa.horizon, Moa.gridspacing)
end

# Increment SL stage:

for i in 1:100
    for bc in Moa.boundaryConditions
        if bc isa Moa.BoundaryConditions.StagedLoadingBC
            Moa.BoundaryConditions.apply_bc(bc)
        end
    end
end
#  Relax setup
using Profile
Profile.clear()
# kehis = Vector{Float64}()
print(Moa.MoaUtil.GetDamageVector(Moa.nodes))
Moa.TimeIntegration.adr(Moa.nodes, Moa.bonds, Moa.boundaryConditions, 9999., true)

## Relax
Moa.TimeIntegration.relax(Moa.nodes, Moa.bonds, Moa.boundaryConditions, 1.0e-20)



## PYVISTA VISUALIZE
plot(Moa.nodes, 2.0)
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

    while !anyBrokenBonds
        Threads.@threads for bond in Moa.bonds
            if !bond.isBroken && Moa.Bonds.should_break(bond)
                Moa.Bonds.break!(bond)
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


