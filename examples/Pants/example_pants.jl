# START
include("../../Moa.jl")
include("../../Visualize.jl")
Moa.parse_input("examples/Pants/example_pants.toml")

# need to put this into a setup function or have it at the end of parse_input
for node in Moa.nodes
    node.stableMass = Moa.Nodes.stableMass(node, 9999., Moa.horizon, Moa.gridspacing)
end

## Split pants
# repeating for now because it doesnt remove everthing in the first pass for an unknow reason
for i in 1:10
    for bond in Moa.bonds
        if (bond.from.position[1] >= 0 && bond.to.position[1] <= 0) ||
            (bond.from.position[1] <= 0 && bond.to.position[1] >= 0)
            if bond.from.position[2] > 0 || bond.to.position[2] > 0
                Moa.Bonds.delete(bond)
                # Moa.Bonds.break!(bond)
                deleteat!(Moa.bonds, findfirst(x->x==bond, Moa.bonds))
            end
        end
    end
end

## Use if you want to have a large initial tab movement
for bc in Moa.boundaryConditions
    if bc isa Moa.BoundaryConditions.StagedLoadingBC
        Moa.BoundaryConditions.apply_bc(bc)
    end
end


## LOOP
# Moa.TimeIntegration.relax(Moa.nodes,Moa.bonds,Moa.boundaryConditions, 1.0e-15)
for i in 1:10
    println("Staged loading step ", i)
    Moa.TimeIntegration.stagedloading(Moa.nodes,Moa.bonds,Moa.boundaryConditions, 1.0e-15)
end

## PYVISTA VISUALIZE
plotDisplacement(Moa.nodes, 1.0, 2)
plotDamage(Moa.nodes, 2.0)

plotDisplacement(Moa.nodes, 10.0, 1)
plotDamage(Moa.nodes, 1.0)

## WRITE OUTPUT
Moa.write_output("example_pants.csv", Moa.nodes)
