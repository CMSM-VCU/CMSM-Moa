# START
include("../../Moa.jl")
include("../../Visualize.jl")
Moa.parse_input("examples/HolePlate/example_holeplate.toml")
dt = 9999.
for node in Moa.nodes
    node.stableMass = Moa.Nodes.stableMass(node, dt, Moa.horizon, Moa.gridspacing)
end


## LOOP

for bc in Moa.boundaryConditions
    if bc isa Moa.BoundaryConditions.StagedLoadingBC
        Moa.BoundaryConditions.apply_bc(bc)
    end
end

Moa.TimeIntegration.stagedloading(Moa.nodes,Moa.bonds,Moa.boundaryConditions, 1.0e-20)

## PYVISTA VISUALIZE
plotDisplacement(Moa.nodes, 1.0, 1)
plotDamage(Moa.nodes, 1.0)

## WRITE OUTPUT
Moa.write_output("output.csv", Moa.nodes)
