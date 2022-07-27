# START
include("../../Moa.jl")
state = Moa.parse_input("examples/StateTest/statetest.toml");

## LOOP
for bc in state.boundaryConditions
    if bc isa Moa.BoundaryConditions.StagedLoadingBC
        Moa.BoundaryConditions.apply_bc(bc)
    end
end

Moa.TimeIntegration.adr(state, true)
for i in 1:500
    # Moa.TimeIntegration.dynamic_integration_no_contact(state, .999)
    Moa.TimeIntegration.adr(state, false)
end
## PYVISTA VISUALIZE
include("../../Visualize.jl")
plotDisplacement(state, 1.0, 1)
plotVelocityMagnitude(state, 1.0)
plotVelocity(state, 1.0, 3)
plotDamage(state, 1.0)

## WRITE OUTPUT
Moa.write_output("output.csv", Moa.nodes, 1)