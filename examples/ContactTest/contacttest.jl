# START
include("../../Moa.jl")
include("../../Visualize.jl")
Moa.parse_input("examples/ContactTest/contacttest.toml")


# Moa.TimeIntegration.dynamic_integration(Moa.nodes,Moa.bonds,Moa.boundaryConditions, 0.88)

Moa.TimeIntegration.dynamic_integration(Moa.nodes,Moa.bonds,Moa.boundaryConditions, 0.736)
plotDisplacement(Moa.nodes, 1.0, 3)

## LOOP
dt = 0.0003
for timestep in 1:3000
    println("t",timestep)
    Moa.TimeIntegration.dynamic_integration(Moa.nodes,Moa.bonds,Moa.boundaryConditions, dt)
end
plotVelocityMagnitude(Moa.nodes, 1.0)


 
## PYVISTA VISUALIZE
plotDisplacement(Moa.nodes, 1000.0, 3)
plotDisplacementMagnitude(Moa.nodes, 1.0)
plotDamage(Moa.nodes, 1.0)
plotMaterialID(Moa.nodes, 1.0)
plotVelocity(Moa.nodes, 1.0, 3)
plotVelocityMagnitude(Moa.nodes, 10.0)

## WRITE OUTPUT
Moa.write_output("examples/ContactTest/output.h5", Moa.nodes, 3000)
