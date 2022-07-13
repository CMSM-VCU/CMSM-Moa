# START
include("../../Moa.jl")
include("../../Visualize.jl")
state = Moa.parse_input("examples/ContactTest/contacttest.toml");
println("DT: ", state.dt)
println("contact: ", state.contact)
println("contact distance: ", state.contactDistance)
println("contact coefficient: ", state.contactCoefficient)

println()

    # Move to before impact
state.dt = 0.736
Moa.TimeIntegration.dynamic_integration(state, 0.999)
# plotDisplacement(state, 1.0, 3)

## LOOP

state.dt = 0.02
for timestep in 1:10
    println("t",timestep)
    Moa.TimeIntegration.dynamic_integration(state, 0.999)
end
plotVelocityMagnitude(state, 1.0)


 
## PYVISTA VISUALIZE
plotDisplacement(Moa.nodes, 1000.0, 3)
plotDisplacementMagnitude(Moa.nodes, 1.0)
plotDamage(state, 1.0)
plotMaterialID(state, 1.0)
plotVelocity(Moa.nodes, 1.0, 3)
plotVelocityMagnitude(Moa.nodes, 10.0)

## WRITE OUTPUT
Moa.write_output("examples/ContactTest/output.h5", Moa.nodes, 3000)
