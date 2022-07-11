# START
include("../../Moa.jl")
state = Moa.parse_input("examples/StateTest/statetest.toml");

## LOOP
for i in 1:5
    Moa.TimeIntegration.stagedloading(state, 1.0e-20)
end
## PYVISTA VISUALIZE
include("../../Visualize.jl")
plotDisplacement(state, 1.0, 1)
plotVelocityMagnitude(state, 1.0)
plotVelocity(state, 1.0, 3)
plotDamage(state, 1.0)

## WRITE OUTPUT
Moa.write_output("output.csv", Moa.nodes, 1)
