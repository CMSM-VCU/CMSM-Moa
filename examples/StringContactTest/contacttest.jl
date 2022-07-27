# START
include("../../Moa.jl")
include("../../Visualize.jl")
state = Moa.parse_input("examples/StringContactTest/contacttest.toml");

## LOOP
for timestep in 1:20
    println("\nt",timestep)
    Moa.TimeIntegration.stagedloading(state, 1e-6, 10000)
    Moa.write_output("examples/StringContactTest/output.h5", state, timestep)
end
