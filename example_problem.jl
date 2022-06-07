using Plots
include("Moa.jl")


Moa.parse_input("example_input.toml")
Moa.nodes[2].displacement[1] = 0.1
disphis_a = Vector{Float64}()
disphis_b = Vector{Float64}()
disphis_c = Vector{Float64}()
pforce = Vector{Float64}()
forces = Vector{Float64}()
print("0")

# dt = 0.02
dt = 9999.
for node in Moa.nodes
    node.stableMass = Moa.Nodes.stableMass(node, dt, Moa.horizon, Moa.gridspacing)
end
for i in 1:300
    print("\r",i)
    # Moa.TimeIntegration.dynamic_integration(Moa.nodes,Moa.bonds, Moa.boundaryConditions, 0.1)
    if i == 1
        Moa.TimeIntegration.adr(Moa.nodes,Moa.bonds, Moa.boundaryConditions, dt, true)
    else
        Moa.TimeIntegration.adr(Moa.nodes,Moa.bonds, Moa.boundaryConditions, dt, false)
    end
    push!(disphis_a, Moa.nodes[1].position[1] + Moa.nodes[1].displacement[1])
    push!(disphis_b, Moa.nodes[2].position[1] + Moa.nodes[2].displacement[1])
    push!(forces, Moa.ForceProbes.measure_force(Moa.forceProbes[1]))
    push!(pforce, Moa.nodes[2].force[1])
end



default(show=true)

plot(disphis_a)
plot!(disphis_b)
# plot(forces)
# plot(pforce)