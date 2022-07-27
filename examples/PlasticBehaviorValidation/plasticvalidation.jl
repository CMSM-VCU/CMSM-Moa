# START
using Plots
include("../../Moa.jl")
state = Moa.parse_input("examples/PlasticBehaviorValidation/plasticvalidation.toml");

# Testing displacements
his_strain = []
his_force = []


for k in 0:0.01:0.1
    for i in -0.01:0.001:k
        state.nodes[2].displacement = [i,0,0]
        push!(his_strain, i)
        push!(his_force, Moa.ForceProbes.measure_force(state.forceProbes[1]))
    end

    for i in k:-0.001:-0.01
        state.nodes[2].displacement = [i,0,0]
        push!(his_strain, i)
        push!(his_force, Moa.ForceProbes.measure_force(state.forceProbes[1]))
    end
end

plot(state.materials[1].customForceResponse[:,1], state.materials[1].customForceResponse[:,2], seriestype = :scatter)
plot!(his_strain, his_force)

##

@gif for i ∈ 1:length(his_strain)
    plot(state.materials[1].customForceResponse[:,1], state.materials[1].customForceResponse[:,2], seriestype = :scatter)
    plot!(his_strain[1:i], his_force[1:i])
    plot!([his_strain[i]], [his_force[i]], seriestype=:scatter)
end