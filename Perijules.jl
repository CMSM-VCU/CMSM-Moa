import Pkg
Pkg.activate(".")

import NearestNeighbors


# Global (for now) varibles
grid_spacing = 2.
horizon = grid_spacing * 3.014


module PD
    include("Materials.jl")
    include("Nodes.jl")
    include("Bonds.jl")
end


# TOOD for tensile bundle:
# Nonlinear force response
# Force plane
# Material interface
# Bond failure
# Staged loading BC


## Initialize grid
# Create Materials
materials = Vector{PD.Material}()
push!(materials, PD.Material(1, 1., 0.2, 75))
push!(materials, PD.Material(2, 1., 0.2, 76))
println("Created ", length(materials), " materials!")

# Create Nodes (this process allows the nodes and position vector of vectors to point to the same thing)
positions = Vector{MVector{3, Float64}}()
push!(positions, [0.0, 0.0, 0.0])
push!(positions, [1.0, 0.0, 0.0])

nodes = Vector{PD.Node}()
push!(nodes, PD.Node(positions[1],materials[1],grid_spacing))
push!(nodes, PD.Node(positions[2],materials[2],grid_spacing))
println("Created ", length(nodes), " nodes!")


# Create Bonds (double bonded)
bonds = Vector{PD.Bond}()
tree = NearestNeighbors.BallTree(Vector{SVector{3, Float64}}(positions))
results = NearestNeighbors.inrange(tree, positions, horizon, false)
for (i, result) in enumerate(results)
    for neighbor in result
        if i != neighbor
            push!(bonds, PD.Bond(nodes[i], nodes[neighbor], false))
        end
    end
end
println("Created ", length(bonds), " bonds!")






## Time iteration
node_a_position = Vector{Float64}() # used to record a node's dispalcement
node_b_position = Vector{Float64}() # used to record a node's dispalcement

nodes[2].displacement[1] = 0.1

dt = grid_spacing / sqrt(maximum([material.bond_constant/material.density for material in materials])) * 0.01
for time_step in 1:10000
    for node in nodes
        # Average velocity for the time step assuming constant acceleration
        node.velocity = node.velocity + (dt*0.5)*(node.force/(node.volume*node.material.density))

        # Update displacement with average velocity
        node.displacement += node.velocity * dt
        node.force *= 0
    end

    # Apply displacement BCs
    nodes[1].displacement[1] = 0.

    # Calculate forces
    for bond in bonds
        # Eventually replace node.force with atomic floats and use threading for this loop
        PD.apply_force(bond)
    end

    for node in nodes
        # Calculate final velocity
        node.velocity = node.velocity + dt*0.5 * (node.force / (node.volume * node.material.density))
    end

    push!(node_a_position, nodes[1].position[1] + nodes[1].displacement[1])
    push!(node_b_position, nodes[2].position[1] + nodes[2].displacement[1])
end



# Plot the two points' deformed position
import Plots
Plots.plot([node_a_position,node_b_position])