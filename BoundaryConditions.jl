module BoundaryConditions
using LinearAlgebra
using ..Materials
using ..Nodes
using ..Bonds
using ..AbstractTypes: ABoundaryCondition


struct DisplacementBC <: ABoundaryCondition
    nodes::Vector{Nodes.Node}
    displacement::Vector{Float64}
end

struct VelocityBC <: ABoundaryCondition
    nodes::Vector{Nodes.Node}
    velocity::Vector{Float64}
end

mutable struct StagedLoadingBC <: ABoundaryCondition
    nodes::Vector{Nodes.Node}
    currentDisplacement::Vector{Float64}
    increment::Vector{Float64}
end

function apply_bc(bc::DisplacementBC)
    Threads.@threads for node in bc.nodes
        node.displacement = bc.displacement
        # Is this faster than creating a new vector of zeros?
        node.velocity *= 0
        @atomic node.force = zeros(3)
    end
end

function apply_bc(bc::VelocityBC)
    Threads.@threads for node in bc.nodes
        node.velocity = bc.velocity
    end
end

function apply_bc(bc::StagedLoadingBC)
    bc.currentDisplacement += bc.increment
    Threads.@threads for node in bc.nodes
        node.displacement += bc.increment
    end
end

"All nodes within the specified volume
(should move this function somewhere else and add more options)"
function get_nodes_within_volume(nodes, volume)
    # Push all relevant nodes
    n::Vector{Nodes.Node} = Vector{Nodes.Node}()
    for node in nodes
        if  node.position[1] > volume[1] &&
            node.position[1] < volume[2] &&
            node.position[2] > volume[3] &&
            node.position[2] < volume[4] &&
            node.position[3] > volume[5] &&
            node.position[3] < volume[6]
            push!(n, node)
        end
    end
    if size(n)[1] == 0
        println("WARNING, CAPTURED VOLUME WITH 0 NODES!!!!")
    end
    return n
end

"Creates boundary conditions from an input dictionary"
function parse_bc(inputDict, nodes)
    @assert haskey(inputDict, "type")
    @assert haskey(inputDict, "volume")

    relevant_nodes::Vector{Nodes.Node} = get_nodes_within_volume(nodes, inputDict["volume"])


    if inputDict["type"] == "Displacement"
        @assert haskey(inputDict, "displacement")
        return DisplacementBC(relevant_nodes, inputDict["displacement"])

    elseif inputDict["type"] == "Velocity"
        @assert haskey(inputDict, "velocity")
        return VelocityBC(relevant_nodes, inputDict["velocity"])

    elseif inputDict["type"] == "Staged Loading"
        @assert haskey(inputDict, "increment")
        return StagedLoadingBC(relevant_nodes, zeros(3), inputDict["increment"])
    else
        println("Unknown boundary condition type: ", inputDict["type"])
        @assert false
    end
end


end