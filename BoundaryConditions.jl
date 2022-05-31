module BoundaryConditions
using ..Materials
using ..Nodes
using ..Bonds

abstract type AbstractBoundaryCondition end

struct DisplacementBoundaryCondition <: AbstractBoundaryCondition
    nodes::Vector{Nodes.Node}
    displacement::Vector{Float64}
end

function apply_bc(bc::DisplacementBoundaryCondition)
    Threads.@threads for node in bc.nodes
        node.displacement = bc.displacement
    end
end



struct VelocityBoundaryCondition <: AbstractBoundaryCondition
    nodes::Vector{Nodes.Node}
    velocity::Vector{Float64}
end

function apply_bc(bc::VelocityBoundaryCondition)
    Threads.@threads for node in bc.nodes
        node.velocity = bc.velocity
    end
end


function get_nodes_within_volume(nodes, volume)
    # Push all relevant nodes
    n::Vector{Nodes.AbstractNode} = Vector{Nodes.AbstractNode}()
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


function parse_bc(inputDict, nodes)
    @assert haskey(inputDict, "type")
    @assert haskey(inputDict, "volume")

    relevant_nodes::Vector{Nodes.AbstractNode} = get_nodes_within_volume(nodes, inputDict["volume"])


    if inputDict["type"] == "Displacement"
        @assert haskey(inputDict, "displacement")
        return BoundaryConditions.DisplacementBoundaryCondition(relevant_nodes, inputDict["displacement"])
    elseif inputDict["type"] == "Velocity"
        @assert haskey(inputDict, "velocity")
        return BoundaryConditions.VelocityBoundaryCondition(relevant_nodes, inputDict["velocity"])
    else
        println("Unknown boundary condition type: ", inputDict["type"])
        @assert false
    end
end


end