module MoaUtil
using ..Nodes

function KineticEnergy(nodes::Vector{Nodes.Node})
    KEs = zeros(Threads.nthreads())
    Threads.@threads for node in nodes
        KEs[Threads.threadid()] += 0.5 * Nodes.mass(node) * sum(node.velocity.^2)
    end
    return sum(KEs)
end

function GetDamageVector(nodes::Vector{Nodes.Node})
    # Create array of zeros
    dmg::Vector{Float64} = fill(0.0, length(nodes))

    # Populate array
    Threads.@threads for i in 1:length(nodes)
        dmg[i] = Nodes.damage(nodes[i])
    end
    
    return dmg
end


end