module MoaUtil
using ..Nodes

function KineticEnergy(nodes::Vector{Nodes.Node})
    KEs = zeros(Threads.nthreads())
    Threads.@threads for node in nodes
        KEs[Threads.threadid()] += 0.5 * Nodes.mass(node) * sum(node.velocity.^2)
    end
    return sum(KEs)
end



end