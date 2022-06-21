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
    dmg::Vector{Float64} = []
    for node in nodes
        numbonds::Int64 = length(node.family)
        numbrokenbonds:: Int64 = 0
        for bond in node.family
            if bond.isBroken
                numbrokenbonds += 1
            end
        end
        if numbonds == 0
            push!(dmg, -1.0)
        else
            push!(dmg, numbrokenbonds / numbonds)
        end
    end
    return dmg
end


end