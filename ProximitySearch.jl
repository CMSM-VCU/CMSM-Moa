module ProximitySearch

using ..Nodes
using LinearAlgebra


struct CellList
    data::Array{Vector{Nodes.Node}, 3}
    division_size::Float64
    data_min::Vector{Float64}
    data_max::Vector{Float64}
end

function create_cell_list(nodes::Vector{Nodes.Node}, radius::Float64)
    # Can be optimized
    positions = [node.position + node.displacement for node in nodes]
    dim_max = copy(positions[1])
    dim_min = copy(positions[1])
    for position in positions
        if position[1] > dim_max[1]
            dim_max[1] = position[1]
        end
        if position[2] > dim_max[2]
            dim_max[2] = position[2]
        end
        if position[3] > dim_max[3]
            dim_max[3] = position[3]
        end
        if position[1] < dim_min[1]
            dim_min[1] = position[1]
        end
        if position[2] < dim_min[2]
            dim_min[2] = position[2]
        end
        if position[3] < dim_min[3]
            dim_min[3] = position[3]
        end
    end
    shape = (dim_max - dim_min) / radius
    size1,size2,size3 = Int64(ceil(shape[1]))+1, Int64(ceil(shape[2]))+1,Int64(ceil(shape[3]))+1
    println("Creating cell list (", size1, ", ", size2, ", ", size3, ")")
    data = Array{Vector{Nodes.Node}}(undef,size1,size2,size3);
    for i in eachindex(data)
        data[i] = Vector{Nodes.Node}()
    end


    for node in nodes
        # Insert node
        pos = node.position + node.displacement
        # index = []
        push!(data[Int64(ceil((pos[1] - dim_min[1]) / radius))+1,
                    Int64(ceil((pos[2] - dim_min[2]) / radius))+1,
                    Int64(ceil((pos[3] - dim_min[3]) / radius))+1], node)
    end

        
    println("Created cell list: ", size(data))
    return CellList(data, radius, dim_min, dim_max)
end

function sample_cell_list(cell_list::CellList, node::Nodes.Node, radius::Float64)
    pos = node.position + node.displacement
    index = [
        Int64(ceil((pos[1] - cell_list.data_min[1]) / radius))+1,
        Int64(ceil((pos[2] - cell_list.data_min[2]) / radius))+1,
        Int64(ceil((pos[3] - cell_list.data_min[3]) / radius))+1
        ]   
    shape = size(cell_list.data)
    result = Vector{Nodes.Node}()

    for i in -1:1
        if index[1]+i > 0 && index[1]+i <= shape[1]
            for j in -1:1
                if index[2]+j > 0 && index[2]+j <= shape[2]
                    for k in -1:1
                        if index[3]+k > 0 && index[3]+k <= shape[3]
                            append!(result, cell_list.data[index[1]+i, index[2]+j, index[3]+k])
                        end
                    end
                end
            end
        end
    end
    return [n for n in result if LinearAlgebra.norm(n.position+n.displacement-node.position-node.displacement) < radius && n != node]
end

end