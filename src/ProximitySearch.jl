module ProximitySearch

using ..Nodes
using LinearAlgebra


struct CellList
    data::Array{Vector{Nodes.Node}, 3}
    division_size::Float64
    data_min::Vector{Float64}
    data_max::Vector{Float64}
end

function find_bounds(positions::Vector{Vector{Float64}})
    return [minimum([row[i] for row in positions]) for i in 1:3], [maximum([row[i] for row in positions]) for i in 1:3]
end


function create_cell_list(nodes::Vector{Nodes.Node}, radius::Float64, reference=false)
    if reference
        dim_min, dim_max = find_bounds([Vector{Float64}(node.position) for node in nodes])
    else
        dim_min, dim_max = find_bounds([Vector{Float64}(node.position + node.displacement) for node in nodes])
    end
    
    shape = (dim_max - dim_min) / radius

    size1,size2,size3 = Int64(ceil(shape[1]))+1, Int64(ceil(shape[2]))+1,Int64(ceil(shape[3]))+1
    # @debug "Creating cell list ($(size1), $(size2), $(size3))"

    # Each thread adds to it's own version
    # threaded_data = fill(fill(Vector{Nodes.Node}(),(size1,size2,size3)), Threads.nthreads());
    data = Array{Vector{Nodes.Node}}(undef,size1,size2,size3);
    for i in eachindex(data)
        data[i] = Vector{Nodes.Node}()
    end


    for node in nodes
        # Insert node
        pos = node.position
        if !reference
            pos += node.displacement
        end
        # index = []
        push!(data[Int64(ceil((pos[1] - dim_min[1]) / radius))+1,
                    Int64(ceil((pos[2] - dim_min[2]) / radius))+1,
                    Int64(ceil((pos[3] - dim_min[3]) / radius))+1], node)
    end

        
    # @debug "Created cell list: $(size(data))"
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