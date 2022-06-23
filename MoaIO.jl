using TOML
using CSV
using .Nodes

function parse_input(path::String)
    input = TOML.parsefile(path)
    
    println("Parsing input...")

    # Parse global scalars
    @show global gridspacing = input["gridspacing"]
    @show global horizon = input["horizon"]

    # Parse materials
    global materials = Vector{Materials.AMaterial}()

    global defaultMaterial = Materials.parse_material(input["MaterialDefault"])
    push!(materials, defaultMaterial)

    if haskey(input, "Material")
        for materialInput in input["Material"]
            mat::Materials.AMaterial = Materials.parse_material(materialInput)
            # Each material should have a unique id
            @assert !(mat.id in [material.id for material in materials])
            push!(materials, mat)
        end
    end

    # Parse grids
    global nodes = Vector{Nodes.Node}()
    for input_grid_object in input["Grid"]
        input_grid = CSV.File(input_grid_object["path"], stripwhitespace=true,  comment="#")
        # Required columns in grid file
        @assert :x in input_grid.names
        @assert :y in input_grid.names
        @assert :z in input_grid.names
        @assert :material in input_grid.names
        for row in input_grid
            mat = defaultMaterial
            material_candidates = [material for material in materials if material.id == row[:material]]
            if size(material_candidates)[1] != 0
                mat = first(material_candidates)
            else
                # Create copy of default material with new material number
                mat = copy(defaultMaterial)
                mat.id = row[:material]
                push!(materials, mat)
            end
            push!(nodes, 
                    Nodes.Node(
                        Float64(row[:x]),
                        Float64(row[:y]),
                        Float64(row[:z]),
                        mat,
                        gridspacing
                    )
                )
        end
    end

    # Create bonds
    global bonds = Vector{AbstractTypes.ABond}()
    cell_list = ProximitySearch.create_cell_list(nodes, horizon)
    println("CREATING BONDS...")
    parallel_bond_lists = fill(Vector{AbstractTypes.ABond}(), Threads.nthreads())

    Threads.@threads for node in nodes
        for other in ProximitySearch.sample_cell_list(cell_list, node, horizon)
            b::Bonds.Bond = Bonds.Bond(node, other)
            push!(node.family, b)
        end
        append!(parallel_bond_lists[Threads.threadid()], node.family)
    end
    bonds = vcat(parallel_bond_lists...)
    println("Created ", length(bonds), " bonds")

    # Parse BCs
    global boundaryConditions = Vector{BoundaryConditions.AbstractBoundaryCondition}()
    if haskey(input, "BC")
        for bc in input["BC"]
            push!(boundaryConditions, BoundaryConditions.parse_bc(bc, nodes))
        end
    end
    println("Created ", length(boundaryConditions), " boundary conditions")

    # NoFail regions TEMPORARY FIX!!!
    if haskey(input, "NoFail")
        for nofail in input["NoFail"]
            @assert haskey(nofail, "volume")
            for node in BoundaryConditions.get_nodes_within_volume(nodes, nofail["volume"])
                node.allowFailure = false
            end
        end
    end

    # Force probes
    global forceProbes = Vector{ForceProbes.AbstractForceProbe}()
    if haskey(input, "ForceProbe")
        for forceProbe in input["ForceProbe"]
            push!(forceProbes, ForceProbes.parse_force_probe_plane(forceProbe, bonds))
        end
    end
    println("Created ", length(forceProbes), " force probes")

    println("Finished reading input!")
end

function write_output(path::String, nodes::Vector{Nodes.Node})
    open(path, "w") do file
        write(file, "x,y,z,ux,uy,uz,dmg\n")
        for node in nodes
            write(file, string(node.position[1]) * ", " *
                        string(node.position[2]) * ", " *
                        string(node.position[3]) * ", " *
                        string(node.displacement[1]) * ", " *
                        string(node.displacement[2]) * ", " *
                        string(node.displacement[3]) * ", " *
                        string(Nodes.damage(node)) * "\n")
        end
    end
end