using TOML
using CSV

function parse_input(path::String)
    input = TOML.parsefile(path)
    
    println("Parsing input...")

    # Parse global scalars
    @show global gridspacing = input["gridspacing"]
    @show global horizon = input["horizon"]

    # Parse materials
    global materials = Vector{Materials.AbstractMaterial}()

    global defaultMaterial = Materials.parse_material(input["MaterialDefault"])
    push!(materials, defaultMaterial)

    if haskey(input, "Material")
        for materialInput in input["Material"]
            mat::Materials.AbstractMaterial = Materials.parse_material(materialInput)
            # Each material should have a unique id
            @assert !(mat.id in [material.id for material in materials])
            push!(materials, mat)
        end
    end

    # Parse grids
    global nodes = Vector{Nodes.AbstractNode}()
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
    for node in nodes
        for other in ProximitySearch.sample_cell_list(cell_list, node, horizon)
            push!(bonds, Bonds.Bond(node, other, false))
        end
    end
    println("Created ", length(bonds), " bonds")


    # Parse BCs
    global boundaryConditions = Vector{BoundaryConditions.AbstractBoundaryCondition}()
    for bc in input["BC"]
        push!(boundaryConditions, BoundaryConditions.parse_bc(bc, nodes))
        # println(typeof(last(boundaryConditions)))
    end
    println("Created ", length(boundaryConditions), " boundary conditions")


    # Other (force planes, etc.)
    global forceProbes = Vector{ForceProbes.AbstractForceProbe}()
    for forceProbe in input["ForceProbe"]
        push!(forceProbes, ForceProbes.parse_force_probe_plane(forceProbe, bonds))
    end
    println("Created ", length(forceProbes), " force probes")

    println("Finished reading input!")
end