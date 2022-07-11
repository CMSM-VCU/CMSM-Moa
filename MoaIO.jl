using TOML
using CSV
using .Nodes
using PyCall

function parse_input(path::String)
    println("Parsing input...")
    input = TOML.parsefile(path)

    # Parse global scalars
    gridspacing = Float64(input["gridspacing"])
    horizon = Float64(input["horizon"])

    # Parse materials
    materials = Vector{Materials.AMaterial}()
    defaultMaterial = Materials.parse_material(input["MaterialDefault"])
    push!(materials, defaultMaterial)

    if haskey(input, "Material")
        for materialInput in input["Material"]
            
            mat::Materials.AMaterial = Materials.parse_material(materialInput)
            
            # Each material should have a unique id
            @assert !(mat.id in [material.id for material in materials])
            push!(materials, mat)
        end
    end

    # Stable timestep
    dt = 0
    if haskey(input, "dt")
        dt = Float64(input["dt"])
    else
        dt = gridspacing / sqrt(maximum([mat.emod for mat in materials]) / minimum([mat.density for mat in materials]))
    end

    if dt <= 0
        throw(Exception)
    end


    # Parse grids
    nodes = Vector{Nodes.Node}()
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
            n ::Nodes.Node = Nodes.Node(
                                        Float64(row[:x]),
                                        Float64(row[:y]),
                                        Float64(row[:z]),
                                        mat,
                                        gridspacing
                                        )
            if :vx in input_grid.names && :vx in input_grid.names && :vx in input_grid.names
                n.velocity = [Float64(row[:vx]), Float64(row[:vy]), Float64(row[:vz])]
            end
            push!(nodes, 
                    n
                )
        end
    end

    # Stable mass for adr
    Threads.@threads for node in nodes
        node.stableMass = Nodes.stableMass(node, 9999., horizon, gridspacing)
    end

    # Create bonds
    println("Creating bonds...")
    bonds = Vector{AbstractTypes.ABond}()
    cell_list = ProximitySearch.create_cell_list(nodes, horizon)
    Threads.@threads for node in nodes
        for other in ProximitySearch.sample_cell_list(cell_list, node, horizon)
            b::Bonds.Bond = Bonds.Bond(node, other)
            push!(node.family, b)
        end
    end
    bonds = vcat([node.family for node in nodes]...)
    println("Created ", length(bonds), " bonds")


    # Parse BCs
    boundaryConditions = Vector{AbstractTypes.ABoundaryCondition}()
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
    forceProbes = Vector{AbstractTypes.AForceProbe}()
    if haskey(input, "ForceProbe")
        for forceProbe in input["ForceProbe"]
            push!(forceProbes, ForceProbes.parse_force_probe_plane(forceProbe, bonds))
        end
    end
    println("Created ", length(forceProbes), " force probes")

    println("Finished parsing input!")
    return Moa.state(gridspacing, horizon, dt, nodes, bonds, materials, boundaryConditions, forceProbes)
end

function write_output(path::String, state, timestep::Int64)
    pd = PyCall.pyimport_conda("pandas", "pandas")
    df = pd.DataFrame(
            data = [(
                node.position.x,
                node.position.y,
                node.position.z,
                node.displacement.x,
                node.displacement.y,
                node.displacement.z,
                node.velocity.x,
                node.velocity.y,
                node.velocity.z,
                node.material.id,
                Nodes.damage(node),
                Nodes.interfaceDamage(node),
                Nodes.materialDamage(node)
                ) for node in state.nodes],
            columns = ["x", "y", "z", "ux", "uy", "uz", "vx", "vy", "vz", "dmg", "dmgi", "dmgm", "mat"]
    )

    df.to_hdf(path, "t"*lpad(timestep, 7, "0"), mode="a")

end