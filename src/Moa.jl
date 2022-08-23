module Moa
include("AbstractTypes.jl")
include("Materials.jl")
include("Nodes.jl")
include("Bonds.jl")
include("ProximitySearch.jl")
include("BoundaryConditions.jl")
include("ForceProbes.jl")
include("MoaUtil.jl")


mutable struct state
    gridSpacing::Float64
    horizon::Float64
    dt::Float64

    nodes::Vector{Nodes.Node}
    bonds::Vector{AbstractTypes.ABond}
    materials::Dict{Int64, AbstractTypes.AMaterial}

    boundaryConditions::Vector{AbstractTypes.ABoundaryCondition}
    forceProbes::Vector{AbstractTypes.AForceProbe}

    contact::Bool
    contactDistance::Float64
    contactCoefficient::Float64
end

include("TimeIntegration.jl")

function HelloWorld()
    println("Hello from Moa!")
end

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
    materials = Dict{Int64, Materials.AMaterial}()
    defaultMaterial = Materials.parse_material(input["MaterialDefault"])
    materials[defaultMaterial.id] = defaultMaterial

    if haskey(input, "Material")
        for materialInput in input["Material"]
            
            mat::Materials.AMaterial = Materials.parse_material(materialInput)
            
            # Each material should have a unique id
            @assert !haskey(materials, mat.id)
            materials[mat.id] = mat
        end
    end

    # Stable timestep
    dt = 0
    if haskey(input, "dt")
        dt = Float64(input["dt"])
    else
        dt = gridspacing / sqrt(maximum([mat.emod for (id,mat) in materials]) / minimum([mat.density for (id,mat) in materials]))
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
            material_candidates = [material for (id,material) in materials if material.id == row[:material]]

            if size(material_candidates)[1] != 0
                mat = first(material_candidates)
            else
                # Create copy of default material with new material number
                mat = copy(defaultMaterial)
                mat.id = row[:material]
                materials[mat.id] = mat
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
    cell_list = ProximitySearch.create_cell_list_reference_configuration(nodes, horizon)
    for node in nodes
        for other in ProximitySearch.sample_cell_list(cell_list, node, horizon)
            neighborhood = collect(Iterators.flatten([(bond.from, bond.to) for bond in node.family]))
            if other ∉ neighborhood
                if node.material isa Materials.CustomPlastic
                    push!(node.family, Bonds.PlasticBond(node, other))
                else
                    push!(node.family, Bonds.Bond(node, other))
                end
            end
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

    # Disconnects
    println("# Bonds before disconnect(s): ", length(bonds))
    if haskey(input, "Disconnect")
        for dc in input["Disconnect"]
            @assert haskey(dc, "ids")

            # All the material ids to be disconnected from each other
            ids = Vector{Int64}(dc["ids"])
            
            # Get all the bonds that are between these materials (distinct)
            toRemove = [bond for bond in bonds if bond.from.material.id != bond.to.material.id &&
                                                    bond.from.material.id ∈ ids &&
                                                    bond.to.material.id ∈ ids]
            println("Disconnecting ", length(toRemove), " bonds")
                                                    
            # Remove bonds from families without effecting damage
            for bond in toRemove
                Bonds.delete(bond)
            end

            # Remove bonds from global list
            filter!(bond->bond∉toRemove, bonds)
        end
    end
    println("# Bonds after disconnect(s): ", length(bonds))

    # Force probes
    forceProbes = Vector{AbstractTypes.AForceProbe}()
    if haskey(input, "ForceProbe")
        for forceProbe in input["ForceProbe"]
            push!(forceProbes, ForceProbes.parse_force_probe_plane(forceProbe, bonds))
        end
    end
    println("Created ", length(forceProbes), " force probes")


    # Contact properties
    useContact = haskey(input, "Contact")
    contactDistance = 0.0
    contactCoefficient = 0.0
    if useContact
        @assert haskey(input["Contact"], "distance")
        @assert haskey(input["Contact"], "coefficient")
        contactDistance = input["Contact"]["distance"]
        contactCoefficient = input["Contact"]["coefficient"]
    end

    println("Finished parsing input!")
    return Moa.state(gridspacing, horizon, dt, nodes, bonds, materials, boundaryConditions, forceProbes, useContact, contactDistance, contactCoefficient)
end

function write_output(path::String, state, timestep::Int64)
    pd = PyCall.pyimport_conda("pandas", "pandas")

    # Node data
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


    # for bc in state.boundaryConditions
    #     if bc isa BoundaryConditions.StagedLoadingBC
    #         # Write current displacement
    #     end
    # end

    # for probe in state.forceProbes
    #     # Write ForceProbes.mreasure_forcce(probe)
    # end

end

end