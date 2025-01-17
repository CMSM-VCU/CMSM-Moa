module Moa
using Logging
include("AbstractTypes.jl")
include("Materials/Materials.jl")
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

    output_folder::String
end

include("TimeIntegration.jl")

function version()::String
    return "0.0.0"
end

using TOML
using CSV
using .Nodes
using PyCall

function parse_input(path::String)
    @info "REMOVED COMPRESSIVE INTERFACE BOND FAILURE!!!"
    @info "Parsing input..."
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
            # Initial material
            mat = defaultMaterial
            material_candidates = [material for (id,material) in materials if material.id == row[:material]]
            if length(material_candidates)[1] > 1
                error("Multiple materials with the same material number")
            elseif length(material_candidates) == 1
                mat = material_candidates[1]
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
            # Initial velocities
            if :vx in input_grid.names
                n.velocity[1] = Float64(row[:vx])
            end
            if :vy in input_grid.names
                n.velocity[2] = Float64(row[:vy])
            end
            if :vz in input_grid.names
                n.velocity[3] = Float64(row[:vz])
            end

            # Initial displacements
            if :ux in input_grid.names
                n.displacement[1] = Float64(row[:ux])
            end
            if :vy in input_grid.names
                n.displacement[2] = Float64(row[:uy])
            end
            if :uz in input_grid.names
                n.displacement[3] = Float64(row[:uz])
            end

            push!(nodes,n)
        end
    end

    # Stable mass for adr
    Threads.@threads for node in nodes
        node.stableMass = Nodes.stableMass(node, 9999., horizon, gridspacing)
    end

    # Create bonds
    @info "Creating bonds..."
    bonds = Vector{AbstractTypes.ABond}()
    cell_list = ProximitySearch.create_cell_list(nodes, horizon, true)
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
    @info "Created $(length(bonds)) bonds"


    # Parse BCs
    boundaryConditions = Vector{AbstractTypes.ABoundaryCondition}()
    if haskey(input, "BC")
        for bc in input["BC"]
            push!(boundaryConditions, BoundaryConditions.parse_bc(bc, nodes))
        end
    end
    @info "Created $(length(boundaryConditions)) boundary conditions"

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
    @info "# Bonds before disconnect(s): $(length(bonds))"
    if haskey(input, "Disconnect")
        for dc in input["Disconnect"]
            @assert haskey(dc, "ids")

            # All the material ids to be disconnected from each other
            ids = Vector{Int64}(dc["ids"])
            
            # Get all the bonds that are between these materials (distinct)
            toRemove = [bond for bond in bonds if bond.from.material.id != bond.to.material.id &&
                                                    bond.from.material.id ∈ ids &&
                                                    bond.to.material.id ∈ ids]
            @info "Disconnecting $(length(toRemove)) bonds"
                                                    
            # Remove bonds from families without effecting damage
            for bond in toRemove
                Bonds.delete(bond)
            end

            # Remove bonds from global list
            filter!(bond->bond∉toRemove, bonds)
        end
    end
    @info "# Bonds after disconnect(s): $(length(bonds))"

    # Force probes
    forceProbes = Vector{AbstractTypes.AForceProbe}()
    if haskey(input, "ForceProbe")
        for forceProbe in input["ForceProbe"]
            push!(forceProbes, ForceProbes.parse_force_probe_plane(forceProbe, bonds))
        end
    end
    @info "Created $(length(forceProbes)) force probes"


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

    output_folder = ""
    if haskey(input, "output")
        output_folder = input["output"]["path"]
        if !isdir(output_folder)
            mkdir(output_folder)
        end
    end


    @info "Finished parsing input!"
    return Moa.state(gridspacing, horizon, dt, nodes, bonds, materials, boundaryConditions, forceProbes, useContact, contactDistance, contactCoefficient, output_folder)
end

function write_output(state, timestep::Int64)
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
                Nodes.materialDamage(node),
                # These two metrics are for BondBasedBundle Material only (could also work with BondBasedTanhElasatic), 
                # if this gives you problems, comment these two out or delete:
                length([bond for bond in node.family if bond.isBroken && bond.from.material.id ∈ bond.to.material.stronglyConnected]),
                length([bond for bond in node.family if bond.isBroken && bond.from.material.id ∉ bond.to.material.stronglyConnected])
                ) for node in state.nodes],
            # columns = ["x", "y", "z", "ux", "uy", "uz", "vx", "vy", "vz", "dmg", "dmgi", "dmgm", "mat"]
            columns = ["x", "y", "z", "ux", "uy", "uz", "vx", "vy", "vz", "mat", "dmg", "dmgi", "dmgm", "num_broken_material_bonds", "num_broken_interface_bonds"]
    )
    df.to_hdf(state.output_folder*"/output.h5", "t"*lpad(timestep, 7, "0"), mode="a")


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