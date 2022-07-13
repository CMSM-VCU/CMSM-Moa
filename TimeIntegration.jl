module TimeIntegration

using ..Materials
using ..Nodes
using ..Bonds
using ..BoundaryConditions
using ..AbstractTypes
using ..MoaUtil
using ..TimeIntegration
using ..ProximitySearch

using LinearAlgebra: norm


function apply_contact_force(state)
    if state.contact
        cell_list = ProximitySearch.create_cell_list(state.nodes, state.contactDistance)
        Threads.@threads for node in state.nodes
            for other in ProximitySearch.sample_cell_list(cell_list, node, state.contactDistance)
                # Dont add force if the node is a familiy member
                isInFamily::Bool = false
                for bond in node.family
                    if bond.to == other
                        isInFamily = true
                        break
                    end
                end
                if isInFamily
                    continue
                end
                
                distance::Float64 = Nodes.distance(node, other)
                direction::Vector{Float64} = (other.position - node.position) / norm(other.position - node.position)
                @atomic node.force += (state.contactCoefficient * ((state.contactDistance - distance) / state.contactDistance) * node.volume * other.volume) * direction
            end
        end
    end
end

function dynamic_integration(state, damping)

    Threads.@threads for node in state.nodes
        # Average velocity for the time step assuming constant acceleration
        node.velocity = node.velocity + (state.dt*0.5)*(node.force/(node.volume*node.material.density))
    end
    
    # Velocity BCs
    for bc in state.boundaryConditions
        if bc isa BoundaryConditions.VelocityBC
            BoundaryConditions.apply_bc(bc)
        end
    end
    
    Threads.@threads for node in state.nodes
        # Update displacement with average velocity
        node.displacement += node.velocity * state.dt

        # Zero out force
        @atomic node.force = zeros(3)
    end

    # Displacement BCs
    for bc in state.boundaryConditions
        if bc isa BoundaryConditions.DisplacementBC
            BoundaryConditions.apply_bc(bc)
        end
    end

    # Break bonds
    Threads.@threads for bond in state.bonds
        if !bond.isBroken && Bonds.should_break(bond)
            Bonds.break!(bond)
        end
    end

    # Apply bond force to nodes
    Threads.@threads for bond in state.bonds
        Bonds.apply_force(bond)
    end

    apply_contact_force(state)

    # Calculate final velocity
    Threads.@threads for node in state.nodes
        node.velocity = node.velocity*damping + state.dt*0.5 * (node.force / (node.volume * node.material.density))
    end

end

function dynamic_integration_no_contact(state, damping)

    Threads.@threads for node in state.nodes
        # Average velocity for the time step assuming constant acceleration
        node.velocity = node.velocity + (state.dt*0.5)*(node.force/(node.volume*node.material.density))
    end
    
    # Velocity BCs
    for bc in state.boundaryConditions
        if bc isa BoundaryConditions.VelocityBC
            BoundaryConditions.apply_bc(bc)
        end
    end
    
    Threads.@threads for node in state.nodes
        # Update displacement with average velocity
        node.displacement += node.velocity * state.dt

        # Zero out force
        @atomic node.force = zeros(3)
    end

    # Displacement BCs
    for bc in state.boundaryConditions
        if bc isa BoundaryConditions.DisplacementBC
            BoundaryConditions.apply_bc(bc)
        end
    end

    # Break bonds
    Threads.@threads for bond in state.bonds
        if !bond.isBroken && Bonds.should_break(bond)
            Bonds.break!(bond)
        end
    end

    # Apply bond force to nodes
    Threads.@threads for bond in state.bonds
        Bonds.apply_force(bond)
    end

    # Calculate final velocity
    Threads.@threads for node in state.nodes
        node.velocity = node.velocity*damping + state.dt*0.5 * (node.force / (node.volume * node.material.density))
    end

end

function adr(state, stale::Bool)
    dt = 9999.0

    # Zero force
    Threads.@threads for node in state.nodes
        @atomic node.force = zeros(3)
    end
    
    # Apply bond force to nodes
    Threads.@threads for bond in state.bonds
        Bonds.apply_force(bond)
    end

    apply_contact_force(state)


    # Calculates factor
    cn = 0.0


    cn1s = zeros(Threads.nthreads())
    cn2s = zeros(Threads.nthreads())
    Threads.@threads for node in state.nodes
        # Need old force and old average velocty
        if node.oldAverageVelocity[1] != 0
            cn1s[Threads.threadid()] -= (node.displacement[1]^2 * (node.force[1] - node.oldForce[1])) / (dt * node.oldAverageVelocity[1] * node.stableMass)
        end
        if node.oldAverageVelocity[2] != 0
            cn1s[Threads.threadid()] -= (node.displacement[2]^2 * (node.force[2] - node.oldForce[2])) / (dt * node.oldAverageVelocity[2] * node.stableMass)
        end
        if node.oldAverageVelocity[3] != 0
            cn1s[Threads.threadid()] -= (node.displacement[3]^2 * (node.force[3] - node.oldForce[3])) / (dt * node.oldAverageVelocity[3] * node.stableMass)
        end
        # cn1 += sum(node.displacement .^ 2 .* (node.force - node.oldForce) / Nodes.mass(node))

        cn2s[Threads.threadid()] += sum(node.displacement .^ 2)
    end

    cn1 = sum(cn1s)
    cn2 = sum(cn2s)


    if cn2 != 0.0
        if cn1 / cn2 > 0.0
            cn = 2.0 * sqrt(cn1/cn2)
        end
    end

    if cn > 2.0
        cn = 1.9
    end
    
    Threads.@threads for node in state.nodes
        # if stale
        #     velhalf = dt * node.force * 0.5 / Nodes.mass(node) 
        # else
        #     velhalf = (2.0 - cn * dt) * node.oldAverageVelocity + 2.0 * dt / Nodes.mass(node)
        # end
        averageVelocity::Vector{Float64} = [0.0, 0.0, 0.0]
        if stale
            averageVelocity = 0.5 * dt * node.force / node.stableMass
        else
            averageVelocity = ((2.0 - cn * dt) * node.oldAverageVelocity + 2.0 * (dt / node.stableMass) * node.force) / (2.0 + cn * dt)
        end
        node.velocity = 0.5 * (node.oldAverageVelocity + averageVelocity)
        node.displacement .+= averageVelocity * dt

        node.oldAverageVelocity = copy(averageVelocity)
        node.oldForce = copy(node.force)
    end

    
    # Displacement BCs
    for bc in state.boundaryConditions
        if bc isa BoundaryConditions.DisplacementBC
            BoundaryConditions.apply_bc(bc)
        elseif bc isa BoundaryConditions.StagedLoadingBC
            Threads.@threads for node in bc.nodes
                node.displacement = bc.currentDisplacement
                node.velocity = zeros(3)
                node.oldAverageVelocity = zeros(3)
                @atomic node.force = zeros(3)
            end
        end
    end
end

function relax(state, kethreshold)
    adr(state, true)
    for i in 1:3
        adr(state ,false)
    end    
    kinetic_energy = MoaUtil.KineticEnergy(state.nodes)

    count = 1
    while kinetic_energy > kethreshold
        adr(state, false)
        kinetic_energy = MoaUtil.KineticEnergy(state.nodes)
        count += 1
        if count % 100 == 0
            print("\r", count, " : ", kinetic_energy)
        end
    end
end

function stagedloading(state, kethreshold::Float64)

    # Advance tabs
    for bc in state.boundaryConditions
        if bc isa BoundaryConditions.StagedLoadingBC
            BoundaryConditions.apply_bc(bc)
        end
    end

    # Relax system
    TimeIntegration.relax(state, kethreshold)


    while true
        # Break bonds
        anybroken = fill(false, Threads.nthreads())
        Threads.@threads for bond in state.bonds
            if !bond.isBroken && Bonds.should_break(bond)
                Bonds.break!(bond)
                anybroken[Threads.threadid()] = true
            end
        end

        # println("Have any bonds broken: ", any(anybroken))

        # Relax system
        TimeIntegration.relax(state, kethreshold*10)


        # repeat until no bonds break
        any(anybroken) || break
    end

end

end