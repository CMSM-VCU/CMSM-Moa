module TimeIntegration

using ..Materials, ..Nodes, ..Bonds, ..BoundaryConditions, ..AbstractTypes

# function dynamic_stable_timestep()

function dynamic_integration(nodes, bonds, bcs, dt)

    Threads.@threads for node in nodes
        # Average velocity for the time step assuming constant acceleration
        node.velocity = node.velocity + (dt*0.5)*(node.force/(node.volume*node.material.density))
    end
    
    # Velocity BCs
    for bc in bcs
        if bc isa BoundaryConditions.VelocityBoundaryCondition
            BoundaryConditions.apply_bc(bc)
        end
    end
    
    Threads.@threads for node in nodes
        # Update displacement with average velocity
        node.displacement += node.velocity * dt

        # Zero out force
        @atomic node.force = zeros(3)
    end

    # Displacement BCs
    for bc in bcs
        if bc isa BoundaryConditions.DisplacementBoundaryCondition
            BoundaryConditions.apply_bc(bc)
        end
    end

    # Break bonds
    Threads.@threads for bond in bonds
        if Bonds.should_break(bond)
            Bonds.break!(bond)
        end
    end

    # Apply bond force to nodes
    Threads.@threads for bond in bonds
        Bonds.apply_force(bond)
    end

    Threads.@threads for node in nodes
        # Calculate final velocity
        node.velocity = node.velocity + dt*0.5 * (node.force / (node.volume * node.material.density))
    end

end

function adr(nodes::Vector{Nodes.AbstractNode},
                bonds::Vector{AbstractTypes.ABond},
                bcs::Vector{BoundaryConditions.AbstractBoundaryCondition},
                dt::Float64,
                stale::Bool)
    

    Threads.@threads for node in nodes
        @atomic node.force = zeros(3)
    end
    
    # Apply bond force to nodes
    Threads.@threads for bond in bonds
        Bonds.apply_force(bond)
    end


    # Calculates factor
    cn = 0.0
    cn1 = 0.0
    cn2 = 0.0

    for node in nodes
        # Need old force and old average velocty
        if node.oldAverageVelocity[1] != 0
            cn1 -= (node.displacement[1]^2 * (node.force[1] - node.oldForce[1])) / (dt * node.oldAverageVelocity[1] * node.stableMass)
        end
        if node.oldAverageVelocity[2] != 0
            cn1 -= (node.displacement[2]^2 * (node.force[2] - node.oldForce[2])) / (dt * node.oldAverageVelocity[2] * node.stableMass)
        end
        if node.oldAverageVelocity[3] != 0
            cn1 -= (node.displacement[3]^2 * (node.force[3] - node.oldForce[3])) / (dt * node.oldAverageVelocity[3] * node.stableMass)
        end
        # cn1 += sum(node.displacement .^ 2 .* (node.force - node.oldForce) / Nodes.mass(node))

        cn2 += sum(node.displacement .^ 2)
    end

    if cn2 != 0.0
        if cn1 / cn2 > 0.0
            cn = 2.0 * sqrt(cn1/cn2)
        end
    end

    if cn > 2.0
        cn = 1.9
    end
    
    Threads.@threads for node in nodes
        # if stale
        #     velhalf = dt * node.force * 0.5 / Nodes.mass(node) 
        # else
        #     velhalf = (2.0 - cn * dt) * node.oldAverageVelocity + 2.0 * dt / Nodes.mass(node)
        # end
        averageVelocity = zeros(3)
        if stale
            averageVelocity = 0.5 * dt * node.force / node.stableMass
        else
            averageVelocity = ((2.0 - cn * dt) * node.oldAverageVelocity + 2.0 * (dt / node.stableMass) * node.force) / (2.0 + cn * dt)
        end
        node.velocity = 0.5 * (node.oldAverageVelocity + averageVelocity)
        node.displacement += averageVelocity * dt

        # Displacement BCs
        for bc in bcs
            if bc isa BoundaryConditions.DisplacementBoundaryCondition
                BoundaryConditions.apply_bc(bc)
            end
        end
        
        node.oldAverageVelocity = copy(averageVelocity)
        node.oldForce = copy(node.force)
    end
end


end