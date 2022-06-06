module TimeIntegration

using ..Materials, ..Nodes, ..Bonds, ..BoundaryConditions

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

end