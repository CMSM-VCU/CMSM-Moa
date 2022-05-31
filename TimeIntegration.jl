module TimeIntegration

using ..Materials, ..Nodes, ..Bonds

function dynamic_integration(nodes, bonds, dt)

    Threads.@threads for node in nodes
        # Average velocity for the time step assuming constant acceleration
        node.velocity = node.velocity + (dt*0.5)*(node.force/(node.volume*node.material.density))

        # Update displacement with average velocity
        node.displacement += node.velocity * dt

        # Zero out force
        @atomic node.force = zeros(3)
    end

    # Break bonds
    Threads.@threads for bond in bonds
        if PD.should_break(bond)
            PD.break!(bond)
        end
    end

    # Apply bond force to nodes
    Threads.@threads for bond in bonds
        PD.apply_force(bond)
    end

    Threads.@threads for node in nodes
        # Calculate final velocity
        node.velocity = node.velocity + dt*0.5 * (node.force / (node.volume * node.material.density))
    end

end

end