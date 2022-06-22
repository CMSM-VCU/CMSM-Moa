module TimeIntegration

using ..Materials, ..Nodes, ..Bonds, ..BoundaryConditions, ..AbstractTypes, ..MoaUtil, ..TimeIntegration

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

function adr(nodes::Vector{Nodes.Node},
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


    cn1s = zeros(Threads.nthreads())
    cn2s = zeros(Threads.nthreads())
    Threads.@threads for node in nodes
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
    
    Threads.@threads for node in nodes
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
    for bc in bcs
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

function relax(nodes, bonds, bcs, kethreshold)
    println("Relaxing system...")

    adr(nodes, bonds,bcs, 9999., true)
    
    kinetic_energy = MoaUtil.KineticEnergy(nodes)
    count = 1

    while kinetic_energy > kethreshold
        TimeIntegration.adr(nodes, bonds, bcs, 9999., false)
        kinetic_energy = MoaUtil.KineticEnergy(nodes)
        count += 1
        print("\r", count, " : ", kinetic_energy)
    end
    println("\nFinished realxation!")
end

function stagedloading(nodes, bonds, bcs, kethreshold::Float64)

    # Advance tabs
    for bc in bcs
        if bc isa BoundaryConditions.StagedLoadingBC
            BoundaryConditions.apply_bc(bc)
        end
    end

    # Relax system
    TimeIntegration.relax(nodes, bonds, bcs, kethreshold)


    while true
        # Break bonds
        anybroken = fill(false, Threads.nthreads())
        Threads.@threads for bond in bonds
            if !bond.isBroken && Bonds.should_break(bond)
                Bonds.break!(bond)
                anybroken[Threads.threadid()] = true
            end
        end

        println("Have any bonds broken: ", any(anybroken))

        # Relax system
        TimeIntegration.relax(nodes, bonds, bcs, kethreshold)


        # repeat until no bonds break
        any(anybroken) || break
    end

end

end