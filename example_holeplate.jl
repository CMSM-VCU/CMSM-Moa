# START
include("Moa.jl")
include("Visualize.jl")
Moa.parse_input("example_holeplate.toml")
dt = 9999.
for node in Moa.nodes
    node.stableMass = Moa.Nodes.stableMass(node, dt, Moa.horizon, Moa.gridspacing)
end


## LOOP
for lpnum in 1:3
    # Increment SL stage:
    for trashvar in 1:1
        for bc in Moa.boundaryConditions
            if bc isa Moa.BoundaryConditions.StagedLoadingBC
                Moa.BoundaryConditions.apply_bc(bc)
            end
        end
    end

    # Relax
    Moa.TimeIntegration.relax(Moa.nodes, Moa.bonds, Moa.boundaryConditions, 1.0e-20)


    # Fail and relax until no more failure
    while true
        anybroken = fill(false, Threads.nthreads())

        # Break bonds
        Threads.@threads for bond in Moa.bonds
            if !bond.isBroken && Moa.Bonds.should_break(bond)
                Moa.Bonds.break!(bond)
                anybroken[Threads.threadid()] = true
            end
        end

        println("Have any bonds broken: ", any(anybroken))

        # Relax
        Moa.TimeIntegration.relax(Moa.nodes, Moa.bonds, Moa.boundaryConditions, 1.0e-20)


        # Exit condition
        any(anybroken) || break
    end

    plot(Moa.nodes, 1.0)
end

## PYVISTA VISUALIZE
# plot(Moa.nodes, 1.0)

## WRITE OUTPUT
Moa.write_output("output.csv", Moa.nodes)
