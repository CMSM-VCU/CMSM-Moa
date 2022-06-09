include("Moa.jl")

# For generalization should read input file name from command line
Moa.parse_input("example_holeplate.toml")


## Initialize stable mass (need to internalize in future)

dt = 9999.
for node in Moa.nodes
    node.stableMass = Moa.Nodes.stableMass(node, dt, Moa.horizon, Moa.gridspacing)
end


## Assign staged loading tabs

tab_positive = [node for node in Moa.nodes if node.position[1] > 59]
tab_negative = [node for node in Moa.nodes if node.position[1] < -59]

## Increment

for newdisp in 0:0.1:30
    for node in tab_positive
        node.displacement[1] = newdisp
    end

    for node in tab_negative
        node.displacement[1] = -newdisp
    end

    

end