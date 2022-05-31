include("Moa.jl")


Moa.parse_input("example_input.toml")

for i in 1:100
    Moa.TimeIntegration.dynamic_integration(Moa.nodes,Moa.bonds, 0.001)
    println(i)
end