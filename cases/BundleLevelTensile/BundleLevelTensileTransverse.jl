# Run as julia -t [NUM THREADS] [PATH TO THIS FILE] [PATH TO INPUT FILE]
include("../../Moa.jl")
state = Moa.parse_input("./cases/BundleLevelTensile/HAHBTransverse.toml");
state.dt = 1.0

matIDs = [id for (id,material) in state.materials]

# Read SPLITS
file = open("./cases/BundleLevelTensile/record_splits.txt", "r")
lines = split(read(file, String), "\n")
close(file)
for line in lines
    if line != ""
        root = parse(Int64, strip(split(line, ">")[1]))
        left = parse(Int64, strip(split(strip(split(line, ">")[2]), ",")[1], [' ', '(']))
        right = parse(Int64, strip(split(strip(split(line, ">")[2]), ",")[2],[' ', ')']))
        if root ∈ matIDs && left ∈ matIDs && right ∈ matIDs
            push!(state.materials[root].stronglyConnected, left)
            push!(state.materials[root].stronglyConnected, right)
            push!(state.materials[left].stronglyConnected, root)
            push!(state.materials[right].stronglyConnected, root)
        end
    end
end

# Read MERGES
file = open("./cases/BundleLevelTensile/record_merges.txt", "r")
lines = split(read(file, String), "\n")
close(file)
# println(lines)
for line in lines
    if line != ""
        root = parse(Int64,strip(split(line, "<")[1]))
        left = parse(Int64, strip(split(strip(split(line, "<")[2]), ",")[1], [' ', '(']))
        right = parse(Int64, strip(split(strip(split(line, "<")[2]), ",")[2],[' ', ')']))
        if root ∈ matIDs &&  left ∈ matIDs && right ∈ matIDs
            push!(state.materials[root].stronglyConnected, left)
            push!(state.materials[root].stronglyConnected, right)
            push!(state.materials[left].stronglyConnected, root)
            push!(state.materials[right].stronglyConnected, root)
        end
    end
end

for id in matIDs
    push!(state.materials[id].stronglyConnected, id)
end


output_folder = "./cases/BundleLevelTensile/outputTransverse"
if !isdir(output_folder)
    mkdir(output_folder)
end
open(output_folder*"/HAHBTransverse.csv", "w") do file
    write(file, "Step, Half Displacement [nm], Halfway Force [nN]\n")
end;

for i in 1:75
    Moa.TimeIntegration.stagedloading(state, 0.5e-13, 15000, 400)

    Moa.write_output(output_folder*"/HAHB.h5", state, i)
    
    force = Moa.ForceProbes.measure_force(state.forceProbes[1])
    halfDisplacement = state.boundaryConditions[2].currentDisplacement[3]
    open(output_folder*"/HAHBTransverse.csv", "a") do file
        write(file, string(i) * ", " * string(halfDisplacement) * ", " * string(force) * "\n")
    end;

end