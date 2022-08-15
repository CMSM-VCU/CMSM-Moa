# Run as julia -t [NUM THREADS] [PATH TO THIS FILE] [PATH TO INPUT FILE] [PATH TO BASE FOLDER]
include("../../Moa.jl")
state = Moa.parse_input(ARGS[1])
state.dt = 1.0



using DelimitedFiles
matIDs = [id for (id,material) in state.materials]
base_folder = ARGS[2]


# Read SPLITS
file = open(base_folder*"record_splits.txt", "r")
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

# Anialiates
file = open(base_folder*"record_annihilates.txt", "r")
ids = readdlm(file, ',', Int64)
close(file)
for row in eachrow(ids)
    if row[1] ∈ matIDs && row[2] ∈ matIDs
        push!(state.materials[row[1]].stronglyConnected, row[2])
        push!(state.materials[row[2]].stronglyConnected, row[1])
    end
end


# Dual spawns
file = open(base_folder*"record_dual_spawns.txt", "r")
ids = readdlm(file, ',', Int64)
close(file)
for row in eachrow(ids)
    if row[1] ∈ matIDs && row[2] ∈ matIDs
        push!(state.materials[row[1]].stronglyConnected, row[2])
        push!(state.materials[row[2]].stronglyConnected, row[1])
    end
end


# Read MERGES
file = open(base_folder*"record_merges.txt", "r")
lines = split(read(file, String), "\n")
close(file)
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

# Add themselves to strongly connected
for id in matIDs
    push!(state.materials[id].stronglyConnected, id)
end

output_folder = base_folder*"output_transverse"
if !isdir(output_folder)
    mkdir(output_folder)
end
open(output_folder*"/ForcePlane.csv", "w") do file
    write(file, "Step, Half Displacement [nm], Halfway Force [nN]\n")
end;

for i in 1:75
    Moa.TimeIntegration.stagedloading(state, 10.0e-13, 15000, 400)

    Moa.write_output(output_folder*"/NodeSnapshots.h5", state, i)
    
    force = Moa.ForceProbes.measure_force(state.forceProbes[1])
    halfDisplacement = state.boundaryConditions[2].currentDisplacement[3]
    open(output_folder*"/ForcePlane.csv", "a") do file
        write(file, string(i) * ", " * string(halfDisplacement) * ", " * string(force) * "\n")
    end;

end