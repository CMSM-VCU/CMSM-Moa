

# Run as julia [PATH TO THIS FILE] [PATH TO INPUT FILE] [PATH TO BASE FOLDER]
using DelimitedFiles
include("../../Moa.jl")
base_folder = ARGS[2]
output_folder = base_folder*"output_XSEDE"

function setup(inputfilepath=ARGS[1], base_folder=ARGS[2], output_folder=ARGS[2]*"output_XSEDE")
    state = Moa.parse_input(inputfilepath)
    state.dt = 1.0


    matIDs = [id for (id,material) in state.materials]


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

    # if !isdir(output_folder)
    #     mkdir(output_folder)
    # end
    # open(output_folder*"/ForcePlane.csv", "w") do file
    #     write(file, "Step, Half Displacement [nm], Halfway Force [nN]\n")
    # end;

    return state
end

function main_stagedloading(state)
    Moa.TimeIntegration.stagedloading(state, 14.0e-13, 102, 100, false)
end

function main_dynamic(state)
    Moa.TimeIntegration.dynamic_integration_no_contact(state, 0.0)
end

function main_dynamic_n(state, n)
    for i in 1:n
        Moa.TimeIntegration.dynamic_integration_no_contact(state, 0.0)
    end
end

function output(state, i, output_folder=ARGS[2]*"output_XSEDE")
    Moa.write_output(output_folder*"/NodeSnapshots.h5", state, i)
        
    force = Moa.ForceProbes.measure_force(state.forceProbes[1])
    halfDisplacement = state.boundaryConditions[2].currentDisplacement[3]
    open(output_folder*"/ForcePlane.csv", "a") do file
        write(file, string(i) * ", " * string(halfDisplacement) * ", " * string(force) * "\n")
    end;

end

println("Running Setup..")
state = setup()
println("Running Main() with "*string(Threads.nthreads())*" threads:")
@time main_dynamic(state)

@time main_dynamic(state)
@time main_dynamic(state)
@time main_dynamic(state)
@time main_dynamic(state)
@time main_dynamic(state)
@time main_dynamic(state)
@time main_dynamic(state)
@time main_dynamic(state)
@time main_dynamic(state)
@time main_dynamic(state)


println("Running Main() 200 times with "*string(Threads.nthreads())*" threads:")
@time main_dynamic_n(state, 200)
