

# Run as julia [PATH TO THIS FILE] [PATH TO INPUT FILE] [PATH TO BASE FOLDER]
include("../../Moa.jl")

function setup(inputfilepath=ARGS[1])
    state = Moa.parse_input(inputfilepath)
    state.dt = 10^-6

    return state
end

function main_dynamic(state)
    Moa.TimeIntegration.dynamic_integration_no_contact(state, 0.0)
end

function main_dynamic_n(state, n)
    for i in 1:n
        Moa.TimeIntegration.dynamic_integration_no_contact(state, 0.0)
    end
end



println("Running Setup..")
@time state = setup()
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
# 
# 
# println("Running Main() 200 times with "*string(Threads.nthreads())*" threads:")
@time main_dynamic_n(state, 100)
