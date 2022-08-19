module Materials

using ..AbstractTypes: AMaterial
using DelimitedFiles


mutable struct LinearElastic <: AMaterial
    id::Int64
    density::Float64
    critical_strain::Float64
    bond_constant::Float64
    emod::Float64
end

Base.copy(material::LinearElastic) = LinearElastic(
        material.id,
        material.density,
        material.critical_strain,
        material.bond_constant,
        material.emod)


mutable struct TanhElastic <: AMaterial
    id::Int64

    # Not very relevant because it is 
    density::Float64

    interface_stiffness_coeff::Float64

    # This is taken from MD data
    critical_strain::Float64

    # These two are parameters that were fit to the md PolyData
    # stress = a * tanh(b * strain)
    # Initial stiffness is a * b
    a::Float64
    b::Float64
    emod::Float64

    # Strongly Connected To
    stronglyConnected::Vector{TanhElastic}
    
end

Base.copy(material::TanhElastic) = TanhElastic(material.id,
    material.density,
    material.interface_stiffness_coeff,
    material.critical_strain,
    material.a,
    material.b,
    material.emod,
    [])

TanhElastic(id::Int64,density::Float64,interface_stiffness_coeff::Float64,critical_strain::Float64,a::Float64,b::Float64,emod::Float64)= TanhElastic(id,density,interface_stiffness_coeff,critical_strain,a,b,emod, [])

struct CustomPlastic <: AMaterial
    id::Int64
    density::Float64
    critical_strain::Float64

    # initial bond constant (only needed for stable timestep estimation, if yields a lot timestep will be severe underestimation)
    bond_constant::Float64
    emod::Float64
    # Stress-strain response first column is strain, second column is force density [N/m^6]
    customForceResponse::Matrix{Float64}

    # # Overrides default constructor to sort stress-strain response
    # CustomPlastic(id::Int64,
    #         density::Float64,
    #         critical_strain::Float64,
    #         bond_constant::Float64,
    #         emod::Float64,
    #         customForceResponse::Matrix{Float64}) = 
    #         new(id,
    #         density,
    #         critical_strain,
    #         bond_constant,
    #         emod,
    #         customForceResponse[sortperm(customForceResponse[:,1]), :])
end



Base.copy(material::CustomPlastic) = CustomPlastic(
        material.id,
        material.density,
        material.critical_strain,
        material.bond_constant,
        material.emod,
        material.customForceResponse)



"Creates material from input dictionary"
function parse_material(input)
    # A material needs at least these 3 properties
    @assert haskey(input, "type")
    @assert haskey(input, "density")
    @assert haskey(input, "id")

    if input["type"] == "LinearElastic"
        return Materials.LinearElastic(
                                        input["id"],
                                        input["density"],
                                        input["critical_strain"],
                                        input["bond_constant"],
                                        Float64(input["bond_constant"]*6)
                                    )

    elseif input["type"] == "TanhElastic"
        @assert haskey(input, "interface_stiffness_coeff")
        return Materials.TanhElastic(
                                        input["id"],
                                        input["density"],
                                        input["interface_stiffness_coeff"],
                                        input["critical_strain"],
                                        input["a"],
                                        input["b"],
                                        Float64(input["a"]) * Float64(input["b"])
                                    )
    elseif input["type"] == "CustomPlastic"
        @assert haskey(input, "path")

        forceResponse = readdlm(input["path"], ',', Float64, skipstart=1)

        return CustomPlastic(input["id"],input["density"],input["critical_strain"], input["bond_constant"], Float64(input["bond_constant"])*6, forceResponse)

        # Parse path
        # Create struct
        
    else
        # Material type not known
        println("Unknown type from material: ", input["id"])
        throw(Exception)
    end
end




end