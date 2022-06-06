module Materials

abstract type AbstractMaterial end

struct LinearElastic <: AbstractMaterial
    id::Int64
    density::Float64
    critical_strain::Float64
    bond_constant::Float64
end

struct CustomMaterial <: AbstractMaterial
    id::Int64
    density::Float64
end

function parse_material(inputDict)
    # A material needs at least these 3 properties
    @assert haskey(inputDict, "type")
    @assert haskey(inputDict, "density")
    @assert haskey(inputDict, "id")

    if inputDict["type"] == "LinearElastic"
        return Materials.LinearElastic(inputDict["id"], inputDict["density"], inputDict["critical_strain"], inputDict["bond_constant"])
    elseif inputDict["type"] == "Custom"
        return Materials.CustomMaterial(inputDict["id"], inputDict["density"])
    else
        # Material type not known
        println("#### Unknown type from material: ", inputDict["id"])
        throw(Exception)
    end
end


end
