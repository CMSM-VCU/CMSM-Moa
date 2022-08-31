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
