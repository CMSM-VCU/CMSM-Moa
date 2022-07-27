module Bonds

using StaticArrays: SVector
using LinearAlgebra: norm
using ..Nodes
using ..Materials
using ..AbstractTypes: ABond, AMaterial

mutable struct Bond{M, N} <: ABond
    from::Nodes.Node{M}
    to::Nodes.Node{N}
    isBroken::Bool
end

Bond(from::Nodes.Node{<:AMaterial}, to::Nodes.Node{<:AMaterial}) = Bond(from, to, false)

function get_strain(bond::Bond)
    initial_bond_length::Float64 = norm(bond.to.position - bond.from.position)
    deformed_bond_vector::SVector{3, Float64} = bond.to.position + bond.to.displacement - bond.from.position - bond.from.displacement
    deformed_bond_length::Float64 = norm(deformed_bond_vector)
    return (deformed_bond_length - initial_bond_length) / initial_bond_length
end

"Returns the force of the bond with the minimum material properties"
function get_force(bond::Bond)
    if bond.isBroken
        return zeros(3)
    end

    # Duplicated code here from get_strain because it uses intermediate calculations
    initial_bond_length::Float64 = norm(bond.to.position - bond.from.position)
    deformed_bond_vector::SVector{3, Float64} = bond.to.position + bond.to.displacement - bond.from.position - bond.from.displacement
    deformed_bond_length::Float64 = norm(deformed_bond_vector)
    strain::Float64 = (deformed_bond_length - initial_bond_length) / initial_bond_length

    direction::SVector{3,Float64} = deformed_bond_vector ./ deformed_bond_length

    return  direction * (min(bond.from.material.bond_constant, bond.to.material.bond_constant) * strain * bond.to.volume * bond.from.volume)
end

function get_force(bond::Bond{Materials.TanhElastic, Materials.TanhElastic})
    if bond.isBroken
        return zeros(3)
    end

    # Duplicated code here from get_strain because it uses intermediate calculations
    initial_bond_length::Float64 = norm(bond.to.position - bond.from.position)
    deformed_bond_vector::SVector{3, Float64} = bond.to.position + bond.to.displacement - bond.from.position - bond.from.displacement
    deformed_bond_length::Float64 = norm(deformed_bond_vector)
    strain::Float64 = (deformed_bond_length - initial_bond_length) / initial_bond_length

    direction::SVector{3,Float64} = deformed_bond_vector ./ deformed_bond_length

    if bond.from.material.id == bond.to.material.id
        return  direction * (bond.from.material.a*tanh(bond.from.material.b*strain) * bond.to.volume * bond.from.volume)
    else
        return  direction * (bond.from.material.a*bond.from.material.b*strain * bond.to.volume * bond.from.volume * min(bond.from.material.interface_stiffness_coeff, bond.to.material.interface_stiffness_coeff))
    end
end

"Applies the bond's force to its from node (ATOMIC OPERATION, THREAD SAFE)"
function apply_force(bond::Bond)
    @atomic bond.from.force += get_force(bond)
end

"Returns whether or not a bond should break"
function should_break(bond::Bond)
    return (bond.from.allowFailure && bond.to.allowFailure) && get_strain(bond) > min(bond.from.material.critical_strain, bond.to.material.critical_strain)
end

"Breaks the bond"
function break!(bond::Bond)
    # bond.isBroken && println("WARNING! BREAKING BOND THAT IS ALREADY BROKEN!")
    # println("SNAPPPPP: ", bond.from.position[1], ", ", bond.to.position[1])
    bond.isBroken = true
end

"Removes the bond from nodes families"
function delete(bond::Bond)
    deleteat!(bond.from.family, findall(x->x==bond, bond.from.family))
    deleteat!(bond.to.family, findall(x->x==bond, bond.to.family))
end

end