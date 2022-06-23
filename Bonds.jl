module Bonds

using StaticArrays: SVector
using LinearAlgebra: norm
using ..Nodes
using ..AbstractTypes

mutable struct Bond <: AbstractTypes.ABond
    from::Nodes.Node
    to::Nodes.Node
    isBroken::Bool
end

Bond(from::Nodes.Node, to::Nodes.Node) = Bond(from, to, false)

function get_strain(bond::Bond)
    initial_bond_length::Float64 = norm(bond.to.position - bond.from.position)
    deformed_bond_vector::SVector{3, Float64} = bond.to.position + bond.to.displacement - bond.from.position - bond.from.displacement
    deformed_bond_length::Float64 = norm(deformed_bond_vector)
    return (deformed_bond_length - initial_bond_length) / initial_bond_length
end

"Returns the force of the bond with the minimum material properties"
function get_force(bond::Bond)
    if bond.isBroken
        return zeros(3,)
    end

    # Duplicated code here from get_strain because it uses intermediate calculations
    initial_bond_length::Float64 = norm(bond.to.position - bond.from.position)
    deformed_bond_vector::SVector{3, Float64} = bond.to.position + bond.to.displacement - bond.from.position - bond.from.displacement
    deformed_bond_length::Float64 = norm(deformed_bond_vector)
    strain::Float64 = (deformed_bond_length - initial_bond_length) / initial_bond_length

    direction::SVector{3,Float64} = deformed_bond_vector ./ deformed_bond_length

    return  direction * (min(bond.from.material.bond_constant, bond.to.material.bond_constant) * strain * bond.to.volume * bond.from.volume)
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
    bond.isBroken && println("WARNING! BREAKING BOND THAT IS ALREADY BROKEN!")
    bond.isBroken = true
end

"Removes the bond from nodes families"
function delete(bond::Bond)
    deleteat!(bond.from.family, findall(x->x==bond, bond.from.family))
    deleteat!(bond.to.family, findall(x->x==bond, bond.to.family))
end

end