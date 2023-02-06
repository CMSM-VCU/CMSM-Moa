module Bonds

using StaticArrays: SVector
using LinearAlgebra: norm
using ..Nodes: Nodes, Node
using ..Materials
using ..AbstractTypes: ABond, AMaterial

mutable struct Bond{M, N} <: ABond
    from::Nodes.Node{M}
    to::Nodes.Node{N}
    bondConstant::Float64
    isBroken::Bool
    max_strain::Float64
end
Bond(from::Nodes.Node{<:AMaterial}, to::Nodes.Node{<:AMaterial}) = Bond(from, to, bondconstant(to, from), false, 0.0)

"This is the number strain gets multiplied by to get force"
function bondconstant(to::Node{Materials.LinearElastic}, from::Node{Materials.LinearElastic})
    return (min(from.material.bond_constant, to.material.bond_constant)) * to.volume * from.volume
end

"This isn't used by the getforce function, but it is the initial slope of the force response"
function bondconstant(to::Node{Materials.TanhElastic}, from::Node{Materials.TanhElastic})
    return min(to.material.a * to.material.b, from.material.a * from.material.b) * to.volume * from.volume
end

"This isn't used by the getforce function, but it is the initial slope of the force response"
function bondconstant(to::Node{Materials.Bundle}, from::Node{Materials.Bundle})
    return min(to.material.a * to.material.b, from.material.a * from.material.b) * to.volume * from.volume
end


"The force of the bond with the minimum material properties between two LinearElastic 
Nodes"
function getforce(bond::Bond{Materials.LinearElastic, Materials.LinearElastic})
    # Duplicated code here from get_strain because it uses intermediate calculations
    initial_bond_length::Float64 = norm(bond.to.position - bond.from.position)
    deformed_bond_vector::SVector{3, Float64} = bond.to.position + bond.to.displacement - bond.from.position - bond.from.displacement
    deformed_bond_length::Float64 = norm(deformed_bond_vector)
    strain::Float64 = (deformed_bond_length - initial_bond_length) / initial_bond_length
    
    if bond.isBroken && strain < 0.0
        return zeros(3)
    end

    direction::SVector{3,Float64} = deformed_bond_vector ./ deformed_bond_length

    return  direction * (bond.bondConstant * strain)
end




mutable struct PlasticBond{M, N} <: ABond
    from::Nodes.Node{M}
    to::Nodes.Node{N}
    bondConstant::Float64
    isBroken::Bool
    maxStrain::Float64
end
PlasticBond(from::Nodes.Node{Materials.CustomPlastic}, to::Nodes.Node{Materials.CustomPlastic}) = PlasticBond(from, to, bondconstant(to, from),false, 0.0)

"The strain of a bond as the change in length of the bond relative to its initial 
length"
function getstrain(bond::Bond)
    initial_bond_length::Float64 = norm(bond.to.position - bond.from.position)
    deformed_bond_vector::SVector{3, Float64} = bond.to.position + bond.to.displacement - bond.from.position - bond.from.displacement
    deformed_bond_length::Float64 = norm(deformed_bond_vector)
    return (deformed_bond_length - initial_bond_length) / initial_bond_length
end


"The force enacted by the bond's to node onto the bond's from node for a bond between 
two TanhElastic Nodes"
function getforce(bond::Bond{Materials.TanhElastic, Materials.TanhElastic})
    # Duplicated code here from get_strain because it uses intermediate calculations
    initial_bond_length::Float64 = norm(bond.to.position - bond.from.position)
    deformed_bond_vector::SVector{3, Float64} = bond.to.position + bond.to.displacement - bond.from.position - bond.from.displacement
    deformed_bond_length::Float64 = norm(deformed_bond_vector)
    strain::Float64 = (deformed_bond_length - initial_bond_length) / initial_bond_length

    if strain > 0.0 && bond.isBroken
        return zeros(3)
    end
    direction::SVector{3,Float64} = deformed_bond_vector ./ deformed_bond_length

    # NON-GENERALIZED IMPLEMENTATION: CHECKS IF IT IS IN THE STRONGLY CONNECTED LIST
    if bond.from.material.id ∈ bond.to.material.stronglyConnected
        return  direction *
        (
            bond.from.material.a*tanh(bond.from.material.b * strain) *
            bond.to.volume *
            bond.from.volume
        )
    else
        bond.isBroken && return zeros(3)
        return  direction * 
        (
            min(bond.from.material.a*bond.from.material.b*strain,bond.to.material.a*bond.to.material.b) *
            bond.to.volume *
            bond.from.volume *
            min(bond.from.material.interface_stiffness, bond.to.material.interface_stiffness)
        )
    end
end


"The force enacted by the bond's to node onto the bond's from node for a bond between 
two TanhElastic Nodes"
function getforce(bond::Bond{Materials.Bundle, Materials.Bundle})
    # Duplicated code here from get_strain because it uses intermediate calculations
    initial_bond_length::Float64 = norm(bond.to.position - bond.from.position)
    deformed_bond_vector::SVector{3, Float64} = bond.to.position + bond.to.displacement - bond.from.position - bond.from.displacement
    deformed_bond_length::Float64 = norm(deformed_bond_vector)
    strain::Float64 = (deformed_bond_length - initial_bond_length) / initial_bond_length
    if strain > bond.max_strain
        bond.max_strain = strain
    end
    if strain > 0.0 && bond.isBroken
        return zeros(3)
    end
    direction::SVector{3,Float64} = deformed_bond_vector ./ deformed_bond_length

    
    if bond.from.material.id ∈ bond.to.material.stronglyConnected        # If not an interface bond...
        a = bond.from.material.a
        b = bond.from.material.b
        c = bond.from.material.c
        d = bond.from.material.d
        e = bond.from.material.e
        e_soften = bond.from.material.e_soften
        if bond.max_strain > bond.from.material.e_soften && strain < bond.max_strain        # If the bond is beyond the softening strain but less than the max strain...
            # Lineary interpolate down
            magnitude =  max(0,
                (((a*tanh(b * e_soften) + c*tanh(d*(bond.max_strain - e_soften)) - e * (e_soften - strain)) * 0.1666666666666666666666) +
                ((strain - bond.max_strain) * a * b)) * bond.from.volume * bond.to.volume)
            return direction * magnitude

        end
        if strain > bond.from.material.e_soften         # If the bond is softened, but at it's peak strain
            return  direction *
            (
                (a*tanh(b * e_soften) + c*tanh(d*(strain - e_soften)) - e * (e_soften - strain)) *
                0.1666666666666666666666 *
                bond.to.volume *
                bond.from.volume
            )
        else        # If the bond hasn't softened (Elastic, follows tanh function)
            return  direction *
            (
                bond.from.material.a*tanh(bond.from.material.b * strain) *
                bond.to.volume *
                bond.from.volume * 0.1666666666666666666666
            )
        end
    else
        bond.isBroken && return zeros(3)
        return  direction * 
        (
            min(bond.from.material.a*bond.from.material.b*strain, bond.to.material.a*bond.to.material.b) *
            bond.to.volume *
            bond.from.volume *
            min(bond.from.material.interface_stiffness, bond.to.material.interface_stiffness)
        )
    end
end


function custombondconstant(strain::Float64, material::Materials.CustomPlastic)
    for i in 1:size(material.customForceResponse)[1]-1
        if strain < material.customForceResponse[i+1, 1]
            # Linearlly interpolate using this strain and the next strain
            return material.customForceResponse[i, 2] +
                (strain - material.customForceResponse[i, 1]) *
                (material.customForceResponse[i+1, 2] - material.customForceResponse[i, 2]) /
                (material.customForceResponse[i+1, 1] - material.customForceResponse[i, 1])
        end
    end

    # If reached here, linearly extrapolate to this point
    i = size(material.customForceResponse)[1]-1
    return material.customForceResponse[i, 2] +
        (strain - material.customForceResponse[i, 1]) *
        (material.customForceResponse[i+1, 2] - material.customForceResponse[i, 2]) /
        (material.customForceResponse[i+1, 1] - material.customForceResponse[i, 1])
end

"The force of a bond bewtween plastic materials"
function getforce(bond::PlasticBond{Materials.CustomPlastic, Materials.CustomPlastic})
    # Find strain
    # Duplicated code here from get_strain because it uses intermediate calculations
    initial_bond_length::Float64 = norm(bond.to.position - bond.from.position)
    deformed_bond_vector::SVector{3, Float64} = bond.to.position + bond.to.displacement - bond.from.position - bond.from.displacement
    deformed_bond_length::Float64 = norm(deformed_bond_vector)
    strain::Float64 = (deformed_bond_length - initial_bond_length) / initial_bond_length

    force = zeros(3)
    
    if bond.isBroken && strain > 0.0
        return force
    end
    
    direction::SVector{3,Float64} = deformed_bond_vector ./ deformed_bond_length
    
    # Evaluate maxStrain
    if strain > bond.maxStrain
        bond.maxStrain = strain
        # Use force from custom force response
        return direction * (max(0.0, custombondconstant(strain, bond.from.material)) * bond.to.volume * bond.from.volume)
    else
        scalar_force = 0.0
        if strain > 0
            scalar_force = max(0.0, custombondconstant(bond.maxStrain, bond.from.material) + (strain - bond.maxStrain)*bond.from.material.bond_constant)
        else
            scalar_force = strain*bond.from.material.bond_constant
        end
        # Find the force at the max strain and interpolate down using the material's initial bond constant and clamp at 0
        return direction*(scalar_force * bond.to.volume * bond.from.volume)
    end


end

"Applies the bond's force to its from node (ATOMIC OPERATION, THREAD SAFE)"
function applyforce!(bond::Bond)
    @atomic bond.from.force += getforce(bond)
    @atomic bond.to.force -= getforce(bond)
end

"Returns whether or not a bond should break"
function shouldbreak(bond::Bond)
    return bond.from.allowFailure &&
            bond.to.allowFailure &&
            getstrain(bond) > min(bond.from.material.critical_strain, bond.to.material.critical_strain)
end


"Returns whether or not a bond should break"
function shouldbreak(bond::Bond{Materials.TanhElastic, Materials.TanhElastic})
    if bond.from.material.id ∈ bond.to.material.stronglyConnected
        return bond.from.allowFailure &&
                bond.to.allowFailure &&
                getstrain(bond) > min(bond.from.material.critical_strain, bond.to.material.critical_strain)
    else
        return bond.from.allowFailure &&
                bond.to.allowFailure &&
                getstrain(bond) > min(bond.from.material.interface_critical_stretch, bond.to.material.interface_critical_stretch)
    end
end

"Returns whether or not a bond should break"
function shouldbreak(bond::Bond{Materials.Bundle, Materials.Bundle})
    if bond.from.material.id ∈ bond.to.material.stronglyConnected
        return bond.from.allowFailure &&
                bond.to.allowFailure &&
                getstrain(bond) > min(bond.from.material.critical_strain, bond.to.material.critical_strain)
    else
        return bond.from.allowFailure &&
                bond.to.allowFailure &&
                getstrain(bond) > min(bond.from.material.interface_critical_stretch, bond.to.material.interface_critical_stretch)
    end
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