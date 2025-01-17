module Nodes

using StaticArrays
using ..Materials
using ..AbstractTypes
using LinearAlgebra: norm

mutable struct Node{M}
    position::MVector{3, Float64}
    displacement::MVector{3, Float64}
    velocity::MVector{3, Float64}
    oldAverageVelocity::MVector{3, Float64}
    stableMass::Float64

    @atomic force::MVector{3, Float64}
    oldForce::MVector{3, Float64}

    volume::Float64
    material::M

    allowFailure::Bool

    family::Vector{AbstractTypes.ABond}
end

# Some constructors
Node(x::Float64, y::Float64, z::Float64, m::Materials.AMaterial, grid_spacing::Float64) = Node{typeof(m)}(
    [x,y,z],
    [0,0,0],
    [0,0,0],
    [0,0,0],
    0.0,
    [0,0,0],
    [0,0,0],
    grid_spacing^3,
    m,
    true,
    [])

Node(pos::MVector{3, Float64}, m::Materials.AMaterial, grid_spacing::Float64) = Node(
    pos,
    [0,0,0],
    [0,0,0],
    [0,0,0],
    0,
    [0,0,0],
    [0,0,0],
    grid_spacing^3,
    m,
    true,
    [])


"Mass of the node"
function mass(node::Node)
    return node.volume * node.material.density
end

"Stable mass calculation for Adaptive Dynamic Relaxation"
function stableMass(node::Node, dt, horizon, gridspacing)
    if node.material isa Materials.TanhElastic || node.material isa Materials.Bundle
        return 0.25 * dt^2 * 4.0/3.0 * pi * horizon^3 * node.material.a * node.material.b / gridspacing
    end
    return 0.25 * dt^2 * 4.0/3.0 * pi * horizon^3 * node.material.bond_constant / gridspacing
end

"Peridynamic damage as ratio of broken bonds to initial bonds.
Nodes with no initial bonds have a damage of -1"
function damage(node::Node)
    numbonds::Int64 = length(node.family)
    numbrokenbonds:: Int64 = 0
    damage::Float64 = 0

    # Count number of broken bonds
    for bond in node.family
        if bond.isBroken
            numbrokenbonds += 1
        end
    end

    # If no family members, damage is -1
    numbonds == 0 ? damage = -1.0 : damage = numbrokenbonds / numbonds
    return damage
end

"Ratio of broken bonds between materials with the same id to total bonds between 
materials with the same id"
function materialDamage(node::Node)
    bonds_relevant = [bond for bond in node.family if bond.from.material.id ∈ bond.to.material.stronglyConnected]
    numbonds::Int64 = length(bonds_relevant)
    numbrokenbonds::Int64 = 0

    # Count number of broken bonds
    for bond in bonds_relevant
        if bond.isBroken
            numbrokenbonds += 1
        end
    end

    damage::Float64 = 0
    numbonds == 0 ? damage = 0 : damage = numbrokenbonds / numbonds
    return damage
end

"Ratio of broken bonds between materials with the same id to total bonds between 
materials with the same id"
function materialDamage(node::Node{Materials.Bundle})
    bonds_relevant = [bond for bond in node.family if bond.from.material.id ∈ bond.to.material.stronglyConnected]
    numbonds::Int64 = length(bonds_relevant)
    damageACC::Float64 = 0

    # Count number of broken bonds
    for bond in bonds_relevant
        if bond.isBroken
            damageACC += 1.0
        else
            if bond.max_strain > node.material.e_soften
                # Max force / Force at max strain
                a = bond.from.material.a
                b = bond.from.material.b
                c = bond.from.material.c
                d = bond.from.material.d
                e = bond.from.material.e
                e_soften = bond.from.material.e_soften
                damageACC += 
                    a*tanh(b*e_soften) /
                    ((a*tanh(b * e_soften) + c*tanh(d*(bond.max_strain - e_soften)) - e * (e_soften - bond.max_strain)) * 0.1666666666666666666666)
            end
        end

    end

    damage::Float64 = 0
    numbonds == 0 ? damage = 0 : damage = damageACC / numbonds
    return damage
end


"Ratio of broken bonds between materials with different ids to total bonds between 
materials with different ids"
function interfaceDamage(node::Node)
    bonds_relevant = [bond for bond in node.family if bond.from.material.id != bond.to.material.id]
    numbonds::Int64 = length(bonds_relevant)
    numbrokenbonds::Int64 = 0

    # Count number of broken bonds
    for bond in bonds_relevant
        if bond.isBroken
            numbrokenbonds += 1
        end
    end

    damage::Float64 = 0
    numbonds == 0 ? damage = 0 : damage = numbrokenbonds / numbonds
    return damage
end

function reference_distance(a::Node, b::Node)
    return norm(b.position - a.position)
end

function deformed_distance(a::Node, b::Node)
    return norm((b.position + b.displacement) .- (a.position + a.displacement))
end


end