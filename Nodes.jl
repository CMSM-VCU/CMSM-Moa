module Nodes

using StaticArrays
using ..Materials
using ..AbstractTypes
using LinearAlgebra: norm
# abstract type AbstractNode end

# mutable struct Node <: AbstractNode
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


function mass(node::Node)
    return node.volume * node.material.density
end

function stableMass(node::Node, dt, horizon, gridspacing)
    if node.material isa Materials.LinearElastic
        return 0.25 * dt^2 * 4.0/3.0 * pi * horizon^3 * node.material.bond_constant / gridspacing
    elseif node.material isa Materials.TanhElastic
        return 0.25 * dt^2 * 4.0/3.0 * pi * horizon^3 * node.material.a * node.material.b / gridspacing
    end
end

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

function materialDamage(node::Node)
    bonds_relevant = [bond for bond in node.family if bond.from.material.id == bond.to.material.id]
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


function distance(a::Node, b::Node)
    return norm((a.position[1] - b.position[1])^2 +
                (a.position[2] - b.position[2])^2 +
                (a.position[3] - b.position[3])^2)
end

end