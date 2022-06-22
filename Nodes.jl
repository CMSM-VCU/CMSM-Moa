module Nodes

using StaticArrays
using ..Materials
using ..AbstractTypes

# abstract type AbstractNode end

# mutable struct Node <: AbstractNode
mutable struct Node
    position::MVector{3, Float64}
    displacement::MVector{3, Float64}
    velocity::MVector{3, Float64}
    oldAverageVelocity::MVector{3, Float64}
    stableMass::Float64

    @atomic force::MVector{3, Float64}
    oldForce::MVector{3, Float64}

    volume::Float64
    material::Materials.AbstractMaterial

    family::Vector{AbstractTypes.ABond}
end

# Some constructors
Node(x::Float64, y::Float64, z::Float64, m::Materials.AbstractMaterial, grid_spacing::Float64) = Node(
    [x,y,z],
    [0,0,0],
    [0,0,0],
    [0,0,0],
    0,
    [0,0,0],
    [0,0,0],
    grid_spacing^3,
    m,
    [])

Node(pos::MVector{3, Float64}, m::Materials.AbstractMaterial, grid_spacing::Float64) = Node(
    pos,
    [0,0,0],
    [0,0,0],
    [0,0,0],
    0,
    [0,0,0],
    [0,0,0],
    grid_spacing^3,
    m,
    [])


function mass(node::Node)
    return node.volume * node.material.density
end

function stableMass(node::Node, dt, horizon, gridspacing)
    0.25 * dt^2 * 4.0/3.0 * pi * horizon^3 * node.material.bond_constant / gridspacing
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

end