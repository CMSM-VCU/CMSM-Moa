module ForceProbes

using ..Nodes
using ..Bonds
using LinearAlgebra

abstract type AbstractForceProbe end


struct ForceProbePlane <: AbstractForceProbe
    bondsPositive::Vector{Bonds.AbstractBond}
    bondsNegative::Vector{Bonds.AbstractBond}
    pointOnPlane::Vector{Float64}
    normal::Vector{Float64}
end

function signed_distance_from_plane(point::Vector{Float64}, pointOnPlane::Vector{Float64}, normal::Vector{Float64})
    return LinearAlgebra.dot(point, normal) - LinearAlgebra.dot(pointOnPlane, normal)
end

function parse_force_probe_plane(inputDict, bonds)
    pointOnPlane::Vector{Float64} = inputDict["pointOnPlane"]
    normal::Vector{Float64} = inputDict["normal"]

    bondsPositive = Vector{Bonds.AbstractBond}()
    bondsNegative = Vector{Bonds.AbstractBond}()

    for bond in bonds
        fromDist = signed_distance_from_plane(convert(Vector{Float64}, bond.from.position), pointOnPlane, normal)
        toDist = signed_distance_from_plane(convert(Vector{Float64}, bond.to.position), pointOnPlane, normal)
        if fromDist < 0 && toDist > 0
            push!(bondsPositive, bond)
        elseif fromDist > 0 && toDist < 0
            push!(bondsNegative, bond)
        end
    end

    return ForceProbePlane(bondsPositive,bondsNegative,pointOnPlane,normal)
end


end