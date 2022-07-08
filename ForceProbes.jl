module ForceProbes

using ..Bonds
using LinearAlgebra: dot
using ..AbstractTypes: AForceProbe

struct ForceProbePlane <: AForceProbe
    bondsPositive::Vector{Bonds.Bond}
    bondsNegative::Vector{Bonds.Bond}
    pointOnPlane::Vector{Float64}
    normal::Vector{Float64}
end

function signed_distance_from_plane(point::Vector{Float64}, pointOnPlane::Vector{Float64}, normal::Vector{Float64})
    return dot(point, normal) - dot(pointOnPlane, normal)
end

"The sum of all normal to the surface of the plane in the reference configuration"
function measure_force(fp::ForceProbePlane)
    force = 0
    for bond in fp.bondsPositive
        force += dot(Bonds.get_force(bond), fp.normal)
    end
    for bond in fp.bondsNegative
        force -= dot(Bonds.get_force(bond), fp.normal)
    end
    return force * 0.5
end

"Creates force probes from input dictionary
(should generalize this)"
function parse_force_probe_plane(inputDict, bonds)
    pointOnPlane::Vector{Float64} = inputDict["pointOnPlane"]
    normal::Vector{Float64} = inputDict["normal"]

    bondsPositive = Vector{Bonds.Bond}()
    bondsNegative = Vector{Bonds.Bond}()

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