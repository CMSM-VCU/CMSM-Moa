module Moa
include("AbstractTypes.jl")
include("Materials.jl")
include("Nodes.jl")
include("Bonds.jl")
include("ProximitySearch.jl")
include("BoundaryConditions.jl")
include("ForceProbes.jl")
include("MoaUtil.jl")


mutable struct state
    gridSpacing::Float64
    horizon::Float64
    dt::Float64

    nodes::Vector{Nodes.Node}
    bonds::Vector{AbstractTypes.ABond}
    materials::Dict{Int64, AbstractTypes.AMaterial}

    boundaryConditions::Vector{AbstractTypes.ABoundaryCondition}
    forceProbes::Vector{AbstractTypes.AForceProbe}

    contact::Bool
    contactDistance::Float64
    contactCoefficient::Float64
end

include("TimeIntegration.jl")
include("MoaIO.jl")
end