mutable struct Bundle <: AMaterial
    id::Int64

    # Not very relevant because it is 
    density::Float64

    interface_stiffness::Float64
    interface_critical_stretch::Float64

    # This is taken from MD data
    critical_strain::Float64

    # These two are parameters that were fit to the md PolyData
    # stress = a * tanh(b * strain) if strain hasn't passed e_crit
    # Initial stiffness is a * b
    a::Float64
    b::Float64
    c::Float64
    d::Float64
    e::Float64
    e_soften::Float64
    area::Float64

    emod::Float64

    # Strongly Connected To
    stronglyConnected::Vector{Int64}
    
end

Base.copy(material::Bundle) = Bundle(material.id,
    material.density,
    material.interface_stiffness,
    material.interface_critical_stretch,
    material.critical_strain,
    material.a,
    material.b,
    material.c,
    material.d,
    material.e,
    material.e_soften,
    material.area,
    material.emod,
    [])

# TanhElastic(id::Int64,
#             density::Float64,
#             interface_stiffness::Float64,
#             interface_critical_stretch::Float64,
#             critical_strain::Float64,
#             a::Float64,
#             b::Float64,
#             emod::Float64) = TanhElastic(id,
#                                         density,
#                                         interface_stiffness,
#                                         interface_critical_stretch,critical_strain,
#                                         a,
#                                         b,
#                                         emod,
#                                         [])