using Pkg
Pkg.add(path="S:/CMSM-Moa")
using Moa
using Test

@testset "Moa.jl" begin
    @test Moa.HelloWorld()
end
