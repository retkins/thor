
using Test 

cd(@__DIR__)
include("test_loop.jl")

@testset "all tests" begin 
    @test test_loop_axis(; MAX_ERR=1e-2, phi=1e-3, verbose=true)
    # @test loop_inplane()
end

