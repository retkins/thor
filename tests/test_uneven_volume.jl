""" Test: source points randomly distributed over a volume; distribution is 
    uneven, leading to unbalanced octree.
"""


cd(@__DIR__) 

include("../src/thor.jl")
using Printf, Statistics, LinearAlgebra

const N = 1000
const M = 1000

"""
    randvec(N::Integer, M::Integer, minvalue::AbstractFloat, maxvalue::AbstractFloat)
    
    Return a vector of random floating point numbers, each in the range 
        `[minvalue, maxvalue]` (inclusive)
"""
function randvec(N::Integer, minvalue::Real, maxvalue::Real) 
    values = zeros(N) 
    for i in 1:N
        values[i] = minvalue + (maxvalue - minvalue) * rand()
    end

    return values 
end 

"""
    randmat(N::Integer, M::Integer, minvalue::AbstractFloat, maxvalue::AbstractFloat)
    
    Return a matrix of random floating point numbers, each in the range 
        `[minvalue, maxvalue]` (inclusive)
"""
function randmat(N::Integer, M::Integer, minvalue::Real, maxvalue::Real) 
    values = zeros(N, M) 
    for i in 1:N*M 
        values[i] = minvalue + (maxvalue - minvalue) * rand()
    end

    return values 
end 

function lorentz(J, B, vol) 

    n = size(J)[1] 
    F = zeros(n, 3) 

    for i in 1:n 
        F[i,:] = cross(J[i,:], B[i,:]) .* vol[i]
    end 
    return F
end



"""
    test_single(source_pts, vol, current_density, target_pts, phi)

"""
function test_single(
    source_pts::Matrix{<:Real}, vol::Vector{<:Real}, 
    current_density::Matrix{<:Real}, target_pts::Matrix{<:Real}, phi::Real
    )

    B_direct = bfield(source_pts, vol, current_density, target_pts, direct)
    B_octree = bfield(source_pts, vol, current_density, target_pts, octree; phi=phi)

    F_direct = lorentz(current_density, B_direct, vol) 
    F_octree = lorentz(current_density, B_octree, vol) 
    Fsum_direct = sum(F_direct, dims=1)
    Fsum_octree = sum(F_octree, dims=1)
    @printf "Direct F sum: \n"
    display(Fsum_direct) 
    @printf "Octree F sum: \n"
    display(Fsum_octree) 

    Ferr = norm((Fsum_octree .- Fsum_direct) ./ Fsum_direct) * 1e2
    @printf "F error: %.2f %%\n" Ferr


    delta = B_octree .- B_direct
    err = abs.(delta ./ B_direct)

    @printf "Max diff: %.6f\n" maximum(delta) 
    @printf "Min diff: %.6f\n" minimum(delta) 
    @printf "Mean diff: %.6f\n" mean(delta)

    @printf "Max error: %.6f\n" maximum(err) 
    @printf "Min error: %.6f\n" minimum(err) 
    @printf "Mean error: %.6f\n" mean(err)
end

function test_multiple(source_pts, vol, current_density, target_pts, phimin, phimax, phistep)

    phi_vals = [x for x in phimin:phistep:phimax]
    nphi = length(phi_vals) 
    maxerr = zeros(nphi) 
    meanerr = zeros(nphi) 
    B_direct = bfield(source_pts, vol, current_density, target_pts, direct)
    for i in 1:nphi 

        B_octree = bfield(source_pts, vol, current_density, target_pts, octree; phi=phi_vals[i])

        delta = B_octree .- B_direct
        err = abs.(delta ./ B_direct)
        maxerr[i] = maximum(err) 
        meanerr[i] = mean(err)

    end

    return phi_vals, maxerr, meanerr

end

source_pts = randmat(N, 3, 0.0, 1.0)
vol = randvec(N, 1e-3, 2e-3)[:,1]
current_density = randmat(N, 3, -1e8, 1e8)
target_pts = randmat(M, 3, 0.0, 1.0)
# target_pts = source_pts

phi_vals, maxerr, meanerr = test_multiple(source_pts, vol, current_density, target_pts, 1e-3, 1e0, 1e-2)
display(meanerr)

p = plot()
plot!(p, phi_vals, 100.0 .* meanerr, label="Mean Error", xscale=:log10)
# plot!(p, phi_vals, 100.0 .* maxerr, label="Max Error")
xlabel!(p, "phi")
ylabel!(p, "Error [%]")
ylims!(p, 0.0, 100.0)

# test_single(source_pts, vol, current_density, target_pts, 5e-3)
