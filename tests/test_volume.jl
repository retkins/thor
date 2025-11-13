cd(@__DIR__) 



include("../src/thor.jl")
using Printf, Statistics, LinearAlgebra

const N = 1000 
const M = 1000 

function randrange(N, M, minvalue, maxvalue) 
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

source_pts = hcat(randrange(N, 3, 0, 1))
vol = rand(N) 
current_density = 1e8 * rand(N, 3)
target_pts = hcat(randrange(N, 3, 1.5, 2.5))

function test_single(source_pts, vol, current_density, target_pts, phi)

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
    for i in 1:nphi 
        B_direct = bfield(source_pts, vol, current_density, target_pts, direct)
        B_octree = bfield(source_pts, vol, current_density, target_pts, octree; phi=phi_vals[i])

        delta = B_octree .- B_direct
        err = abs.(delta ./ B_direct)
        maxerr[i] = maximum(err) 
        meanerr[i] = mean(err)

    end

    return phi_vals, maxerr, meanerr

end

# phi_vals, maxerr, meanerr = test_multiple(source_pts, vol, current_density, target_pts, 1e-3, 1e0, 5e-3)

# p = plot(phi_vals, meanerr, label="Mean Error")
# plot!(p, phi_vals, maxerr, label="Max Error")

test_single(source_pts, vol, current_density, target_pts, 5e-3)
