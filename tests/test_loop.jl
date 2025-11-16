"""
    Test Loop 

Create a current loop around the z-axis. Compare direct sum vs. octree methods. 
Compare both to the theoretical solution along the axis of the loop.
"""

cd(@__DIR__)
include("../src/thor.jl")

using LinearAlgebra, Plots, Statistics

# Select which axis to evaluate over
@enum EvalAxis X=1 Y=2 Z=3

# Test constants 
const side_min = 1e-3       # Minimum length of an element 'side'

function make_loop_problem(N::Integer, M::Integer, R::Real, I::Real, ax::EvalAxis, linemin::Real, linemax::Real)

    # Create a loop around the z-axis 
    # To avoid creating extremely small elements, create a sequence of rings
    # that occupy the same volume in space
    source_pts = zeros(N, 3) 
    vol = ones(N) 
    current_density = zeros(N, 3)

    # Determine how many elements per each ring 
    n_per_ring = floor(R / side_min) 
    nrings = ceil(N / n_per_ring)
    side = side_min
    theta_step = 2 * pi / n_per_ring
    if n_per_ring > N 
        side = R / N 
        if N > 1 
            theta_step = 2 * pi / (N - 1) 
        else 
            theta_step = R/10
        end
    end 

    @printf "N rings : %i\n" nrings

    # Configure each element 
    side = R*theta_step
    area = side*side
    vol .*= area*side
    j = (I / area) / nrings
    theta = 0 
    for i in 1:N 
        source_pts[i,1] = R*cos(theta) 
        source_pts[i,2] = R*sin(theta)
        current_density[i,1] = -j*sin(theta) 
        current_density[i,2] = j*cos(theta) 
        theta += theta_step 
    end

    # Create a line on which to evaluate the field
    target_pts = zeros(M, 3) 
    lineval = linemin
    linestep = (linemax - linemin) / (M - 1)

    for i in 1:M
        target_pts[i,Integer(ax)] = lineval
        lineval += linestep
    end 

    element_radius = (area * side / (4.0 * pi / 3.0)) ^ 1/3.0
    @printf "Element radius: %.6e\n" element_radius

    return source_pts, vol, current_density, target_pts
end

function test_loop_axis(;MAX_ERR::Real=1e-3, phi::Real=1e-3, verbose::Bool=false, make_plot::Bool=false) 

    # Problem parameters
    N = 10000
    M = 1000
    R = 2.0
    I = 3e6
    ax = Z
    linemin = -10*R
    linemax = 10*R

    source_pts, vol, current_density, target_pts = make_loop_problem(N, M, R, I, ax, linemin, linemax)

    B_exact = bfield_loop_axis(I, R, target_pts[:,3])
    B_direct = bfield(source_pts, vol, current_density, target_pts, direct)
    B_octree = bfield(source_pts, vol, current_density, target_pts, octree; phi=phi)

    # Check the field in the z-direction (other components are zero)
    direct_err = abs.(B_direct[:,3] .- B_exact) ./ B_exact
    octree_err = abs.(B_octree[:,3] .- B_exact) ./ B_exact
    mean_direct_err = mean(direct_err)
    max_direct_err = maximum(direct_err)
    mean_octree_err = mean(octree_err)
    max_octree_err = maximum(octree_err)

    if verbose
        @printf "Mean direct error: %.3f %%\n" mean_direct_err * 100.0
        @printf "Max direct error:  %.3f %%\n" max_direct_err * 100.0
        @printf "Mean octree error: %.3f %%\n" mean_octree_err * 100.0
        @printf "Max octree error:  %.3f %%\n" max_octree_err * 100.0
    end

    if make_plot
        p = plot() 
        plot!(p, target_pts[:,Integer(ax)], B_exact, label="Exact")
        scatter!(p, target_pts[:,Integer(ax)], B_direct[:,3], label="Direct", markershape=:circle)
        scatter!(p,target_pts[:,Integer(ax)], B_octree[:,3], label="Octree", markershape=:rect)
        display(p)
    end

    if (mean_direct_err < MAX_ERR) && (max_direct_err < MAX_ERR) && (mean_octree_err < MAX_ERR) && (max_octree_err < MAX_ERR)
        return true 
    else 
        return false 
    end

end 


function test_phi_values(phimin=1e-3, phimax=1e-1, phistep=1e-3)

    # Problem parameters
    N = 10000
    M = 100
    R = 1.0
    I = 1e6
    ax = Z
    linemin = -1.0*R
    linemax = 1.0*R

    phivals = [x for x in phimin:phistep:phimax]
    nphi = length(phivals)
    mean_err = zeros(nphi)
    max_err = zeros(nphi)

    source_pts, vol, current_density, target_pts = make_loop_problem(N, M, R, I, ax, linemin, linemax)

    B_exact = bfield_loop_axis(I, R, target_pts[:,3])
    for i in 1:nphi
        B_octree = bfield(source_pts, vol, current_density, target_pts, octree; phi=phivals[i])

        @printf "Bexact @ min z = %.3f T\n" B_exact[1]
        @printf "Bz @ min z     = %.3f T\n" B_octree[1,3]
        # Check the field in the z-direction (other components are zero)
        octree_err = abs.(B_octree[:,3] .- B_exact) ./ B_exact
        mean_err[i] = mean(octree_err)
        max_err[i] = maximum(octree_err)
    end

    p = plot() 
    plot!(p, phivals, 100.0 .* mean_err, label="Mean Error") #, xscale=:log10)
    plot!(p, phivals, 100.0 .* max_err, label="Max Error") #, yscale=:log10)
    xlabel!(p, "Phi")
    ylabel!(p, "Error [%]")
    display(p)

end

test_phi_values()
