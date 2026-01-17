""" Test the octree code on the solenoid problem 
"""

cd(@__DIR__)
include("../../src/thor.jl")
using Plots, Printf

# Problem parameters
const MU0_OVER_4PI = 1e-7 
const MU0 = 4*pi *MU0_OVER_4PI 
const h = 1e-2                      # [m], element size
const n = 1/h                       # [turns/m] and [elements/turn]
const I = 1000                      # [A] 
const r = 1/(2pi)                   # [m], solenoid radius 
const L = 20*r                      # [m], solenoid length
const Bz_expected = MU0 * n * I     # [T], Bz 

# Setup functions 

# Mesh a solenoid given the problem parameters 
function make_solenoid(N::Int64) 

    if N < 32000 
        N = 32000 
    end

    pts = zeros(N,3) 
    vol = (h^3) .* ones(N) 
    J = zeros(N,3) 
    m = convert(Int64, ceil(n))     # shadows `n` above
    area = h*h 
    j = I/area

    if N % m != 0 
        N = ceil(N/m) * m 
    end

    theta_step = 2*pi/n
    zstep = h 
    z = L/2 

    # Start from z=0 and mesh each turn of the solenoid individually
    k = 1 
    for _ in 2:N/m 
        theta = 0 
        for _ in 1:m 
            pts[k,1] = r*cos(theta) 
            pts[k,2] = r*sin(theta)
            pts[k,3] = z
            J[k,1] = -j*sin(theta) 
            J[k,2] = j*cos(theta)
            theta += theta_step 
            k+=1
        end
        # z += zstep 
    end

    return pts, vol, J
end

function make_target_pts(M::Int64) 
    target_pts = zeros(M,3)
    Lstep = L / (M - 1.0) 
    for i in 2:M 
        target_pts[i,3] = target_pts[i-1,3] + Lstep 
    end 

    return target_pts
end

# Measure the RMS error of y against known x 
# Note: ignores x[i] < 1e-8
function rms_error(y::Vector{<:Real}, x::Vector{<:Real})
    err = 0.0
    if length(y) != length(x) 
        @printf "Error in rms_error(): y and x are not the same length!\n"
        return 1e32
    end
    for i in 1:length(y) 
        if x[i] > 1e-8 
            err += ((y[i] - x[i]) / x[i])^2
        end
    end
    return sqrt(err/length(y))
end

N = 20000 
M = 1000
phi = 1e-1

source_pts, vol, J = make_solenoid(N) 
target_pts = make_target_pts(M) 

Bexact = MU0 * n * I 
# time_octree = @elapsed begin 
    Boctree = bfield(source_pts, vol, J, target_pts, octree; phi=phi)
# end
# time_direct = @elapsed begin 
    Bdirect = bfield(source_pts, vol, J, target_pts, direct) 
# end

# Errors at the midpoint (to compare against analytical)
i = convert(Int64, floor(M/2.0))
rms_error_octree = (Boctree[i,3] - Bexact) / Bexact
rms_error_direct = (Bdirect[i,3] - Bexact) / Bexact
@printf "Comparing octree and direct methods to exact solution: \n---\n"
@printf "Expected value: %.3f T \n" Bexact 
@printf "Octree result:  %.3f T (%.2f %% error)\n" Boctree[i,3] rms_error_octree*100.0
@printf "Direct result:  %.3f T (%.2f %% error)\n" Bdirect[i,3] rms_error_direct*100.0
@printf "\n"
@printf "Comparing octree to direct method: \n---\n"
rms_error_octree = rms_error(Boctree[:,3], Bdirect[:,3])
@printf "Octree MRS error relative to direct: %.3f %% \n" rms_error_octree*100.0

@printf "\n"
@printf "Time for direct solution: %.3f s\n" time_direct 
@printf "Time for octree solution: %.3f s\n" time_octree 
@printf "Octree speedup: %.1fx \n" time_direct/time_octree

p = plot()
plot!(p, target_pts[:,3], Boctree[:,3], label="octree")
plot!(p, target_pts[:,3], Bdirect[:,3], label="direct", linestyle=:dash)
xlabel!(p, "Distance Along Solenoid Axis [m]")
ylabel!(p, "Magnetic Flux Density (Bz) [T]")
title!(p, "Solenoid Test: Octree vs Direct Integration (\$\\phi\$="*string(phi)*")")



display(p)

# p2 = plot() 
# scatter!(source_pts[:,1], source_pts[:,2], source_pts[:,3])

