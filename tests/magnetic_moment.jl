"""
    test using the magnetic moment and dipole expansion method
"""

using LinearAlgebra, Printf

const R = 100.0 
const I = 1e6 
const mu0 = 4*pi*1e-7 
const mu0_over_4pi = mu0/(4*pi)

const N = 1000 

function make_sources(N, R, I) 

    source_pts = zeros(N, 3) 
    vol = ones(N) 
    current_density = zeros(N, 3) 

    r = R/10 
    area = 4*r*r 
    theta = 0 
    theta_step = 2*pi/(N-1)
    step = theta_step * R 
    v = area*step
    vol .*= v 
    j = I/area

    for i in 1:N 
        source_pts[i,1] = R*cos(theta) 
        source_pts[i,2] = R*sin(theta) 
        current_density[i,1] = -j*sin(theta)
        current_density[i,2] = j*cos(theta)
        theta += theta_step 
    end

    return source_pts, vol, current_density

end

function bfield_direct(source_pts, vol, current_density, target_pts) 

    n = size(source_pts)[1]
    m = size(target_pts)[1]
    B = zeros(m,3)

    for i in 1:n 
        for j in 1:m 
            rp = target_pts[j,:] - source_pts[i,:]

            B[j,:] .+= mu0_over_4pi * vol[i] * cross(current_density[i,:], rp) / norm(rp)^3

        end
    end

    return B 

end 

function bfield_loop_axis(R, I, z) 

    return mu0 * I * R^2 / (2*(z^2 + R^2)^1.5)

end 

function magnetic_moment(source_pts, vol, current_density) 

    n = size(source_pts)[1]
    M = zeros(3) 

    for i in 1:n 
        M += 0.5*vol[i] .* cross(source_pts[i,:], current_density[i,:])
    end
    return M
end 

function bfield_dipole(m, target_pts) 

    n = size(target_pts)[1] 
    B = zeros(n,3) 
    for i in 1:n 
        rhat = target_pts[i,:] ./ norm(target_pts[i,:])
        B[i,:] = mu0_over_4pi * (3*rhat .* dot(m, rhat) .- m) ./ norm(target_pts[i,:])^3
    end

    return B
end

z = 1.0
source_pts, vol, current_density = make_sources(1000, R, I) 
target_pts = zeros(3,3)
target_pts[:,3] .= z
B = bfield_direct(source_pts, vol, current_density, target_pts)
@printf "B direct: %.6f\n" B[1,3]

m = magnetic_moment(source_pts, vol, current_density) 

Baxis = bfield_loop_axis(R, I, z)
@printf "B exact:  %.6f\n" Baxis

Bdipole = bfield_dipole(m, target_pts)
@printf "B dipole: %.6f\n" Bdipole[1,3]

error = (Bdipole[1,3] - Baxis) / Baxis
@printf "Error: %.6f %% \n" error*100.0
