""" Test the multipole expansion method 

Bmono = mu0/4pi * vol * J x r / |r|^3 
Bdip = mu0/4pi * (3*r * (m*r) - r^2 * m) / r^5
where: m = 0.5 * sum( (re - rc) x (Je*Ve))
re -> position of element 
rc -> position of cluster centroid 

Tests show that the multipole expansion might work in Barnes Hut, though unclear
if the tests generalize to actual element distributions.
"""

using Statistics, LinearAlgebra, Printf, Plots

const MU0_OVER_4PI = 1e-7 
const MU0 = MU0_OVER_4PI * 4*pi

""" Compute the magnetic moment about the centroid of a cluster of points 
"""
function magnetic_moment(source_pts, vol, current_density) 
    rc = mean(source_pts, dims=1)[:]       # centroid of cluster, as a vector
    m = zeros(3) 
    for i in 1:size(source_pts)[1] 
        re = source_pts[i,:]
        rp = re .- rc 
        Je = current_density[i,:]
        m .+= cross(rp, vol[i] .* Je)
    end
    return 0.5 .* m
end

function magnetic_monopole(vol::Vector{<:Real}, current_density) 
    P = zeros(3) 
    for i in 1:length(vol) 
        P .+= vol[i] * current_density[i,:]
    end
    return P 
end

function bfield_direct(source_pts, vol, current_density, target_pts) 
    B = zeros(size(target_pts))
    for i in 1:size(source_pts)[1]
        Je = current_density[i,:] 
        re = source_pts[i,:]
        for j in 1:size(target_pts)[1] 
            rp = target_pts[j,:] .- re
            B[j,:] .+= MU0_OVER_4PI * vol[i] * cross(Je, rp) / norm(rp)^3
        end
    end
    return B
end

function bfield_monopole(source_pts, vol, current_density, target_pts) 
    p = magnetic_monopole(vol, current_density)
    rc = mean(source_pts, dims=1)[:]
    B = zeros(size(target_pts))
    for j in 1:size(target_pts)[1] 
        rp = target_pts[j,:] - rc
        n = rp ./ norm(rp)
        B[j,:] += MU0_OVER_4PI * cross(p, rp) / norm(rp)^3
        # B[j,:] += MU0_OVER_4PI * (3 .*n .*(dot(m,n)) .- m) / norm(rp)^3 
        # B[j,:] += MU0_OVER_4PI * (3 * rp .* dot(m, rp) .- dot(rp,rp) .* m) ./ norm(rp)^5
    end
    return B
end

function bfield_multipole(source_pts, vol, current_density, target_pts) 
    m = magnetic_moment(source_pts, vol, current_density)
    p = magnetic_monopole(vol, current_density)
    rc = mean(source_pts, dims=1)[:]
    B = zeros(size(target_pts))
    for j in 1:size(target_pts)[1] 
        rp = target_pts[j,:] - rc
        n = rp ./ norm(rp)
        B[j,:] += MU0_OVER_4PI * cross(p, rp) / norm(rp)^3
        B[j,:] += MU0_OVER_4PI * (3 .*n .*(dot(m,n)) .- m) / norm(rp)^3 
        # B[j,:] += MU0_OVER_4PI * (3 * rp .* dot(m, rp) .- dot(rp,rp) .* m) ./ norm(rp)^5
    end
    return B
end

# create a box of source points around (0,0,0)
# box has side length 1.0
n = 1000
source_pts = rand(n,3) .- 0.5
vol = rand(n) .+ 0.01
current_density = (rand(n,3) .-0.5) .* 1e8
# current_density[:,2] = 1e8 * rand(n) 
current_density[:,1] = 1e8 .* (rand(n) .+ 1.0)
current_density[:,2] = 1e8 .* (rand(n) .+ 1.0)
current_density[:,3] = 1e8 .* (rand(n) .+ 1.0)


# target points at distance 
m = 100
theta = [x for x in LinRange(0.001,0.1, m)]
# theta = s / r; s = 1.0 
target_pts = zeros(m,3)
target_pts[:,1] .= 1 ./theta

# Solution
Bdirect = bfield_direct(source_pts, vol, current_density, target_pts)
@printf "B direct:    (%.6f, %.6f, %.6f) \n" Bdirect[1,1] Bdirect[1,2] Bdirect[1,3]

Bmonopole = bfield_monopole(source_pts, vol, current_density, target_pts)
@printf "B monopole:  (%.6f, %.6f, %.6f) \n" Bmonopole[1,1] Bmonopole[1,2] Bmonopole[1,3]
Bmultipole = bfield_multipole(source_pts, vol, current_density, target_pts) 
@printf "B multipole: (%.6f, %.6f, %.6f) \n" Bmultipole[1,1] Bmultipole[1,2] Bmultipole[1,3]

merrx = (Bmultipole[:,1] .- Bdirect[:,1]) ./ Bdirect[:,1]
merry = (Bmultipole[:,2] .- Bdirect[:,2]) ./ Bdirect[:,2]
merrz = (Bmultipole[:,3] .- Bdirect[:,3]) ./ Bdirect[:,3]
perrx = (Bmonopole[:,1] .- Bdirect[:,1]) ./ Bdirect[:,1]
perry = (Bmonopole[:,2] .- Bdirect[:,2]) ./ Bdirect[:,2]
perrz = (Bmonopole[:,3] .- Bdirect[:,3]) ./ Bdirect[:,3]

mnorm_erry = zeros(m)
mnorm_errz = zeros(m)
mnorm_err = zeros(m)
pnorm_erry = zeros(m)
pnorm_errz = zeros(m)
pnorm_err = zeros(m)
for i in 1:m
    # norm_errx[i] = norm(errx[i])
    mnorm_erry[i] = norm(merry[i])
    mnorm_errz[i] = norm(merrz[i])
    mnorm_err[i] = norm([mnorm_erry[i], mnorm_errz[i]])
    pnorm_erry[i] = norm(perry[i])
    pnorm_errz[i] = norm(perrz[i])
    pnorm_err[i] = norm([pnorm_erry[i], pnorm_errz[i]])
end
@printf "Multipole error at theta = %.3f is \n" theta[1] 
# @printf "\tx: %.3f %%\n" mnorm_errx[1]*100.0
@printf "\ty: %.3f %%\n" mnorm_erry[1]*100.0
@printf "\tz: %.3f %%\n" mnorm_errz[1]*100.0
@printf "Multipole error at theta = %.3f is \n" theta[end] 
# @printf "\tx: %.3f %%\n" mnorm_errx[end]*100.0
@printf "\ty: %.3f %%\n" mnorm_erry[end]*100.0
@printf "\tz: %.3f %%\n" mnorm_errz[end]*100.0
# plot(theta, Bmultipole[:,3], label="multipole")
# plot!(theta, Bdirect[:,3], label="direct")

# p = plot(theta,norm_errx, label="x")
p = plot()
plot!(p, theta,pnorm_err, label="monopole")
plot!(p, theta,mnorm_err, label="multipole")
# ylims!(p, 0, 0.1)
# plot(theta, norm_errz)