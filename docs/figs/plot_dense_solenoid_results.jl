cd(@__DIR__)

using Plots 

# phi = 0.1
# Test performed on 12/4/2025
# Test machine: AMD Ryzen 9 9950X, single-threaded execution
# Compared direct SIMD vs. recursive octree
nvals = [99300, 206494, 726600]
ninteractions = nvals.^2

# Loop times (note we had to use phi = 0.005 to keep error down)
l_octree_times = [0.41, 0.844, 2.99]
l_direct_times = [4.23, 18.56, 218]

# Solenoid times
s_octree_times = [3.8, 8.71, 36.6]
s_direct_times = [4.17, 18.51, 231.1]
# Dense Solenoid times 
ds_octree_times = [2.079, 4.04, 16.025]
ds_direct_times = [3.96, 17.75, 216.4]
ds_speedup = ds_direct_times ./ ds_octree_times


p = plot(grid=:true, minorgrid=:true) 
plot!(p, ninteractions, octree_vs_direct, xscale=:log10)
xlabel!(p, "Interactions (\$N \\times N\$)")
ylabel!(p, "Speedup (Direct Time / Octree Time)")
title!(p, "Octree vs. Direct Solution, Speedup (\$\\phi = 0.1\$)")
display(p)
savefig(p,"dense_solenoid_speedup.svg")

p2 = plot(grid=:true, minorgrid=:true) 
plot!(p2, ninteractions, octree_times, scale=:log10, label="Octree")
plot!(p2, ninteractions, direct_times, scale=:log10, label="Direct")
xlabel!(p2, "Interactions (\$N \\times N\$)")
ylabel!(p2, "Execution Time [s]")
title!(p2, "Octree vs. Direct Solution Times (\$\\phi = 0.1\$)")
display(p2)
savefig(p2, "dense_solenoid_times.svg")