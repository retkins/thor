cd(@__DIR__)

using Plots 

# phi = 0.1
# Test performed on 11/17/2025
# Test machine: AMD Ryzen 9 9950X, single-threaded execution
# Note that all of the octree times were recorded, but the last two direct times are estimated
nvals = [1_000, 10_000, 100_000, 1_000_000, 10_000_000]
ninteractions = nvals.^2
octree_times = [0.003600, 0.127102, 3.224143, 32.136709, 323.883919]
direct_times = [0.003215, 0.216577, 21.698355, 21.698355*100, 21.698355*10000]
octree_vs_direct = direct_times ./ octree_times

# Solution increments to support O(N) claim
octree_increment = zeros(4)
direct_increment = zeros(4)
for i in 2:5
    octree_increment[i-1] = octree_times[i]/octree_times[i-1]
    direct_increment[i-1] = direct_times[i]/direct_times[i-1]
end

p = plot(grid=:true, minorgrid=:true) 
plot!(p, nvals, octree_vs_direct, scale=:log10)
xlabel!(p, "Interactions (\$N \\times N\$)")
ylabel!(p, "Speedup (Direct Time / Octree Time)")
title!(p, "Octree vs. Direct Solution, Speedup (\$\\phi = 0.1\$)")
display(p)
savefig(p,"loop_speedup.svg")

p2 = plot(grid=:true, minorgrid=:true) 
plot!(p2, nvals, octree_times, scale=:log10, label="Octree")
plot!(p2, nvals, direct_times, scale=:log10, label="Direct")
xlabel!(p2, "Interactions (\$N \\times N\$)")
ylabel!(p2, "Execution Time [s]")
title!(p2, "Octree vs. Direct Solution Times (\$\\phi = 0.1\$)")
display(p2)
savefig(p2, "loop_times.svg")