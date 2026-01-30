""" Use the helmholtz coil problem as a benchmarking example
""" 

import thor 
import numpy as np 
import matplotlib.pyplot as plt 
from time import perf_counter

nbenches: int = 10
theta: float = 0.5

mesh_sizes = np.linspace(33.0, 2.0, nbenches)
direct_times = []
direct_interactions = []
est_direct_times = []
est_direct_interactions = []
i_est = 0
octree_times = np.zeros(nbenches)
interactions = np.zeros(nbenches)

for (i, mesh_size) in enumerate(mesh_sizes):

    centroids, vol, jdensity = thor.make_helmholtz(mesh_size)
    n = centroids.shape[0] 
    interactions[i] = n*n

    if n < 50000:
        start = perf_counter()
        bdirect = thor.bfield_direct(centroids, vol, jdensity, centroids)
        end = perf_counter() 
        direct_times.append(end - start)
        est_direct_times = [end - start ]
        direct_interactions.append(n*n)
        est_direct_interactions = [n*n]
        i_est = i
    else:
        m = (direct_times[i_est] - direct_times[0]) / (interactions[i_est] - interactions[0])
        b = direct_times[i_est] - m*interactions[i_est]
        est_direct_interactions.append(n*n)
        est_direct_times.append(m*n*n + b)
    

    start = perf_counter()
    boctree = thor.bfield_octree(centroids, vol, jdensity, centroids, theta=theta, leaf_threshold=1)
    end = perf_counter() 
    octree_times[i] = end - start 

# Plot solution times
fig, ax = plt.subplots()
ax.plot(direct_interactions, direct_times, 'r', label='direct')
ax.plot(interactions, octree_times, 'k', label='octree')
ax.plot(est_direct_interactions, est_direct_times, 'r--', label='direct trend')
ax.set_xlabel("Interactions ($N^2$)")
ax.set_ylabel("Solution time [s]")
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_title(f"Thor Benchmarks: Helmholtz Coil Problem\n$\\theta={theta:.2}$")
ax.legend()
fig.savefig("tests/benchmarks.png")

