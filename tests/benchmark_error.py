from thor.testing import make_helmholtz

import thor 
import numpy as np 
import matplotlib.pyplot as plt 

# Benchmark error
nthetas = 10
theta_vals = np.linspace(0.1, 1.0, nthetas) 
errs = np.zeros(nthetas)
centroids, vol, jdensity = make_helmholtz(6.0)
bdirect = thor.bfield_direct(centroids, vol, jdensity, centroids)
for (i, theta) in enumerate(theta_vals): 
    boctree = thor.bfield_octree(centroids, vol, jdensity, centroids, theta=theta, leaf_threshold=1)
    bmag_direct = np.linalg.norm(bdirect, axis=1) 
    bmag_octree = np.linalg.norm(boctree, axis=1) 
    # errs[i] = thor.mean_relative_error(bmag_direct, bmag_octree)
    errs[i] = thor.testing.smape(bmag_direct, bmag_octree)

print(errs)
fig, ax = plt.subplots()
ax.plot(theta_vals, errs)
ax.set_xlabel("Theta (B-H Acceptance Criteria)")
ax.set_ylabel("B-Field Error Over Mesh")
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_title("Thor Benchmarks: Helmholtz Coil Problem\nAcceptance Criteria vs. Error")
fig.savefig("tests/benchmarks_error.png")
