from oersted.testing import make_helmholtz

import oersted
import numpy as np
import matplotlib.pyplot as plt

# Benchmark error
nthetas = 1
size = 15.0
theta_vals = np.linspace(0.5, 1.0, nthetas)
errs = np.zeros(nthetas)
centroids, vol, jdensity = make_helmholtz(size)
bdirect = oersted.bfield_direct(centroids, vol, jdensity, centroids)
for i, theta in enumerate(theta_vals):
    boctree = oersted.bfield_octree(centroids, vol, jdensity, centroids, theta=theta, leaf_threshold=1)
    bmag_direct = np.linalg.norm(bdirect, axis=1)
    bmag_octree = np.linalg.norm(boctree, axis=1)
    # errs[i] = oersted.mean_relative_error(bmag_direct, bmag_octree)
    errs[i] = oersted.testing.smape(bmag_direct, bmag_octree)

print(errs)
fig, ax = plt.subplots()
ax.plot(theta_vals, errs)
ax.set_xlabel("Theta (B-H Acceptance Criteria)")
ax.set_ylabel("B-Field Error Over Mesh")
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_title("oersted Benchmarks: Helmholtz Coil Problem\nAcceptance Criteria vs. Error")
fig.savefig("tests/fig/benchmarks_error.svg")
